import os
import subprocess as subproc
import numpy as np
import sktools.hsd.converter as conv
import sktools.common as sc
from sktools.taggedfile import TaggedFile
import sktools.compressions
import sktools.radial_grid as oc


AVAILABLE_FUNCTIONALS = [ sc.XC_FUNCTIONAL_LDA, sc.XC_FUNCTIONAL_PBE ]
INPUT_FILE = "slateratom.in"
STDOUT_FILE = "output"
DEFAULT_BINARY = "slateratom"


def register_onecenter_calculator():
    """Returns data for calculator registration"""
    calc = sc.ClassDict()
    calc.settings = SlaterAtomSettings
    calc.calculator = SlaterAtom
    return calc


def register_hsd_settings():
    return SlaterAtomSettings


class SlaterAtomSettings(sc.ClassDict):
    """Specific settings for slateratom program.

    Attributes
    ----------
    exponents : list
        [ exp_s, exp_p, ... ] list, where each exp_* is a list of the exponents
        for the given angular momentum.
    maxpowers : list
        Maximal power for every angular momentum.
    """

    def __init__(self, exponents, maxpowers):
        super().__init__()
        self.exponents = exponents
        self.maxpowers = maxpowers

    @classmethod
    def fromhsd(cls, root, query):
        node = query.getchild(root, "exponents")
        exponents = sc.get_shellvalues_list(node, query, conv.float1)
        node = query.getchild(root, "maxpowers")
        maxpowers = sc.get_shellvalues_list(node, query, conv.int0)
        return cls(exponents, maxpowers)

    def __eq__(self, other):
        if not isinstance(other, SlaterAtomSettings):
            return False
        if len(self.exponents) != len(other.exponents):
            return False
        if len(self.maxpowers) != len(other.maxpowers):
            return False
        for ll in range(len(self.exponents)):
            if self.maxpowers[ll] != other.maxpowers[ll]:
                return False
            myexps = self.exponents[ll]
            otherexps = other.exponents[ll]
            if len(myexps) != len(otherexps):
                return False
            for ii in range(len(myexps)):
                if abs(myexps[ii] - otherexps[ii]) > sc.INPUT_FLOAT_TOLERANCE:
                    return False
        return True


class SlaterAtom:

    def __init__(self, workdir):
        self._workdir = workdir

    def set_input(self, settings, atomconfig, functional, compression):
        myinput = SlateratomInput(settings, atomconfig, functional, compression)
        myinput.write(self._workdir)

    def run(self, binary=DEFAULT_BINARY):
        runner = SlateratomCalculation(binary, self._workdir)
        runner.run()

    def get_result(self):
        return SlateratomResult(self._workdir)


class SlateratomInput:
    """Represents the input of the slateratom program.

    Parameters
    ----------
    atomconfig : AtomConfig
        Configuration of the atom to be calculated.
    functional : str
        DFT functional type ('lda' or 'pbe')
    compressions : list
        List of PowerCompression objects. Either empty (no compression applied)
        or has a compression object for every angular momentum of the atom.
    settings : SlaterAtom
        Further detailed settings of the program.
    """

    _XCFUNCTIONALS = { sc.XC_FUNCTIONAL_LDA: 2, sc.XC_FUNCTIONAL_PBE: 3 }

    _LOGICALSTRS = { True: ".true.", False: ".false." }

    _COMMENT = "#"


    def __init__(self, settings, atomconfig, functional, compressions):
        self._settings = settings
        self._atomconfig = atomconfig
        znuc = self._atomconfig.atomicnumber
        if abs(znuc - int(znuc)) > 1e-12:
            msg = "Slateratom: Only integer atomic numbers are allowed"
            raise sc.SkgenException(msg)
        if len(settings.exponents) != atomconfig.maxang + 1:
            msg = "Slateratom: Missing STO exponents for some shells"
            raise sc.SkgenException(msg)
        if len(settings.maxpowers) != atomconfig.maxang + 1:
            msg = "Slateratom: Missing STO max. powers for some shells"
            raise sc.SkgenException(msg)
        myxcfuncs = sc.XC_FUNCTIONAL_LDA, sc.XC_FUNCTIONAL_PBE
        if functional not in myxcfuncs:
            msg = "Invalid xc-functional type for slateratom"
            raise sc.SkgenException(msg)
        self._functional = self._XCFUNCTIONALS[functional]

        if compressions is None:
            compressions = []
        for comp in compressions:
            if not isinstance(comp, sktools.compressions.PowerCompression):
                msg = "Invalid compressiont type {} for slateratom".format(
                    comp.__class__.__name__)
                raise sc.SkgenException(msg)
        maxang = atomconfig.maxang
        ncompr = len(compressions)
        if ncompr and ncompr != maxang + 1:
            msg = "Invalid number of compressions" \
                "(expected {:d}, got {:d})".format(maxang + 1, ncompr)
            raise sc.SkgenException(msg)
        self._compressions = compressions
        myrelativistics = sc.RELATIVISTICS_NONE, sc.RELATIVISTICS_ZORA
        if atomconfig.relativistics not in myrelativistics:
            raise sc.SkgenException("Invalid relativistics type for slateratom")
        self._relativistic = atomconfig.relativistics == sc.RELATIVISTICS_ZORA


    def write(self, workdir):
        """Writes a valid input for the program.

        Parameters
        ----------
        workdir : str
            Existing working directory where the input should be written to.
        """
        maxang = self._atomconfig.maxang
        out = [
            "{:d} {:d} {:d} {:s} \t{:s} znuc maxang nscc relativistic".format(
                int(self._atomconfig.atomicnumber), maxang, 120,
                self._LOGICALSTRS[self._relativistic], self._COMMENT),

            "{:d}\t{:s} functional: 0=HF, 1=X-Alpha, 2=PW-LDA, 3=PBE ".format(
                self._functional, self._COMMENT)
        ]

        # Compressions
        if not len(self._compressions):
            out += [ "{:g} {:d} \t{:s} Compr. radius and power ({:s})".format(
                1e30, 0, self._COMMENT, sc.ANGMOM_TO_SHELL[ll])
                for ll in range(maxang + 1) ]
        else:
            out += [ "{:g} {:g} \t{:s} Compr. radius and power ({:s})".format(
                compr.radius, compr.power, self._COMMENT,
                sc.ANGMOM_TO_SHELL[ll])
                for ll, compr in enumerate(self._compressions) ]

        out += [ "{:d} \t{:s} nr. of occupied shells ({:s})".format(
            len(occ), self._COMMENT, sc.ANGMOM_TO_SHELL[ll])
            for ll, occ in enumerate(self._atomconfig.occupations) ]

        # STO powers and exponents
        exponents = self._settings.exponents
        maxpowers = self._settings.maxpowers
        out += [ "{:d} {:d} \t{:s} nr. of exponents, max. power ({:s})".format(
            len(exponents[ll]), maxpowers[ll], self._COMMENT,
            sc.ANGMOM_TO_SHELL[ll])
            for ll in range(maxang + 1) ]
        out.append("{:s} \t{:s} automatic exponent generation".format(
            self._LOGICALSTRS[False], self._COMMENT))
        for ll, skexp_ang in enumerate(exponents):
            for ii, skexp in enumerate(skexp_ang):
                out.append("{:10f} \t{:s} exponent {:d} ({:s})".format(
                    skexp, self._COMMENT, ii + 1, sc.ANGMOM_TO_SHELL[ll]))

        out.append("{:s} \t{:s} write eigenvectors".format(
            self._LOGICALSTRS[False], self._COMMENT))
        out.append("{} {:g} \t{:s} broyden mixer, mixing factor".format(
            self._LOGICALSTRS[True], 0.1, self._COMMENT))

        # Occupations
        for ll, occperl in enumerate(self._atomconfig.occupations):
            for ii, occ in enumerate(occperl):
                nn = ii + 1 + ll  # principal quantum number
                out.append("{:g} {:g} \t{:s} occupations ({:d}{:s})".format(
                    occ[0], occ[1], self._COMMENT, nn, sc.ANGMOM_TO_SHELL[ll]))

        # Valence shell range
        valenceqns = [[ sc.MAX_PRINCIPAL_QN, 0 ], ] * (maxang + 1)
        for nn, ll in self._atomconfig.valenceshells:
            valenceqns[ll][0] = min(valenceqns[ll][0], nn)
            valenceqns[ll][1] = max(valenceqns[ll][1], nn)
        for ll, vqns in enumerate(valenceqns):
            out.append("{:d} {:d} \t{:s} valence shells from to ({:s})".format(
                vqns[0], vqns[1], self._COMMENT, sc.ANGMOM_TO_SHELL[ll]))

        fp = open(os.path.join(workdir, INPUT_FILE), "w")
        fp.write("\n".join(out))
        fp.close()


class SlateratomCalculation:
    """Represents a program run.

    Parameters
    ----------
    binary : str
        Binary to use.
    workdir : str
        Working directory with valid input in it.
    """

    def __init__(self, binary, workdir):
        self._binary = binary
        self._workdir = workdir

    def run(self):
        """Run the code."""
        fpin = open(os.path.join(self._workdir, INPUT_FILE), "r")
        fpout = open(os.path.join(self._workdir, STDOUT_FILE), "w")
        proc = subproc.Popen([ self._binary ], cwd=self._workdir,
                             stdin=fpin, stdout=fpout, stderr=subproc.STDOUT)
        proc.wait()
        fpin.close()
        fpout.close()


class SlateratomResult:
    """Represents the output of a run.

    Parameters
    ----------
    workdir : str
        Working directory with the output of a run.
    """

    def __init__(self, workdir):
        self._workdir = workdir
        fp = open(os.path.join(self._workdir, "energies.tag"), "r")
        self._energiestag = TaggedFile.fromfile(fp, transpose=True)
        fp.close()

    def get_homo_or_lowest_nl(self, ss):
        """Returns homo. If spin channel has no electrons, lowest level.
        """
        tagname = "eigenlevels_dn" if ss else "eigenlevels_up"
        energies = self._energiestag[tagname]
        tagname = "occupations_dn" if ss else "occupations_up"
        occupations = self._energiestag[tagname].flat
        sorted_energy_inds = np.argsort(energies.flat)
        if np.all(occupations < 1e-8):
            # No electrons (e.g. spin down in H) -> lowest level as homo
            homo = sorted_energy_inds[0]
        else:
            for homo in sorted_energy_inds[::-1]:
                if occupations[homo] >= 1e-8:
                    break
            else:
                raise sc.SkgenException("Homo not found!")
        homo_ll = homo // energies.shape[1]
        mm = homo % energies.shape[1]
        homo_nn = mm + homo_ll + 1
        return homo_nn, homo_ll

    def get_eigenvalue(self, ss, nn, ll):
        """Returns an eigenvalue.

        Parameters
        ----------
        ss : int
            Spin channel (0, 1, ...).
        nn : int
            Principal quantum number (1, 2, ...).
        ll : int
            Angular momentum (0, 1, ...).

        Returns
        -------
        eigenvalue : float
           Required eigenvalue.
        """
        if ss:
            tagname = "eigenlevels_dn"
        else:
            tagname = "eigenlevels_up"
        return self._energiestag[tagname][ll, nn - ll - 1]

    def get_occupation(self, ss, nn, ll):
        """Returns an occupation.

        Parameters
        ----------
        ss : int
            Spin channel (0, 1, ...).
        nn : int
            Principal quantum number (1, 2, ...).
        ll : int
            Angular momentum (0, 1, ...).

        Returns
        -------
        occupation : float
            Required occupation number.
        """
        if ss:
            tagname = "occupations_dn"
        else:
            tagname = "occupations_up"
        return self._energiestag[tagname][ll, nn - ll - 1]

    def get_energy(self):
        """Returns the total energy.

        Returns
        -------
        energy: float
            Required total energy.
        """
        return self._energiestag["total_energy"]

    def get_potentials(self):
        """Returns various potential components of the atom

        Returns
        -------
        potentials : GridData
           Grid data with following potentials:
           nuclear, coulomb, xc-spinup, xc-spindown.
        """
        fp = open(os.path.join(self._workdir, "pot.dat"), "r")
        fp.readline()
        fp.readline()
        ngrid = int(fp.readline())
        # noinspection PyNoneFunctionAssignment,PyTypeChecker
        pots = np.fromfile(fp, dtype=float, count=ngrid * 6, sep=" ")
        fp.close()
        pots.shape = (ngrid, 6)
        grid = oc.RadialGrid(pots[:, 0], pots[:, 1])
        potentials = pots[:,2:6]
        return oc.GridData(grid, potentials)

    def get_density012(self):
        """Returns the radial density and its first and second derivatives.

        Returns
        -------
        density : GridData
           Grid data with the density and its first and second derivatives.
        """
        fp = open(os.path.join(self._workdir, "dens.dat"), "r")
        fp.readline()
        fp.readline()
        fp.readline()
        fp.readline()
        fp.readline()
        ngrid = int(fp.readline())
        # noinspection PyNoneFunctionAssignment,PyTypeChecker
        dens = np.fromfile(fp, dtype=float, count=ngrid * 7, sep=" ")
        fp.close()
        dens.shape = (ngrid, 7)
        grid = oc.RadialGrid(dens[:,0], dens[:,1])
        density = dens[:,2:5]
        return oc.GridData(grid, density)

    def get_wavefunction012(self, ss, nn, ll):
        """Returns radial wave function and its first and second derivatives.

        Returns
        -------
        density : GridData
           Grid data with the wavefunction and its first and second derivatives.
        """
        if ss == 0:
            formstr = "wave_{:02d}{:s}_up.dat"
        else:
            formstr = "wave_{:02d}{:s}_dn.dat"
        wavefile = formstr.format(nn, sc.ANGMOM_TO_SHELL[ll])
        wavefile = os.path.join(self._workdir, wavefile)
        if not os.path.exists(wavefile):
            raise sc.SkgenException("Missing wave function file " + wavefile)
        fp = open(wavefile, "r")
        fp.readline()
        fp.readline()
        ngrid = int(fp.readline())
        fp.readline()
        # noinspection PyNoneFunctionAssignment,PyTypeChecker
        wavefunc = np.fromfile(fp, dtype=float, count=5 * ngrid, sep=" ")
        wavefunc.shape = (ngrid, 5)
        grid = oc.RadialGrid(wavefunc[:, 0], wavefunc[:, 1])
        wfcs = wavefunc[:,2:5]
        return oc.GridData(grid, wfcs)
