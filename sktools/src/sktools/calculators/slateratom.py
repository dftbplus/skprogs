'''
Module to perform or find atomic DFT calculations, using slateratom.
'''

import os
import sys
import logging
import subprocess as subproc
import numpy as np
import sktools.hsd.converter as conv
import sktools.common as sc
from sktools.taggedfile import TaggedFile
import sktools.compressions as skcomp
import sktools.radial_grid as oc
import sktools.xcfunctionals as xc


LOGGER = logging.getLogger('slateratom')

SUPPORTED_FUNCTIONALS = {'lda' : 2, 'pbe' : 3, 'blyp' : 4, 'lcy-pbe' : 5,
                         'lcy-bnl' : 6, 'pbe0' : 7, 'b3lyp' : 8,
                         'camy-b3lyp' : 9, 'camy-pbeh' : 10}

INPUT_FILE = "slateratom.in"
STDOUT_FILE = "output"
DEFAULT_BINARY = "slateratom"


def fatalerror(msg, errorcode=-1):
    '''Issue error message and exit.

    Args:

        msg (str): error message
        errorcode (int): error code to raise

    '''

    LOGGER.critical(msg)
    sys.exit(errorcode)


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

    def __init__(self, exponents, maxpowers, scftol, maxscfiter):
        super().__init__()
        self.exponents = exponents
        self.maxpowers = maxpowers
        self.scftol = scftol
        self.maxscfiter = maxscfiter

    @classmethod
    def fromhsd(cls, root, query):
        node = query.getchild(root, "exponents")
        exponents = sc.get_shellvalues_list(node, query, conv.float1)
        node = query.getchild(root, "maxpowers")
        maxpowers = sc.get_shellvalues_list(node, query, conv.int0)
        scftol = query.getvalue(
            root, "scftolerance", converter=conv.float0, defvalue=1.0e-10)
        maxscfiter = query.getvalue(
            root, "maxscfiterations", converter=conv.int0, defvalue=120)
        return cls(exponents, maxpowers, scftol, maxscfiter)

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
        if (abs(self.scftol - other.scftol)
                > sc.INPUT_FLOAT_TOLERANCE):
            return False
        if self.maxscfiter != other.maxscfiter:
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
        DFT functional type ('lda', 'pbe', 'blyp', 'lcpbe', 'lcbnl')
    compressions : list
        List of Compression objects. Either empty (no compression applied)
        or has a compression object for every angular momentum of the atom.
    settings : SlaterAtom
        Further detailed settings of the program.
    """

    _LOGICALSTRS = {True: ".true.", False: ".false."}

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

        if self._settings.scftol <= 0.0:
            msg = "Slateratom: SCF tolerance must be >0.0 a.u."
            raise sc.SkgenException(msg)

        if self._settings.maxscfiter < 1:
            msg = "Slateratom: Maximum number of SCF iterations must be >=1"
            raise sc.SkgenException(msg)

        if self.isXCFunctionalSupported(functional):
            xcfkey = functional.type
            self._functional = SUPPORTED_FUNCTIONALS[xcfkey]

            if xcfkey in ('lcy-pbe', 'lcy-bnl', 'camy-b3lyp', 'camy-pbeh'):
                self._omega = functional.omega
            else:
                self._omega = None

            if xcfkey in ('camy-b3lyp', 'camy-pbeh'):
                self._alpha = functional.alpha
                self._beta = functional.beta
            elif xcfkey == 'pbe0':
                self._alpha = functional.alpha
                self._beta = None
            else:
                self._alpha = None
                self._beta = None

        else:
            msg = 'Invalid xc-functional type for slateratom'
            raise sc.SkgenException(msg)

        if compressions is None:
            compressions = []
        for comp in compressions:
            if not isinstance(comp, (skcomp.PowerCompression,
                                     skcomp.WoodsSaxonCompression)):
                msg = "Invalid compression type {} for slateratom".format(
                    comp.__class__.__name__)
                raise sc.SkgenException(msg)
        maxang = atomconfig.maxang
        ncompr = len(compressions)
        if ncompr > 0 and ncompr != maxang + 1:
            msg = "Invalid number of compressions" \
                "(expected {:d}, got {:d})".format(maxang + 1, ncompr)
            raise sc.SkgenException(msg)

        # block different compression types of different shells for now
        if ncompr > 0:
            compids = [comp.compid for comp in compressions]
            if not len(set(compids)) == 1:
                msg = "At the moment, shells may only be compressed by the " \
                    + "same type of potential."
                raise sc.SkgenException(msg)

        self._compressions = compressions
        myrelativistics = sc.RELATIVISTICS_NONE, sc.RELATIVISTICS_ZORA
        if atomconfig.relativistics not in myrelativistics:
            raise sc.SkgenException("Invalid relativistics type for slateratom")
        self._relativistic = atomconfig.relativistics == sc.RELATIVISTICS_ZORA


    def isXCFunctionalSupported(self, functional):
        '''Checks if the given xc-functional is supported by the calculator,
           in particular: checks if AVAILABLE_FUNCTIONALS intersect with
           xc.XCFUNCTIONALS

        Args:
            functional: xc-functional, defined in xcfunctionals.py

        Returns:
            true, if xc-functional is supported, otherwise false
        '''

        tmp = []
        for xx in SUPPORTED_FUNCTIONALS:
            if xx in xc.XCFUNCTIONALS:
                tmp.append(xc.XCFUNCTIONALS[xx])

        return bool(functional.__class__ in tmp)


    def write(self, workdir):
        """Writes a valid input for the program.

        Parameters
        ----------
        workdir : str
            Existing working directory where the input should be written to.
        """
        maxang = self._atomconfig.maxang
        out = [
            "{:d} {:d} {:d} {:g} {:s} \t{:s}".format(
                int(self._atomconfig.atomicnumber), maxang,
                self._settings.maxscfiter, self._settings.scftol,
                self._LOGICALSTRS[self._relativistic], self._COMMENT) + \
            " znuc maxang nscc scftol relativistic",

            "{:d} \t\t\t{:s} functional".format(
                self._functional, self._COMMENT)
        ]

        # range-separated functionals
        xctype = list(SUPPORTED_FUNCTIONALS.keys())[
            list(SUPPORTED_FUNCTIONALS.values()).index(self._functional)]
        if xctype in ('lcy-pbe', 'lcy-bnl'):
            out += [
                "{:g} \t{:s} range-separation parameter (omega)".format(
                    self._omega, self._COMMENT),

                # numerical interator
                # hardcoded parameters for the Becke integration
                # --> should be moved to skdef.hsd!
                "2000 194 11 1.0 \t{:s} Becke integrator settings"
                .format(self._COMMENT)]
        # B3LYP
        elif xctype == 'b3lyp':
            out += [
                # numerical interator
                # hardcoded parameters for the Becke integration
                # --> should be moved to skdef.hsd!
                "2000 194 11 1.0 \t{:s} Becke integrator settings"
                .format(self._COMMENT)]
        # PBE0
        elif xctype == 'pbe0':
            out += [
                "{:g} \t{:s} ".format(self._alpha, self._COMMENT) + \
                "Global portion of HFX",

                # numerical interator
                # hardcoded parameters for the Becke integration
                # --> should be moved to skdef.hsd!
                "2000 194 11 1.0 \t{:s} Becke integrator settings"
                .format(self._COMMENT)]
        # CAM functionals
        elif xctype in ('camy-b3lyp', 'camy-pbeh'):
            out += [
                "{:g} {:g} {:g} \t{:s} ".format(
                    self._omega, self._alpha, self._beta, self._COMMENT) + \
                "range-separation parameter (omega), CAM alpha, CAM beta",

                # numerical interator
                # hardcoded parameters for the Becke integration
                # --> should be moved to skdef.hsd!
                "2000 194 11 1.0 \t{:s} Becke integrator settings"
                .format(self._COMMENT)]

        # Compressions
        if len(self._compressions) == 0:
            # no compression for all shells
            out += ["{:d} \t\t\t{:s} Compression ID".format(
                skcomp.SUPPORTED_COMPRESSIONS['nocompression'],
                self._COMMENT) + " ({:s})".format(sc.ANGMOM_TO_SHELL[ll])
                    for ll in range(maxang + 1)]
        else:
            # define the type of compression for each shell
            for ll, compr in enumerate(self._compressions):
                if compr.compid == \
                   skcomp.SUPPORTED_COMPRESSIONS['powercompression']:
                    out += ["{:d} \t\t{:s} Compression ID".format(
                        skcomp.SUPPORTED_COMPRESSIONS['powercompression'],
                        self._COMMENT) \
                            + " ({:s})".format(sc.ANGMOM_TO_SHELL[ll])]
                elif compr.compid == \
                     skcomp.SUPPORTED_COMPRESSIONS['woodssaxoncompression']:
                    out += ["{:d} \t\t{:s} Compression ID".format(
                        skcomp.SUPPORTED_COMPRESSIONS['woodssaxoncompression'],
                        self._COMMENT) \
                            + " ({:s})".format(sc.ANGMOM_TO_SHELL[ll])]
                else:
                    msg = 'Invalid compression type.'
                    raise sc.SkgenException(msg)
            # provide the compression parametrization for each shell
            for ll, compr in enumerate(self._compressions):
                if compr.compid == \
                   skcomp.SUPPORTED_COMPRESSIONS['powercompression']:
                    out += ["{:g} {:g} \t\t{:s} Compr. radius and power ({:s})"
                            .format(compr.radius, compr.power, self._COMMENT,
                                    sc.ANGMOM_TO_SHELL[ll])]
                elif compr.compid == \
                     skcomp.SUPPORTED_COMPRESSIONS['woodssaxoncompression']:
                    out += ["{:g} {:g} {:g} \t\t{:s}".format(
                        compr.onset, compr.cutoff, compr.vmax, self._COMMENT) \
                            + " Compr. onset, cutoff and vmax ({:s})"
                            .format(sc.ANGMOM_TO_SHELL[ll])]

        out += ["{:d} \t\t\t{:s} nr. of occupied shells ({:s})".format(
            len(occ), self._COMMENT, sc.ANGMOM_TO_SHELL[ll])
                for ll, occ in enumerate(self._atomconfig.occupations)]

        # STO powers and exponents
        exponents = self._settings.exponents
        maxpowers = self._settings.maxpowers
        out += ["{:d} {:d} \t\t\t{:s} nr. of exponents, max. power ({:s})"
                .format(len(exponents[ll]), maxpowers[ll], self._COMMENT,
                        sc.ANGMOM_TO_SHELL[ll])
                for ll in range(maxang + 1)]
        out.append("{:s} \t\t{:s} automatic exponent generation".format(
            self._LOGICALSTRS[False], self._COMMENT))
        for ll, skexp_ang in enumerate(exponents):
            for ii, skexp in enumerate(skexp_ang):
                out.append("{:10f} \t\t{:s} exponent {:d} ({:s})".format(
                    skexp, self._COMMENT, ii + 1, sc.ANGMOM_TO_SHELL[ll]))

        out.append("{:s} \t\t{:s} write eigenvectors".format(
            self._LOGICALSTRS[False], self._COMMENT))
        out.append("{} {:g} \t\t\t{:s} broyden mixer, mixing factor".format(
            2, 0.1, self._COMMENT))

        # Occupations
        for ll, occperl in enumerate(self._atomconfig.occupations):
            for ii, occ in enumerate(occperl):
                nn = ii + 1 + ll  # principal quantum number
                out.append("{:g} {:g} \t\t\t{:s} occupations ({:d}{:s})".format(
                    occ[0], occ[1], self._COMMENT, nn, sc.ANGMOM_TO_SHELL[ll]))

        # Occupied shell range
        occqns = [[sc.MAX_PRINCIPAL_QN + 1, 0],] * (maxang + 1)
        for qn, occ in self._atomconfig.occshells:
            nn = qn[0]
            ll = qn[1]
            occqns[ll][0] = min(occqns[ll][0], nn)
            occqns[ll][1] = max(occqns[ll][1], nn)
        for ll, vqns in enumerate(occqns):
            out.append("{:d} {:d} \t\t\t{:s} occupied shells from to ({:s})".format(
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
        self._check_output()
        fp = open(os.path.join(self._workdir, "energies.tag"), "r",
                  encoding="utf8")
        self._energiestag = TaggedFile.fromfile(fp, transpose=True)
        fp.close()

    def _check_output(self):
        """Checks calculation for SCF convergence and energies.tag file."""

        error_str = "SCF is NOT converged, maximal SCF iterations exceeded."
        troubleshoot_str = "Possible troubleshooting steps:" \
            + "\n" \
            + "1) Check skdef.hsd (e.g. an extremely large basis can result" \
            + " in unstable SCF cycles, ...)" \
            + "\n" \
            + "2) Increase MaxSCFIterations" \
            + "\n" \
            + "3) Choose less tight SCFTolerance"
        output_fname = os.path.join(self._workdir, "output")

        # Check for SCF convergence
        with open(output_fname, "r", encoding="utf8") as outfile:
            if error_str in outfile.read():
                fatalerror(error_str + "\n" + "Path: " + self._workdir
                           + "\n\n" + troubleshoot_str)

        # Check for any other reason, causing a missing energies.tag file
        error_str = "energies.tag absent, reason unknown (check yourself)."
        if not os.path.exists(os.path.join(self._workdir, "energies.tag")):
            fatalerror(error_str + "\n" + "Path: " + self._workdir)

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
