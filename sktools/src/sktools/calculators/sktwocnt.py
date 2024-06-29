import os
import shelve
import subprocess as subproc
import numpy as np
from sktools import hsd
import sktools.hsd.converter as conv
import sktools.common as sc
from sktools import twocenter_grids
from sktools import radial_grid


SUPPORTED_FUNCTIONALS = {'lda' : 1, 'pbe' : 2, 'blyp' : 3, 'lcy-pbe' : 4,
                         'lcy-bnl' : 5, 'pbe0' : 6, 'b3lyp' : 7,
                         'camy-b3lyp' : 8, 'camy-pbeh' : 9, 'tpss': 10,
                         'scan': 11, 'r2scan': 12, 'r4scan': 13, 'task': 14,
                         'task+cc': 15}

INPUT_FILE = "sktwocnt.in"
STDOUT_FILE = "output"
BASISFUNCTION_FILE = "basisfuncs.dbm"
DEFAULT_BINARY = "sktwocnt"


class SktwocntSettings(sc.ClassDict):
    """Specific settings for sktwocnt program.

    Attributes
    ----------
    integrationpoints : int, int
        Two integers representing the nr. of points for radial and angular
        integration.
    """

    def __init__(self, integrationpoints):
        super().__init__()
        self.integrationpoints = integrationpoints

    @classmethod
    def fromhsd(cls, node, query):
        """Generate the object from HSD tree"""
        integrationpoints, child = query.getvalue(
            node, "integrationpoints", conv.int1, returnchild=True)
        if len(integrationpoints) != 2:
            raise hsd.HSDInvalidTagValueException(
                "Two integration point parameters must be specified", child)
        return cls(integrationpoints)

    def __eq__(self, other):
        if not isinstance(other, SktwocntSettings):
            return False
        if self.integrationpoints != other.integrationpoints:
            return False
        return True


class Sktwocnt:

    def __init__(self, workdir):
        self._workdir = workdir

    def set_input(self, settings, superpos, functional, grid, atom1data,
                  atom2data=None):
        myinput = SktwocntInput(settings, superpos, functional, grid, atom1data,
                                atom2data)
        myinput.write(self._workdir)

    def run(self, binary=DEFAULT_BINARY):
        runner = SktwocntCalculation(binary, self._workdir)
        runner.run()

    def get_result(self):
        result = SktwocntResult(self._workdir)
        return result


class SktwocntInput:

    _INTERACTION_FROM_NTYPES = {
        1: "homo",
        2: "hetero",
    }

    _POTENTIAL_SUPERPOS = "potential"

    def __init__(self, settings, superpos, functional, grid, atom1data,
                 atom2data=None):
        self._settings = settings
        self._atom1data = atom1data
        self._hetero = atom2data is not None
        if self._hetero:
            self._atom2data = atom2data
        else:
            self._atom2data = self._atom1data
        self._check_superposition(superpos)
        self._densitysuperpos = (superpos == sc.SUPERPOSITION_DENSITY)
        self._check_functional(functional.type)
        self._functional = functional
        self._check_grid(grid)
        self._grid = grid

    @staticmethod
    def _check_superposition(superpos):
        if superpos not in \
        [sc.SUPERPOSITION_POTENTIAL, sc.SUPERPOSITION_DENSITY]:
            msg = "Sktwocnt: Invalid superposition type"
            sc.SkgenException(msg)

    @staticmethod
    def _check_functional(functional):
        if functional not in SUPPORTED_FUNCTIONALS:
            raise sc.SkgenException("Invalid functional type")

    @staticmethod
    def _check_grid(grid):
        if not isinstance(grid, twocenter_grids.EquidistantGrid):
            msg = "Sktwocnt can only handle equidistant grids"
            raise sc.SkgenException(msg)

    def write(self, workdir):
        atomfiles1 = self._store_atomdata(workdir, self._atom1data, 1)
        if self._hetero:
            atomfiles2 = self._store_atomdata(workdir, self._atom2data, 2)
        else:
            atomfiles2 = None
        self._store_twocnt_input(workdir, atomfiles1, atomfiles2)
        self._store_basisfunctions(workdir)

    def _store_atomdata(self, workdir, atomdata, iatom):
        atomfiles = sc.ClassDict()
        atomfiles.wavefuncs = self._store_wavefuncs(workdir, atomdata.wavefuncs,
                                                    iatom)
        atomfiles.potential = self._store_potentials(workdir,
                                                     atomdata.potentials, iatom)
        atomfiles.density = self._store_density(workdir, atomdata.density,
                                                iatom)
        xcn = self._functional.type
        if xcn in ('lcy-bnl', 'lcy-pbe', 'pbe0', 'b3lyp', 'camy-b3lyp',
                   'camy-pbeh'):
            atomfiles.dens_wavefuncs = self._store_dens_wavefuncs(
                workdir, atomdata.dens_wavefuncs, iatom)
        atomfiles.occshells = atomdata.occshells
        return atomfiles

    @staticmethod
    def _store_wavefuncs(workdir, wavefuncs, iatom):
        wavefuncfiles = []
        for nn, ll, wfc012 in wavefuncs:
            fname = "wave{:d}_{:d}{:s}.dat".format(iatom, nn,
                                                   sc.ANGMOM_TO_SHELL[ll])
            wfc012.tofile(os.path.join(workdir, fname))
            wavefuncfiles.append((nn, ll, fname))
        return wavefuncfiles

    @staticmethod
    def _store_dens_wavefuncs(workdir, wavefuncs, iatom):
        wavefuncfiles = []
        for nn, ll, wfc012 in wavefuncs:
            fname = "dens_wave{:d}_{:d}{:s}.dat".format(iatom, nn,
                                                        sc.ANGMOM_TO_SHELL[ll])
            wfc012.tofile(os.path.join(workdir, fname))
            wavefuncfiles.append((nn, ll, fname))
        return wavefuncfiles

    @staticmethod
    def _store_potentials(workdir, potentials, iatom):
        fname = "potentials{:d}.dat".format(iatom)
        # Vxc up and down should be equivalent, twocnt reads only one.
        newdata = potentials.data.take((radial_grid.VNUC, radial_grid.VHARTREE,
                                        radial_grid.VXCUP), axis=1)
        newgriddata = radial_grid.GridData(potentials.grid, newdata)
        newgriddata.tofile(os.path.join(workdir, fname))
        return fname

    @staticmethod
    def _store_density(workdir, density, iatom):
        fname = "density{:d}.dat".format(iatom)
        density.tofile(os.path.join(workdir, fname))
        return fname

    def _store_basisfunctions(self, workdir):
        config = shelve.open(
            os.path.join(workdir, BASISFUNCTION_FILE), "n")
        config["basis1"] = [(nn, ll) for nn, ll, wfc012
                            in self._atom1data.wavefuncs]
        config["basis2"] = [(nn, ll) for nn, ll, wfc012
                            in self._atom2data.wavefuncs]
        config.close()

    def _store_twocnt_input(self, workdir, atomfiles1, atomfiles2=None):
        fp = open(os.path.join(workdir, INPUT_FILE), "w")
        self._write_twocnt_header(fp)
        self._write_twocnt_gridinfo(fp)
        self._write_twocnt_integration_parameters(fp)
        self._write_twocnt_atom_block(fp, atomfiles1)
        if self._hetero:
            self._write_twocnt_atom_block(fp, atomfiles2)
        fp.close()

    def _write_twocnt_header(self, fp):
        if self._densitysuperpos:
            superposname = 'density'
        else:
            superposname = 'potential'

        xcfkey = self._functional.type
        ixc = SUPPORTED_FUNCTIONALS[xcfkey]
        fp.write('{} {} {}\n'.format('hetero' if self._hetero else 'homo',
                                     superposname, ixc))

    def _write_twocnt_gridinfo(self, fp):
        '''Writes integration grid info.'''

        # long-range corrected functionals
        if self._functional.type in ('lcy-bnl', 'lcy-pbe'):
            # hardcoded parameters for the Becke integration,
            # -> should probably be moved to skdef.hsd
            becke = '2000 194 11 1.0'
            fp.write("{:f}\n".format(self._functional.omega))
            fp.write("{:s}\n".format(becke))
        # B3LYP
        elif self._functional.type == 'b3lyp':
            # hardcoded parameters for the Becke integration,
            # -> should probably be moved to skdef.hsd
            becke = '2000 194 11 1.0'
            fp.write("{:s}\n".format(becke))
        # PBE0
        elif self._functional.type == 'pbe0':
            becke = '2000 194 11 1.0'
            fp.write("{:f}\n".format(self._functional.alpha))
            fp.write("{:s}\n".format(becke))
        # CAM functionals
        elif self._functional.type in ('camy-b3lyp', 'camy-pbeh'):
            becke = '2000 194 11 1.0'
            fp.write("{:f} {:f} {:f}\n".format(self._functional.omega,
                                               self._functional.alpha,
                                               self._functional.beta))
            fp.write("{:s}\n".format(becke))

        fp.write("{:f} {:f} {:e} {:f}\n".format(
            self._grid.gridstart, self._grid.gridseparation,
            self._grid.tolerance, self._grid.maxdistance))

    def _write_twocnt_integration_parameters(self, fp):
        fp.write("{:d} {:d}\n".format(*self._settings.integrationpoints))

    def _write_twocnt_atom_block(self, fp, atomfiles):
        if self._functional.type in ('lcy-bnl', 'lcy-pbe', 'pbe0', 'b3lyp',
                                     'camy-b3lyp', 'camy-pbeh'):
            fp.write("{:d} {:d}\n".format(len(atomfiles.wavefuncs),
                                          len(atomfiles.dens_wavefuncs)))
        else:
            fp.write("{:d}\n".format(len(atomfiles.wavefuncs)))

        for nn, ll, wavefuncfile in atomfiles.wavefuncs:
            fp.write("'{}' {:d}\n".format(wavefuncfile, ll))

        if self._functional.type in ('lcy-bnl', 'lcy-pbe', 'pbe0', 'b3lyp',
                                     'camy-b3lyp', 'camy-pbeh'):
            occdict = {}
            for xx in atomfiles.occshells:
                occdict[xx[0]] = xx[1]

            for nn, ll, dens_wavefuncfile in atomfiles.dens_wavefuncs:
                fp.write("'{}' {:d} {:f}\n"
                         .format(dens_wavefuncfile, ll, occdict[(nn, ll)]))

        fp.write("'{}'\n".format(atomfiles.potential))
        if self._densitysuperpos:
            fp.write("'{}'\n".format(atomfiles.density))
        else:
            fp.write("'{}'\n".format("nostart"))


class SktwocntCalculation:

    def __init__(self, binary, workdir):
        self._binary = binary
        self._workdir = workdir

    def run(self):
        fpin = open(os.path.join(self._workdir, INPUT_FILE), "r")
        fpout = open(os.path.join(self._workdir, STDOUT_FILE), "w")
        proc = subproc.Popen([self._binary], cwd=self._workdir,
                             stdin=fpin, stdout=fpout, stderr=subproc.STDOUT)
        proc.wait()
        fpin.close()
        fpout.close()


class SktwocntResult:

    def __init__(self, workdir):
        basis1, basis2 = self._read_basis(workdir)
        self._integmap = self._create_integral_mapping(basis1, basis2)
        ninteg = len(self._integmap)
        self._skham = self._read_sktable(
            os.path.join(workdir, "at1-at2.ham.dat"), ninteg)
        self._skover = self._read_sktable(
            os.path.join(workdir, "at1-at2.over.dat"), ninteg)

    @staticmethod
    def _read_basis(workdir):
        config = shelve.open(os.path.join(workdir, BASISFUNCTION_FILE), "r")
        basis1 = list(config["basis1"])
        basis2 = list(config["basis2"])
        config.close()
        return basis1, basis2

    @staticmethod
    def _create_integral_mapping(basis1, basis2):
        ninteg = 0
        integmap = {}
        for n1, l1 in basis1:
            for n2, l2 in basis2:
                for mm in range(min(l1, l2) + 1):
                    ninteg += 1
                    integmap[(n1, l1, n2, l2, mm)] = ninteg
        return integmap

    @staticmethod
    def _read_sktable(fname, ninteg):
        fp = open(fname, "r")
        nline = int(fp.readline())
        # noinspection PyNoneFunctionAssignment,PyTypeChecker
        tmp = np.fromfile(fp, dtype=float, count=ninteg * nline, sep=" ")
        tmp.shape = (nline, ninteg)
        return tmp

    def get_hamiltonian(self):
        return self._skham

    def get_overlap(self):
        return self._skover
