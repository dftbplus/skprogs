"""Contains the representation of the old SK-file."""

import os.path

import numpy as np

from . import common as sc
import sktools.twocenter_grids



# Dummy null spline
NULL_SPLINE = """
Spline
12 0.0553585
112.9353346817185 2.801373701455403 -0.1119994835253462
0.035 0.0375    0.204206 -35.71077211012958 2016.504000000031 24177.93762071238
0.0375 0.04    0.12791 -25.17491577974109 2197.838532155373 -120889.6881035729
0.04 0.0425    0.07682029999999999 -16.45240477090621 1291.165871378576 -57585.58520643491
0.0425 0.045    0.0428593 -11.07630513663398 859.2739823303137 16659.22892930921
0.045 0.04533    0.0207993 -6.467574682557872 984.2181993001326 -2167173.572075024
0.04533 0.045334    0.0186943 -6.526006277016704 -1161.283637054166 353213222.4907721
0.045334 0.046259    0.0186682 -6.518342311433599 3077.275032831984 -1324559.571220061
0.046259 0.047184    0.0142234 -4.225362350069925 -598.3777773036936 561811.1110751317
0.047184 0.0493131    0.0102476 -3.890262342340788 960.6480559297889 -100763.5210502349
0.0493131 0.0503195    0.00534702 -1.169934109375229 317.0412179256228 -143026.9144497911
0.0503195 0.0513259    0.00434492 -0.9663840979460291 -114.7856421811885 10348.58893883691
0.0513259 0.0553585    0.00326664 -1.165980214261954 -83.5411824570522 -5782.515169399558 27636944.82683195 -3877959552.095367

This SPLINE is just a DUMMY-SPLINE!!!!!!!!!!!!!!!
"""

FLOAT_FORMSTR = " {:20.12E}"


class OldSKFile:

    def __init__(self, extended, dr, hamiltonian, overlap, onsites=None,
                 spinpolerror=None, hubbardus=None, occupations=None, mass=None,
                 splinerep=None, polyrep=None, extratag=None):
        self.extended = extended
        self.dr = dr
        self.nr = hamiltonian.shape[0]
        self.hamiltonian = hamiltonian
        self.overlap = overlap
        self.homo = onsites is not None
        self.onsites = onsites
        self.spinpolerror = spinpolerror
        self.hubbardus = hubbardus
        self.occupations = occupations
        self.mass = mass
        self.splinerep = splinerep
        self.polyrep = polyrep
        self.extratag = extratag


    @classmethod
    def fromfile(cls, fname, homo):
        fp = open(fname, "r")
        line = fp.readline()
        extended = line.startswith("@")
        nshell = 4 if extended else 3
        ninteg = 20 if extended else 10
        if extended:
            line = fp.readline()
        words = sc.split_fortran_fields(line)
        dr = float(words[0])
        nr = int(words[1])
        if homo:
            values = sc.convert_fortran_floats(fp.readline())
            onsites = np.array(values[0:nshell], dtype=float)
            spinpolerror = float(values[nshell])
            hubbardus = np.array(values[nshell+1:2*nshell+1], dtype=float)
            occupations = np.array(values[2*nshell+1:3*nshell+1], dtype=float)
        else:
            onsites = spinpolerror = hubbardus = occupations = None
        values = sc.convert_fortran_floats(fp.readline())
        if homo:
            mass = values[0]
        else:
            mass = None
        polyrep = np.array(values[1:10], dtype=float)
        hamiltonian = np.zeros(( nr, ninteg ), dtype=float)
        overlap = np.zeros(( nr, ninteg ), dtype=float)
        for iline in range(nr - 1):
            values = sc.convert_fortran_floats(fp.readline())
            hamiltonian[iline,0:ninteg] = values[0:ninteg]
            overlap[iline,0:ninteg] = values[ninteg:2*ninteg]
        # Currently, everything after SK table is treated as spline repulsive
        splinerep = fp.read()
        fp.close()
        return cls(extended, dr, hamiltonian, overlap, onsites, spinpolerror,
                   hubbardus, occupations, mass, splinerep, polyrep)


    def tofile(self, fname):
        fp = open(fname, "w")
        if self.extended:
            fp.write("@ Data set with f-electrons, for DFTB+ only\n")
        fp.write("{:f} {:d}\n".format(self.dr, self.nr))
        nshell = 4 if self.extended else 3
        ninteg = 20 if self.extended else 10
        shellfloats = FLOAT_FORMSTR * nshell
        if self.homo:
            fp.write(shellfloats.format(*self.onsites))
            fp.write(FLOAT_FORMSTR.format(self.spinpolerror))
            fp.write(shellfloats.format(*self.hubbardus))
            fp.write(shellfloats.format(*self.occupations))
            fp.write("\n")
        if self.homo:
            fp.write(FLOAT_FORMSTR.format(self.mass))
        else:
            fp.write(FLOAT_FORMSTR.format(0.0))
        if self.polyrep is not None:
            polyfloats = FLOAT_FORMSTR * 9
            fp.write(polyfloats.format(self.polyrep))
        else:
            fp.write(" 0.0" * 9)
        fp.write(" 0.0" * 10 + "\n")
        integralfloats = FLOAT_FORMSTR * ninteg
        for ir in range(self.nr):
            fp.write(integralfloats.format(*self.hamiltonian[ir,:]))
            fp.write(integralfloats.format(*self.overlap[ir,:]))
            fp.write("\n")
        if self.splinerep:
            fp.write("\n")
            fp.write(self.splinerep)
            fp.write("\n")
        if self.extratag:
            for xx in self.extratag:
                fp.write(xx)
                fp.write("\n")
        fp.close()



class OldSKFileSet:

    def __init__(self, grid, hamiltonian, overlap, basis1, basis2=None,
                 onsites=None, spinpolerror=None, hubbardus=None,
                 occupations=None, mass=None, dummy_repulsive=False,
                 extraTag=None):

        self._dr, self._nr0 = self._get_grid_parameters(grid)
        self._dummy_repulsive = dummy_repulsive
        self._hamiltonian = hamiltonian
        self._overlap = overlap
        self._basis1 = basis1
        self._homo = basis2 is None
        self._extraTag = extraTag

        if self._homo:
            self._basis2 = self._basis1
            self._onsites = self._get_basis_indexed_dict(basis1, onsites)
            self._hubbardus = self._get_basis_indexed_dict(basis1, hubbardus)
            self._occupations = self._get_basis_indexed_dict(basis1,
                                                             occupations)
            self._spinpolerror = spinpolerror
            self._mass = mass
        else:
            self._basis2 = basis2

        self._SK_for_shell1, self._shells_in_SK1 = self._split_basis(basis1)
        if not self._homo:
            self._SK_for_shell2, self._shells_in_SK2 = self._split_basis(basis2)
        else:
            self._SK_for_shell2 = self._SK_for_shell1
            self._shells_in_SK2 = self._shells_in_SK1
        self._integmap = self.get_integralmap(self._basis1, self._basis2)


    def tofile(self, workdir, elem1name, elem2name):
        skfiles = self._get_skfiles()
        skfilenames = self._write_skfiles(workdir, elem1name, elem2name,
                                          skfiles)
        return skfilenames


    def _get_skfiles(self):
        """Returns array of old SK file object pairs (A, B) representing the
        interaction between two atoms.
        """
        nsk1, nsk2 = len(self._shells_in_SK1), len(self._shells_in_SK2)
        # Loop over the number of SK-files necessary to represent element
        skfiles = []
        for isk1 in range(nsk1):
            skfiles1 = []
            for isk2 in range(nsk2):
                homoskfile = self._homo and isk1 == isk2
                skfile1 = self._get_skfile(isk1, isk2, homoskfile=homoskfile,
                                           reverse=False)
                if homoskfile:
                    skfile2 = None
                else:
                    skfile2 = self._get_skfile(isk1, isk2, homoskfile=False,
                                               reverse=True)
                skfiles1.append(( skfile1, skfile2 ))
            skfiles.append(skfiles1)
        return skfiles


    def _get_skfile(self, isk1, isk2, homoskfile, reverse):
        shells1 = self._shells_in_SK1[isk1]
        shells2 = self._shells_in_SK2[isk2]
        maxang1 = self._get_highest_angmom(shells1)
        maxang2 = self._get_highest_angmom(shells2)
        extended = (maxang1 == 3 or maxang2 == 3)
        maxang = 3 if extended else 2
        if extended:
            skintegmap = self._get_oldsk_integralmap(3)
            ninteg = 20
        else:
            skintegmap = self._get_oldsk_integralmap(2)
            ninteg = 10
        oldsk_ham = self._map_to_oldsk_integral_table(
            self._hamiltonian, shells1, shells2, skintegmap, ninteg, reverse)
        oldsk_over = self._map_to_oldsk_integral_table(
            self._overlap, shells1, shells2, skintegmap, ninteg, reverse)
        padding = np.zeros(( self._nr0 - 1, ninteg ), dtype=float)
        oldsk_ham = np.vstack(( padding, oldsk_ham ))
        oldsk_over = np.vstack(( padding, oldsk_over ))
        if self._dummy_repulsive:
            repulsive = NULL_SPLINE
        else:
            repulsive = None
        if homoskfile:
            onsites = self._map_to_oldsk_shell_values(self._onsites, shells1,
                                                      maxang)
            hubbus = self._map_to_oldsk_shell_values(self._hubbardus, shells1,
                                                     maxang)
            occupations = self._map_to_oldsk_shell_values(self._occupations,
                                                          shells1, maxang)
            skfile = OldSKFile(
                extended, self._dr, oldsk_ham, oldsk_over, onsites=onsites,
                spinpolerror=self._spinpolerror, hubbardus=hubbus,
                occupations=occupations, mass=self._mass, splinerep=repulsive,
                extratag=self._extraTag)
        else:
            skfile = OldSKFile(extended, self._dr, oldsk_ham, oldsk_over,
                               splinerep=repulsive, extratag=self._extraTag)
        return skfile


    @staticmethod
    def _write_skfiles(workdir, elem1, elem2, skfiles):
        elem1_capital = sc.capitalize_elem_name(elem1)
        elem2_capital = sc.capitalize_elem_name(elem2)
        form_elem1 = "{elem:s}:{ind:d}" if len(skfiles) > 1 else "{elem:s}"
        form_elem2 = "{elem:s}:{ind:d}" if len(skfiles[0]) > 1 else "{elem:s}"
        skfilenames = []
        for isk1, skfiles1 in enumerate(skfiles):
            elem1name = form_elem1.format(elem=elem1_capital, ind=isk1+1)
            for isk2, skfile12 in enumerate(skfiles1):
                elem2name = form_elem2.format(elem=elem2_capital, ind=isk2+1)
                skfile_ab, skfile_ba = skfile12
                fname = "{}-{}.skf".format(elem1name, elem2name)
                skfilenames.append(fname)
                skfile_ab.tofile(os.path.join(workdir, fname))
                if skfile_ba is not None:
                    fname = "{}-{}.skf".format(elem2name, elem1name)
                    skfilenames.append(fname)
                    skfile_ba.tofile(os.path.join(workdir, fname))
        return skfilenames


    @staticmethod
    def get_integralmap(basis1, basis2):
        """Gives column index for integral <n1,l1,m|n2,l2,m>."""
        integmap = {}
        ind = 0
        for n1, l1 in basis1:
            for n2, l2 in basis2:
                for mm in range(min(l1, l2) + 1):
                    integmap[n1, l1, n2, l2, mm] = ind
                    ind += 1
        return integmap


    @staticmethod
    def _get_oldsk_integralmap(lmax):
        skintegmap = {}
        ind = 0
        for l1 in range(lmax, -1, -1):
            for l2 in range(lmax, l1 - 1, -1):
                for mm in range(min(l1, l2) + 1):
                    skintegmap[l1, l2, mm] = ind
                    ind += 1
        return skintegmap


    @staticmethod
    def _get_highest_angmom(shells):
        maxang = 0
        for nn, ll in shells:
            maxang = max(maxang, ll)
        return maxang


    def _map_to_oldsk_integral_table(self, mytable, shells1, shells2,
                                     skintegmap, ninteg, reverse):
        oldsk_table = np.zeros(( mytable.shape[0], ninteg ))
        for n1, l1 in shells1:
            for n2, l2 in shells2:
                for mm in range(min(l1, l2) + 1):
                    if reverse:
                        ioldsk = skintegmap.get(( l2, l1, mm ), None)
                        prefac = float(1 - 2 * ((l1 + l2) % 2))
                    else:
                        ioldsk = skintegmap.get(( l1, l2, mm), None)
                        prefac = 1.0
                    if ioldsk is not None:
                        imy = self._integmap[n1, l1, n2, l2, mm]
                        oldsk_table[:,ioldsk] = prefac * mytable[:,imy]
        return oldsk_table

    @staticmethod
    def _map_to_oldsk_shell_values(myvalues, shells, maxang):
        oldsk_values = np.zeros(maxang + 1, dtype=float)
        for nn, ll in shells:
            oldsk_values[maxang - ll] = myvalues[nn, ll]
        return oldsk_values


    @staticmethod
    def _split_basis(basis):
        # Max angular momentum
        lmax = 0
        for nn, ll in basis:
            lmax = max(lmax, ll)

        # separete valence shells by angular momentum
        basis_per_l = [ None, ] * (lmax + 1)
        for nn, ll in basis:
            if basis_per_l[ll] is None:
                basis_per_l[ll] = [ nn, ]
            else:
                basis_per_l[ll].append(nn)

        # How many sk table compatible atoms can represent a given basis
        nskatom = 0
        for lbasis in basis_per_l:
            if lbasis is not None:
                # noinspection PyTypeChecker
                nskatom = max(nskatom, len(lbasis))

        lastshells = [ None, ] * (lmax + 1)
        # Gives the atom number for given shell (nn, ll)
        iskatom_shell = {}
        # Gives the shells of a given iskatom ii.
        shells_iskatom = []
        for iskatom in range(nskatom):
            shells = []
            hasshell = False
            for ll in range(lmax, -1, -1):
                lbasis = basis_per_l[ll]
                # noinspection PyTypeChecker
                if len(lbasis):
                    nn = lbasis.pop(0)
                    lastshells[ll] = nn
                    iskatom_shell[nn, ll] = iskatom
                    hasshell = True
                elif hasshell:
                    nn = lastshells[ll]
                else:
                    continue
                shells.insert(0, ( nn, ll ))
            shells_iskatom.append(shells)

        return iskatom_shell, shells_iskatom


    @staticmethod
    def _get_grid_parameters(grid):
        if not isinstance(grid, sktools.twocenter_grids.EquidistantGrid):
            raise sc.SkgenException(
                "Can not handle grid type " + grid.__class__.__name__)
        dr = grid.gridseparation
        nr0 = int(np.rint(grid.gridstart / grid.gridseparation))
        if np.abs(nr0 * dr - grid.gridstart) > 1e-12:
            msg = "Start distance incommensurable with grid separation"
            raise sc.SkgenException(msg)
        return dr, nr0


    @staticmethod
    def _get_basis_indexed_dict(basis, values):
        mydict = {shell: value for shell, value in zip(basis, values)}
        return mydict
