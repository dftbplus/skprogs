import logging
import numpy as np
import sktools.oldskfile
import sktools.common as sc
from .atom import run_atom
from .twocnt import run_twocnt

logger = logging.getLogger("skgen.sktable")


def run_sktable(skdefs, elem1, elem2, builddir, searchdirs, onecnt_binary,
                twocnt_binary, workdir, add_dummy_repulsive):
    logger.info("Started for {}-{}".format(
        sc.capitalize_elem_name(elem1), sc.capitalize_elem_name(elem2)))
    hetero = (elem1.lower() != elem2.lower())
    prereq_atom1 = _get_sktable_atom_prereq(elem1, skdefs, builddir, searchdirs,
                                            onecnt_binary)
    if hetero:
        prereq_atom2 = _get_sktable_atom_prereq(elem2, skdefs, builddir,
                                                searchdirs, onecnt_binary)
    else:
        prereq_atom2 = None
    prereq_twocnt = _get_sktable_twocnt_prereq(
        elem1, elem2, skdefs, builddir, searchdirs, onecnt_binary,
        twocnt_binary)
    calculator = SkgenSktable(builddir, searchdirs)
    calculator.set_input(skdefs, elem1, elem2, prereq_atom1, prereq_atom2,
                         prereq_twocnt)
    skfiles_written = calculator.write_sktables(workdir, add_dummy_repulsive)
    logger.info("Finished")
    return skfiles_written


def _get_sktable_atom_prereq(elem, skdefs, builddir, searchdirs, onecnt_binary):
    logger.info("Creating free atom prerequisite for {}".format(
        sc.capitalize_elem_name(elem)))
    calc_atom = run_atom(skdefs, elem, builddir, searchdirs, onecnt_binary)
    dir_atom = calc_atom.get_result_directory()
    result_atom = calc_atom.get_result()
    return SkgenSktableAtomPrereq(dir_atom, result_atom)


def _get_sktable_twocnt_prereq(elem1, elem2, skdefs, builddir, searchdirs,
                               onecnt_binary, twocnt_binary):
    logger.info("Creating twocnt prerequisite for {}-{}".format(
        sc.capitalize_elem_name(elem1), sc.capitalize_elem_name(elem2)))
    calc_twocnt = run_twocnt(skdefs, elem1, elem2, builddir, searchdirs,
                             onecnt_binary, twocnt_binary)
    dir_twocnt = calc_twocnt.get_result_directory()
    result_twocnt = calc_twocnt.get_result()
    return SkgenSktableTwocntPrereq(dir_twocnt, result_twocnt)


class SkgenSktableAtomPrereq:

    def __init__(self, directory, result):
        self.directory = directory
        self.result = result



class SkgenSktableTwocntPrereq:

    def __init__(self, directory, result):
        self.directory = directory
        self.result = result



class SkgenSktable:

    def __init__(self, builddir, searchdirs):
        self._builddir = builddir
        self._searchdirs = searchdirs
        self._skdefs = None
        self._elem1 = None
        self._elem2 = None
        self._input = None
        self._atom_prereqs = None
        self._twocnt_prereq = None


    def set_input(self, skdefs, elem1, elem2, atom_prereq1, atom_prereq2,
                  twocnt_prereq):
        elem1 = elem1.lower()
        elem2 = elem2.lower()
        self._elem1 = min(elem1, elem2)
        self._elem2 = max(elem1, elem2)
        _elements_reversed = (self._elem1 != elem1)
        self._skdefs = skdefs
        if _elements_reversed:
            self._atom_prereqs = ( atom_prereq2, atom_prereq1 )
        else:
            self._atom_prereqs = ( atom_prereq1, atom_prereq2 )
        self._twocnt_prereq = twocnt_prereq
        self._input = SkgenSktableInput(self._skdefs, self._elem1, self._elem2)


    def write_sktables(self, workdir, add_dummy_repulsive):
        assembly = SkgenSktableAssembly(self._input, self._atom_prereqs,
                                        self._twocnt_prereq)
        skfiles_written = assembly.write_sktables(workdir, add_dummy_repulsive)
        return skfiles_written



class SkgenSktableInput:

    def __init__(self, skdefs, elem1, elem2):
        self.elem1 = elem1
        self.elem2 = elem2
        self.homo = (elem1 == elem2)
        atomparam1 = skdefs.atomparameters[elem1]
        atomparam2 = skdefs.atomparameters[elem2]
        self.atomconfig1 = atomparam1.atomconfig
        self.atomconfig2 = atomparam2.atomconfig
        self.xcf = skdefs.globals.xcf
        if self.homo:
            dftbatom = atomparam1.dftbatom
            self.shellresolved = dftbatom.shellresolved
            self.custom_onsites = dftbatom.customizedonsites
            self.custom_hubbards = dftbatom.customizedhubbards
            self.custom_occupations = dftbatom.customizedoccupations
        else:
            self.shellresolved = None
            self.custom_onsites = None
            self.custom_hubbards = None
        twocntpars = skdefs.twocenterparameters[(elem1, elem2)]
        self.grid = twocntpars.grid



class SkgenSktableAssembly:

    def __init__(self, myinput, atom_prereqs, twocnt_prereq):
        self._input = myinput
        self._atom_prereq1, self._atom_prereq2 = atom_prereqs
        self._twocnt_prereq = twocnt_prereq


    def write_sktables(self, workdir, add_dummy_repulsive):
        result_twocnt = self._twocnt_prereq.result
        ham = result_twocnt.get_hamiltonian()
        over = result_twocnt.get_overlap()
        myinput = self._input
        valshells1 = myinput.atomconfig1.valenceshells
        valshells2 = myinput.atomconfig2.valenceshells
        grid = myinput.grid

        # if range-separated hybrid is used, add the RangeSep tag
        extra_tag = None
        xcn = myinput.xcf.type
        if xcn in ('lcy-bnl', 'lcy-pbe'):
            rsh_tag = "RangeSep\nLC {:f}".format(myinput.xcf.omega)
            extra_tag = [rsh_tag]
        if xcn in ('camy-b3lyp', 'camy-pbeh'):
            rsh_tag = "RangeSep\nCAM {:f} {:f} {:f}".format(myinput.xcf.omega,
                                                            myinput.xcf.alpha,
                                                            myinput.xcf.beta)
            extra_tag = [rsh_tag]
        if xcn in ('pbe0', 'b3lyp'):
            rsh_tag = "GlobalHybrid"
            extra_tag = [rsh_tag]

        if self._input.homo:
            onsites, occs, hubbus, spinpolerr, mass = self._get_atomic_data()
            if not myinput.shellresolved:
                hubbus = self._override_with_homo_value(
                    myinput.atomconfig1, self._atom_prereq1.result, hubbus)
            skfiles = sktools.oldskfile.OldSKFileSet(
                grid, ham, over, valshells1, None, onsites=onsites,
                spinpolerror=spinpolerr, hubbardus=hubbus, occupations=occs,
                mass=mass, dummy_repulsive=add_dummy_repulsive,
                extraTag=extra_tag)
        else:
            skfiles = sktools.oldskfile.OldSKFileSet(
                grid, ham, over, valshells1, valshells2,
                dummy_repulsive=add_dummy_repulsive, extraTag=extra_tag)
        files_written = skfiles.tofile(workdir, myinput.elem1, myinput.elem2)
        return files_written


    def _get_atomic_data(self):
        myinput = self._input
        shells = myinput.atomconfig1.valenceshells
        atomresult = self._atom_prereq1.result
        # Occupation can not be overriden by the users -> only defaults supplied
        occs = self._get_shell_value_or_default(
            shells, myinput.custom_occupations, atomresult.get_occupations())
        onsites = self._get_shell_value_or_default(
            shells, myinput.custom_onsites, atomresult.get_eigenvalues())
        # SkgenAtom returns Hubbard U matrix
        hubbus = atomresult.get_hubbardus()
        diag_hubbus = np.diagonal(hubbus)
        hubbus = self._get_shell_value_or_default(
            shells, myinput.custom_hubbards, diag_hubbus)
        spinpolerror = atomresult.get_spinpolarization_energy()
        mass = myinput.atomconfig1.mass
        return onsites, occs, hubbus, spinpolerror, mass


    @staticmethod
    def _get_shell_value_or_default(shells, shellvalues, defaults):
        result = []
        for ii in range(len(shells)):
            nn, ll = shells[ii]
            result.append(shellvalues.get(( nn, ll ), defaults[ii]))
        return result


    @staticmethod
    def _override_with_homo_value(atomconfig, atomresult, values):
        shells = atomconfig.valenceshells
        homo_nl_up, homo_nl_down = atomresult.get_homo_nl()
        if np.any(homo_nl_up != homo_nl_down):
            msg = "Different homo for spin up and down ({} vs. {})".format(
                sc.shell_ind_to_name(*homo_nl_up),
                sc.shell_ind_to_name(*homo_nl_down))
            raise sc.SkgenException(msg)
        # Homo indices may be stored in a numpy array
        homo_nl = tuple(homo_nl_up)
        try:
            ind = shells.index(homo_nl)
        except IndexError:
            msg = "Homo shell {} not among valence shells".format(
                sc.shell_ind_to_name(*homo_nl))
            raise sc.SkgenException(msg)
        return [ values[ind], ] * len(values)
