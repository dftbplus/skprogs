import os
import glob
import logging
from .. import common as sc
from .. import calculators


logger = logging.getLogger("skgen.common")

SHELL_FORMAT = "{:d}{:s}"

ATOM_WORKDIR_PREFIX = "atom."
ATOM_SIGNATURE_FILE = "_atom-inp.db"
ATOM_RESULT_FILE = "_atom-res.db"

COMPRESSION_WORKDIR_PREFIX = "comp."
COMPRESSION_SIGNATURE_FILE = "_comp-inp.db"
DENSCOMP_RESULT_FILE = "_denscomp-res.db"
WAVECOMP_RESULT_FILE = "_wavecomp-res.db"

TWOCNT_WORKDIR_PREFIX = "twocnt."
TWOCNT_SIGNATURE_FILE = "_twocnt_inp.db"
TWOCNT_RESULT_FILE = "_twocnt-res.db"
DIRLINK_POTDENS_PREFIX = "dir_potdens"
DIRLINK_WAVE_PREFIX = "dir_wave"


class OnecenterCalculatorWrapper:

    def __init__(self, calcsettings):
        self._calculatorclass = get_calculator_class(
            calcsettings, calculators.ONECENTER_CALCULATORS)
        self._calculator_name = calcsettings.__class__.__name__
        self._calcsettings = calcsettings


    def do_calculation(self, atomconfig, xcfunc, compressions, binary, workdir):
        sc.create_workdir(workdir, reuse_existing=True)

        calculator = self._calculatorclass(workdir)
        calculator.set_input(self._calcsettings, atomconfig, xcfunc,
                             compressions)

        logger.debug("Running {}".format(binary))
        calculator.run(binary)
        result = calculator.get_result()
        return result


    def get_output(self, workdir):
        calculator = self._calculatorclass(workdir)
        result = calculator.get_result()
        return result



class TwocenterCalculatorWrapper:

    def __init__(self, calcsettings):
        self._calculatorclass = get_calculator_class(
            calcsettings, calculators.TWOCENTER_CALCULATORS)
        self._calculator_name = calcsettings.__class__.__name__
        self._calcsettings = calcsettings


    def do_calculation(self, superpos, functional, grid, atom1data, atom2data,
                       binary, workdir):
        sc.create_workdir(workdir, reuse_existing=True)
        calculator = self._calculatorclass(workdir)
        calculator.set_input(self._calcsettings, superpos, functional, grid,
                             atom1data, atom2data)

        logger.debug("Running {}".format(binary))
        calculator.run(binary)
        result = calculator.get_result()
        return result


    def get_output(self, workdir):
        calculator = self._calculatorclass(workdir)
        result = calculator.get_result()
        return result



class InputWithSignature:

    SIGNATURE_FILE = None

    def store_signature(self, workdir):
        sc.store_as_shelf(os.path.join(workdir, self.SIGNATURE_FILE),
                          self.get_signature())


    def get_first_dir_with_matching_signature(self, search_dirs):
        return sc.find_dir_with_matching_shelf(
            search_dirs, self.SIGNATURE_FILE, **self.get_signature())


    def get_all_dirs_with_matching_signature(self, search_dirs):
        return sc.get_dirs_with_matching_shelf(
            search_dirs, self.SIGNATURE_FILE, **self.get_signature())


    def get_signature(self):
        raise NotImplementedError



def get_matching_subdirectories(dirs, subdirprefix):
    dirglobs = [ os.path.join(mydir, subdirprefix + "*")
                 for mydir in dirs ]
    matching_subdirs = []
    for dirglob in dirglobs:
        matching_subdirs += glob.glob(dirglob)
    return matching_subdirs


def get_onecenter_searchdirs(searchdirs, elem):
    onecenter_searchdirs = [ os.path.join(dirname, get_onecenter_dirname(elem))
                             for dirname in searchdirs ]
    return onecenter_searchdirs


def get_twocenter_searchdirs(searchdirs, elem1, elem2):
    twocenter_searchdirs = [ os.path.join(dirname,
                                          get_twocenter_dirname(elem1, elem2))
                             for dirname in searchdirs ]
    return twocenter_searchdirs


def get_onecenter_dirname(elem):
    return elem


def get_twocenter_dirname(elem1, elem2):
    return "{}-{}".format(elem1, elem2)


def create_onecenter_workdir(builddir, workdir_prefix, elem):
    workroot = os.path.join(builddir, get_onecenter_dirname(elem))
    workdir = _create_workdir(workroot, workdir_prefix)
    return workdir


def create_twocenter_workdir(builddir, workdir_prefix, elem1, elem2):
    workroot = os.path.join(builddir, get_twocenter_dirname(elem1, elem2))
    workdir = _create_workdir(workroot, workdir_prefix)
    return workdir


def _create_workdir(workroot, workdir_prefix):
    sc.create_workdir(workroot, reuse_existing=True)
    workdir = sc.create_unique_workdir(workroot, workdir_prefix)
    return workdir


def get_calculator_class(settings, registered_calculators):
    for curr in registered_calculators:
        if isinstance(settings, curr.settings):
            return curr.calculator
    raise sc.SkgenException("Unknown calculator {}".format(
        settings.__class__.__name__))
