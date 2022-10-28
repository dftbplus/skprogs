'''
Module to carry out wavefunction and density compressions.
'''


import os.path
import logging
import sktools.common as sc
from . import common as ssc


logger = logging.getLogger('skgen.compression')


def run_denscomp(skdefs, elem, builddir, searchdirs, onecnt_binary):

    logger.info('Started for {}'.format(sc.capitalize_elem_name(elem)))
    calculator = SkgenDenscomp(builddir, searchdirs, onecnt_binary)
    calculator.set_input(skdefs, elem)
    calculator.find_or_run_calculation()
    logger.info('Finished')

    return calculator


def run_wavecomp(skdefs, elem, builddir, searchdirs, onecnt_binary):

    logger.info('Started for {}'.format(sc.capitalize_elem_name(elem)))
    calculator = SkgenWavecomp(builddir, searchdirs, onecnt_binary)
    calculator.set_input(skdefs, elem)
    calculator.find_or_run_calculation()
    logger.info('Finished')

    return calculator


def search_wavecomp(skdefs, elem, builddir, searchdirs):

    logger.info('Started for {}'.format(sc.capitalize_elem_name(elem)))
    calculator = SkgenWavecomp(builddir, searchdirs, None)
    calculator.set_input(skdefs, elem)
    calculator.find_calculation()
    logger.info('Finished')

    return calculator


class SkgenDenscomp:

    def __init__(self, builddir, searchdirs, onecenter_binary):
        self._builddir = builddir
        self._searchdirs = searchdirs
        self._onecenter_binary = onecenter_binary
        self._elem = None
        self._input = None
        self._onecenter_searchdirs = None
        self._resultdir = None
        self._occshells = []


    def set_input(self, skdefs, elem):
        elemlow = elem.lower()
        self._elem = elemlow
        atomparams = skdefs.atomparameters[elemlow]
        atomconfig = atomparams.atomconfig
        compression = atomparams.dftbatom.densitycompression
        compressions = [compression,] * (atomconfig.maxang + 1)
        xcfunc = skdefs.globals.xcf
        self._occshells = atomconfig.occshells
        calculator = skdefs.onecenterparameters[elemlow].calculator
        self._input = AtomCompressionInput(elemlow, atomconfig, compressions,
                                           xcfunc, calculator)
        self._onecenter_searchdirs = ssc.get_onecenter_searchdirs(
            self._searchdirs, self._elem)
        self._resultdir = None


    def find_or_run_calculation(self):
        previous_calc_dirs = ssc.get_matching_subdirectories(
            self._onecenter_searchdirs, ssc.COMPRESSION_WORKDIR_PREFIX)
        resultdir = self._input.get_first_dir_with_matching_signature(
            previous_calc_dirs)
        recalculation_need = not resultdir
        if recalculation_need:
            resultdir = ssc.create_onecenter_workdir(
                self._builddir, ssc.COMPRESSION_WORKDIR_PREFIX, self._elem)
            logger.info('Calculating compressed atom ' + sc.log_path(resultdir))
            calculation = AtomCompressionCalculation(self._input)
            calculation.run(resultdir, self._onecenter_binary)
        else:
            logger.info('Matching calculation found ' + sc.log_path(resultdir))
        self._extract_results_if_not_present(self._input, self._occshells,
                                             resultdir)
        if recalculation_need:
            self._input.store_signature(resultdir)
        self._resultdir = resultdir


    @staticmethod
    def _extract_results_if_not_present(myinput, occshells, resultdir):
        resultshelf = os.path.join(resultdir, ssc.DENSCOMP_RESULT_FILE)
        if sc.shelf_exists(resultshelf):
            return
        calculator = AtomCompressionResult(myinput.calculator)
        output = calculator.get_output(resultdir)
        result = {
            'potentials': output.get_potentials(),
            'density': output.get_density012()
        }
        for qn, _ in occshells:
            nn = qn[0]
            ll = qn[1]
            # needs name as shelf allows only strings as keys
            shellname = sc.shell_ind_to_name(nn, ll)
            result[shellname] = output.get_wavefunction012(0, nn, ll)
        sc.store_as_shelf(resultshelf, result)


    def get_result(self):
        if self._resultdir is None:
            self.find_or_run_calculation()
        return SkgenDenscompResult(self._resultdir)


    def get_result_directory(self):
        return self._resultdir


class SkgenDenscompResult:

    def __init__(self, resultdir):
        resultshelf = os.path.join(resultdir, ssc.DENSCOMP_RESULT_FILE)
        self._results = sc.retrive_from_shelf(resultshelf)

    def get_potential(self):
        return self._results['potentials']

    def get_density(self):
        return self._results['density']

    def get_dens_wavefunction(self, nn, ll):
        shellname = sc.shell_ind_to_name(nn, ll)
        try:
            wfc = self._results[shellname]
        except KeyError:
            msg = 'Missing wavefunction {}'.format(shellname)
            raise sc.SkgenException(msg)
        return wfc


class SkgenWavecomp:

    def __init__(self, builddir, searchdirs, onecenter_binary):
        self._builddir = builddir
        self._searchdirs = searchdirs
        self._onecenter_binary = onecenter_binary
        self._elem = None
        self._shells_and_inputs = []
        self._onecenter_searchdirs = None
        self._resultdirs = None


    def set_input(self, skdefs, elem):
        elemlow = elem.lower()
        self._elem = elemlow
        atomparams = skdefs.atomparameters[elemlow]
        atomconfig = atomparams.atomconfig
        xcfunc = skdefs.globals.xcf
        calculator = skdefs.onecenterparameters[elemlow].calculator
        comprcontainer = atomparams.dftbatom.wavecompressions
        atomcompressions = comprcontainer.getatomcompressions(atomconfig)
        self._shells_and_inputs = []
        for compressions, shells in atomcompressions:
            myinput = AtomCompressionInput(elemlow, atomconfig, compressions,
                                           xcfunc, calculator)
            self._shells_and_inputs.append((shells, myinput))
        self._onecenter_searchdirs = ssc.get_onecenter_searchdirs(
            self._searchdirs, self._elem)
        self._resultdirs = None


    def find_or_run_calculation(self):
        resultdirs = []
        resultdir_for_nl = {}
        previous_calc_dirs = ssc.get_matching_subdirectories(
            self._onecenter_searchdirs, ssc.COMPRESSION_WORKDIR_PREFIX)
        for shells, myinput in self._shells_and_inputs:
            shellnames = [sc.shell_ind_to_name(nn, ll) for nn, ll in shells]
            logger.info('Processing compression for shell(s) {}'.format(
                ' '.join(shellnames)))
            resultdir = myinput.get_first_dir_with_matching_signature(
                previous_calc_dirs)
            recalculation_needed = not resultdir
            if recalculation_needed:
                resultdir = ssc.create_onecenter_workdir(
                    self._builddir, ssc.COMPRESSION_WORKDIR_PREFIX, self._elem)
                logger.info(
                    'Calculating compressed atom ' + sc.log_path(resultdir))
                calculation = AtomCompressionCalculation(myinput)
                calculation.run(resultdir, self._onecenter_binary)
            else:
                logger.info(
                    'Matching calculation found ' + sc.log_path(resultdir))
            self._extract_results_if_not_present(myinput, shells, resultdir)
            if recalculation_needed:
                myinput.store_signature(resultdir)
            resultdirs.append(resultdir)
            for nn, ll in shells:
                resultdir_for_nl[(nn, ll)] = resultdir

        self._resultdirs = resultdirs
        self._resultdir_for_nl = resultdir_for_nl


    def find_calculation(self):
        resultdirs = []
        resultdir_for_nl = {}
        previous_calc_dirs = ssc.get_matching_subdirectories(
            self._onecenter_searchdirs, ssc.COMPRESSION_WORKDIR_PREFIX)
        for shells, myinput in self._shells_and_inputs:
            shellnames = [sc.shell_ind_to_name(nn, ll) for nn, ll in shells]
            logger.info('Processing compression for shell(s) {}'.format(
                ' '.join(shellnames)))
            resultdir = myinput.get_first_dir_with_matching_signature(
                previous_calc_dirs)
            recalculation_needed = not resultdir
            if recalculation_needed:
                sc.fatalerror('Could not find wavecomp calculation')
            logger.info(
                'Matching calculation found ' + sc.log_path(resultdir))
            resultdirs.append(resultdir)
            for nn, ll in shells:
                resultdir_for_nl[(nn, ll)] = resultdir

        self._resultdirs = resultdirs
        self._resultdir_for_nl = resultdir_for_nl


    @staticmethod
    def _extract_results_if_not_present(myinput, shells, resultdir):
        resultshelf = os.path.join(resultdir, ssc.WAVECOMP_RESULT_FILE)
        if sc.shelf_exists(resultshelf):
            return
        calculator = AtomCompressionResult(myinput.calculator)
        output = calculator.get_output(resultdir)
        resultdict = {}
        for nn, ll in shells:
            # Needs name as shelf allows only strings as keys
            shellname = sc.shell_ind_to_name(nn, ll)
            resultdict[shellname] = output.get_wavefunction012(0, nn, ll)
        sc.store_as_shelf(resultshelf, resultdict)


    def get_result(self):
        if self._resultdirs is None:
            self.find_or_run_calculation()
        return SkgenWavecompResult(self._resultdirs)


    def get_result_directories(self):
        return self._resultdirs


    def get_result_directory_for_shell(self, nn, ll):
        resdir = self._resultdir_for_nl.get((nn, ll), None)
        if resdir is None:
            msg = 'No result directory for shell {:s}'.format(
                sc.shell_ind_to_name(nn, ll))
            raise sc.SkgenException(msg)
        return resdir


class SkgenWavecompResult:

    def __init__(self, workdirs):
        self._result = {}
        for workdir in workdirs:
            resultshelf = os.path.join(workdir, ssc.WAVECOMP_RESULT_FILE)
            curres = sc.retrive_from_shelf(resultshelf)
            self._result.update(curres)


    def get_wavefunction(self, nn, ll):
        shellname = sc.shell_ind_to_name(nn, ll)
        try:
            wfc = self._result[shellname]
        except KeyError:
            msg = 'Missing wavefunction {}'.format(shellname)
            raise sc.SkgenException(msg)
        return wfc


class AtomCompressionInput(ssc.InputWithSignature):

    SIGNATURE_FILE = ssc.COMPRESSION_SIGNATURE_FILE

    def __init__(self, elem, atomconfig, shell_compressions, xcfunc,
                 calculator):
        self.elem = elem
        self.atomconfig = atomconfig
        self.shell_compressions = shell_compressions
        self.xcfunc = xcfunc
        self.calculator = calculator


    def get_signature(self):
        signature = {
            'atomconfig': self.atomconfig,
            'compressions': self.shell_compressions,
            'xcfunc': self.xcfunc,
            'calculator': self.calculator
        }
        return signature


class AtomCompressionCalculation:

    def __init__(self, myinput):
        self._atomconfig = myinput.atomconfig
        self._atomconfig.make_spinaveraged()

        self._shell_compressions = myinput.shell_compressions
        calculator = myinput.calculator
        self._onecnt_calculator = ssc.OnecenterCalculatorWrapper(calculator)
        self._xcfunc = myinput.xcfunc


    def run(self, workdir, binary):
        self._onecnt_calculator.do_calculation(
            self._atomconfig, self._xcfunc, self._shell_compressions, binary,
            workdir)


class AtomCompressionResult:

    def __init__(self, calculator):
        self._onecnt_calculator = ssc.OnecenterCalculatorWrapper(calculator)


    def get_output(self, workdir):
        return self._onecnt_calculator.get_output(workdir)
