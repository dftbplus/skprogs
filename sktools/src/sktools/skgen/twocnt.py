import os
import glob
import logging

import numpy as np

import sktools.common as sc
import sktools.radial_grid as soc
from . import common as ssc
from .compression import run_denscomp, run_wavecomp


logger = logging.getLogger("skgen.twocnt")


def run_twocnt(skdefs, elem1, elem2, builddir, searchdirs, onecnt_binary,
               twocnt_binary):
    logger.info("Started for {}-{}".format(
        sc.capitalize_elem_name(elem1), sc.capitalize_elem_name(elem2)))
    hetero = (elem1.lower() != elem2.lower())
    prereq1 = _get_compression_prereq(elem1, skdefs, builddir, searchdirs,
                                      onecnt_binary)
    if hetero:
        prereq2 = _get_compression_prereq(elem2, skdefs, builddir, searchdirs,
                                          onecnt_binary)
    else:
        prereq2 = None
    calculator = SkgenTwocnt(builddir, searchdirs, twocnt_binary)
    calculator.set_input(skdefs, elem1, elem2, prereq1, prereq2)
    calculator.find_or_run_calculation()
    logger.info("Finished")
    return calculator


def _get_compression_prereq(elem, skdefs, builddir, searchdirs, onecnt_binary):
    logger.info("Creating compressed atom prerequisite for {}".format(
        sc.capitalize_elem_name(elem)))
    calc_dens = run_denscomp(skdefs, elem, builddir, searchdirs, onecnt_binary)
    dir_dens = calc_dens.get_result_directory()
    result_dens = calc_dens.get_result()
    calc_wave = run_wavecomp(skdefs, elem, builddir, searchdirs, onecnt_binary)
    dirs_wave = calc_wave.get_result_directories()
    result_wave = calc_wave.get_result()
    return SkgenTwocntCompressionPrereq(dir_dens, result_dens, dirs_wave,
                                        result_wave)


class SkgenTwocntCompressionPrereq:

    def __init__(self, dens_dir, dens_result, wave_dirs, wave_result):
        self.dens_dir = dens_dir
        self.dens_result = dens_result
        self.wave_dirs = wave_dirs
        self.wave_result = wave_result


class SkgenTwocnt:

    def __init__(self, builddir, searchdirs, twocnt_binary):
        self._builddir = builddir
        self._searchdirs = searchdirs
        self._twocnt_binary = twocnt_binary
        self._elem1 = None
        self._elem2 = None
        self._hetero = False
        self._skdefs = None
        self._input = None
        self._compression_prereqs = None
        self._twocenter_searchdirs = None
        self._resultdir = None


    def set_input(self, skdefs, elem1, elem2, comp_prereq1, comp_prereq2):
        elem1 = elem1.lower()
        elem2 = elem2.lower()
        self._elem1 = min(elem1, elem2)
        self._elem2 = max(elem1, elem2)
        _elements_reversed = (self._elem1 != elem1)
        self._hetero = (self._elem1 != self._elem2)
        self._skdefs = skdefs
        if _elements_reversed:
            self._compression_prereqs = ( comp_prereq2, comp_prereq1 )
        else:
            self._compression_prereqs = ( comp_prereq1, comp_prereq2 )
        self._input = SkgenTwocntInput(self._skdefs, self._elem1, self._elem2)
        self._twocenter_searchdirs = ssc.get_twocenter_searchdirs(
            self._searchdirs, self._elem1, self._elem2)
        self._resultdir = None


    def find_or_run_calculation(self):
        previous_calc_dirs = ssc.get_matching_subdirectories(
            self._twocenter_searchdirs, ssc.TWOCNT_WORKDIR_PREFIX)
        resultdirs = self._input.get_all_dirs_with_matching_signature(
            previous_calc_dirs)
        calculation_needed = True
        for resultdir in resultdirs:
            if self._check_prereq_dir_links(resultdir):
                calculation_needed = False
                logger.info("Matching twocnt calculation found "
                            + sc.log_path(resultdir))
                break
        if calculation_needed:
            resultdir = ssc.create_twocenter_workdir(
                self._builddir, ssc.TWOCNT_WORKDIR_PREFIX, self._elem1,
                self._elem2)
            logger.info("Doing twocnt calculation " + sc.log_path(resultdir))
            self._create_prereq_dir_links(resultdir,)
            calculation = SkgenTwocntCalculation(self._input,
                                                 self._compression_prereqs)
            calculation.run_and_convert_results(resultdir, self._twocnt_binary)
            self._input.store_signature(resultdir)
        self._resultdir = resultdir


    def get_result_directory(self):
        return self._resultdir


    def get_result(self):
        if self._resultdir is None:
            self.find_or_run_calculation()
        return SkgenTwocntResult(self._resultdir)


    def _check_prereq_dir_links(self, workdir):
        prereq1, prereq2 = self._compression_prereqs
        links_ok = self._check_prereq_dir_links_for_elem(1, workdir, prereq1)
        if links_ok and self._hetero:
            tmp = self._check_prereq_dir_links_for_elem(2, workdir, prereq2)
            links_ok = links_ok and tmp
        return links_ok


    def _check_prereq_dir_links_for_elem(self, ielem, workdir, prereq):
        existing_links = self._get_existing_dir_links_for_elem(ielem, workdir)
        links_to_create = self._get_prereq_dir_links_for_elem(ielem, prereq,
                                                              workdir)
        if len(existing_links) != len(links_to_create):
            return False
        for linkname, linktarget in links_to_create:
            if linkname not in existing_links:
                return False
            linkname = os.path.realpath(os.path.join(workdir, linkname))
            linktarget = os.path.realpath(os.path.join(workdir, linktarget))
            if not os.path.samefile(linkname, linktarget):
                return False
        return True


    @staticmethod
    def _get_prereq_dir_links_for_elem(ielem, prereq, workdir):
        dir_links = []
        densdir = os.path.relpath(prereq.dens_dir, workdir)
        linkname = "{}{:d}".format(ssc.DIRLINK_POTDENS_PREFIX, ielem)
        dir_links.append(( linkname, densdir ))
        for ind, wavedir in enumerate(prereq.wave_dirs):
            wavedir = os.path.relpath(wavedir, workdir)
            linkname = "{}{:d}.{:d}".format(ssc.DIRLINK_WAVE_PREFIX, ielem,
                                            ind + 1)
            dir_links.append(( linkname, wavedir ))
        return dir_links


    @staticmethod
    def _get_existing_dir_links_for_elem(ielem, workdir):
        glob1 = os.path.join(workdir,
                             "{}{:d}*".format(ssc.DIRLINK_POTDENS_PREFIX,
                                              ielem))
        dir_links = glob.glob(glob1)
        glob2 = os.path.join(workdir,
                             "{}{:d}*".format(ssc.DIRLINK_WAVE_PREFIX, ielem))
        dir_links += glob.glob(glob2)
        dir_links_basename = [ os.path.basename(mydir) for mydir in dir_links
                               if os.path.exists(mydir) ]
        return set(dir_links_basename)


    def _delete_existing_dir_links(self, workdir):
        links_to_delete = self._get_existing_dir_links_for_elem(1, workdir)
        links_to_delete.update(
            self._get_existing_dir_links_for_elem(2, workdir))
        for link in links_to_delete:
            os.remove(os.path.join(workdir, link))


    def _create_prereq_dir_links(self, workdir):
        prereq1, prereq2 = self._compression_prereqs
        links_to_create = self._get_prereq_dir_links_for_elem(1, prereq1,
                                                              workdir)
        if self._hetero:
            links_to_create += self._get_prereq_dir_links_for_elem(2, prereq2,
                                                                   workdir)
        for linkname, linktarget in links_to_create:
            os.symlink(linktarget, os.path.join(workdir, linkname))



class SkgenTwocntInput(ssc.InputWithSignature):

    SIGNATURE_FILE = ssc.TWOCNT_SIGNATURE_FILE

    def __init__(self, skdefs, elem1, elem2):
        atomparam1 = skdefs.atomparameters[elem1]
        atomparam2 = skdefs.atomparameters[elem2]
        self.atomconfig1 = atomparam1.atomconfig
        self.atomconfig2 = atomparam2.atomconfig
        twocentpars = skdefs.twocenterparameters[(elem1, elem2)]
        self.calculator = twocentpars.calculator
        self.grid = twocentpars.grid
        self.hetero = elem1 != elem2
        self.superposition = skdefs.globals.superposition
        self.functional = skdefs.globals.xcfunctional


    def get_signature(self):
        signature = {
            "atomconfig1": self.atomconfig1,
            "atomconfig2": self.atomconfig2,
            "calculator": self.calculator,
            "grid": self.grid,
            "hetero": self.hetero,
            "superposition": self.superposition,
            "functional": self.functional
        }
        return signature



class SkgenTwocntCalculation:

    def __init__(self, myinput, prerequisites):
        self._input = myinput
        self._prereq1, self._prereq2 = prerequisites
        self._twocnt_calculator = ssc.TwocenterCalculatorWrapper(
            myinput.calculator)


    def run_and_convert_results(self, workdir, twocnt_binary):
        atom1data = self._get_atomdata(self._input.atomconfig1, self._prereq1)
        if self._input.hetero:
            atom2data = self._get_atomdata(self._input.atomconfig2,
                                           self._prereq2)
        else:
            atom2data = None
        self._twocnt_calculator.do_calculation(
            self._input.superposition, self._input.functional, self._input.grid,
            atom1data, atom2data, twocnt_binary, workdir)
        result = self._twocnt_calculator.get_output(workdir)
        self._store_results(result, workdir)


    def _get_atomdata(self, atomconfig, atomcalcs):
        atomdata = sc.ClassDict()
        atomdata.potentials = atomcalcs.dens_result.get_potential()
        atomdata.density = atomcalcs.dens_result.get_density()
        atomdata.wavefuncs = self._get_standardized_compressed_wfcs(
            atomconfig, atomcalcs.wave_result)
        return atomdata


    @staticmethod
    def _get_standardized_compressed_wfcs(atomconfig, wavecomp_result):
        wavefuncs = []
        waves_found_for_shell = {}
        for nn, ll in atomconfig.valenceshells:
            wfc012 = wavecomp_result.get_wavefunction(nn, ll)
            wfc0_data = wfc012.data[:,0]
            wfc0_grid = wfc012.grid
            norm = wfc0_grid.dot(wfc0_data, wfc0_data)
            logger.debug("Norm for wavefunc {:d}{:s}: {:f}".format(
                nn, sc.ANGMOM_TO_SHELL[ll], norm))
            wfc0 = soc.GridData(wfc0_grid, wfc0_data)
            sign = get_normalized_sign(nn, ll, wfc0)
            logger.debug("Sign for wavefunc {:d}{:s}: {:.1f}".format(
                nn, sc.ANGMOM_TO_SHELL[ll], sign))
            wfc012.data *= sign
            previous_wfcs = waves_found_for_shell.get(ll, [])
            if len(previous_wfcs):
                coeffs = get_expansion_coefficients(wfc012, previous_wfcs)
                msg = "Expansion coeffs of previous wavefuncs:"
                msg += " {:f}" * len(coeffs)
                logger.debug(msg.format(*coeffs))
                wfc012 = orthogonalize_wave_and_derivatives(
                    wfc012, previous_wfcs, coeffs)
            newwave = ( nn, ll, wfc012 )
            if ll in waves_found_for_shell:
                waves_found_for_shell[ll].append(newwave)
            else:
                waves_found_for_shell[ll] = newwave
            wavefuncs.append(newwave)
        return wavefuncs


    @staticmethod
    def _store_results(result, workdir):
        result_file = os.path.join(workdir, ssc.TWOCNT_RESULT_FILE)
        sc.store_as_shelf(result_file, hamiltonian=result.get_hamiltonian(),
                          overlap=result.get_overlap())



class SkgenTwocntResult:

    def __init__(self, workdir):
        self._result_db = sc.retrive_from_shelf(
            os.path.join(workdir, ssc.TWOCNT_RESULT_FILE))


    def get_hamiltonian(self):
        return self._result_db["hamiltonian"]


    def get_overlap(self):
        return self._result_db["overlap"]



def get_normalized_sign(nn, ll, wavefunc):
    # Note: wavefunc data has shape (ngrid, 1)
    rR = wavefunc.grid.rr * wavefunc.data[:,0]
    imax = np.argmax(np.abs(rR))
    sign = np.sign(rR[imax])
    # Note: normalized sign should be n-independent to make sure also the
    # twocenter integration program can check if the conditions are fulfilled.
    normalized_sign = 1
    return float(sign / normalized_sign)


def get_expansion_coefficients(wavefunc, prev_wavefuncs):
    coeffs = []
    for wfcprev in prev_wavefuncs:
        if wavefunc.grid != wfcprev.grid:
            msg = "Incompatible grids found."
            raise sc.SkgenException(msg)
        coeffs.append(wavefunc.grid.dot(wavefunc.data[:,0], wfcprev.data[:,0]))
    return coeffs


def orthogonalize_wave_and_derivatives(wavefunc, prev_wavefuncs, coeffs):
    if len(prev_wavefuncs) == 0:
        return wavefunc
    wfcnew_data = wavefunc.data.copy()
    for coeff, wfcprev in zip(coeffs, prev_wavefuncs):
        wfcnew_data -= coeff * wfcprev.data
    norm = wavefunc.grid.dot(wfcnew_data[:,0], wfcnew_data[:,0])
    wfcnew_data /= norm
    return soc.GridData(wavefunc.grid, wfcnew_data)
