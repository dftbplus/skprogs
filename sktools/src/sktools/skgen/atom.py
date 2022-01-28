import os.path
import copy
import logging
import numpy as np
import sktools.common as sc
from . import common as ssc


logger = logging.getLogger("skgen.atom")


def run_atom(skdefs, elem, builddir, searchdirs, onecnt_binary,
             eigenonly=False, eigenspinonly=False):
    logger.info("Started for {}".format(
        sc.capitalize_elem_name(elem)))
    calculator = SkgenAtom(builddir, searchdirs, onecnt_binary)
    calculator.set_input(skdefs, elem)
    calculator.find_or_run_calculation(eigenonly, eigenspinonly)
    logger.info("Finished")
    return calculator


class SkgenAtom:

    def __init__(self, builddir, searchdirs, onecenter_binary):
        self._builddir = builddir
        self._searchdirs = searchdirs
        self._onecenter_binary = onecenter_binary
        self._elem = None
        self._input = None
        self._onecenter_searchdirs = None
        self._resultdir = None


    def set_input(self, skdefs, elem):
        elemlow = elem.lower()
        self._elem = elemlow
        self._input = SkgenAtomInput(skdefs, elemlow)
        self._onecenter_searchdirs = ssc.get_onecenter_searchdirs(
            self._searchdirs, self._elem)
        self._resultdir = None


    def find_or_run_calculation(self, eigenonly=False, eigenspinonly=False):
        previous_calc_dirs = ssc.get_matching_subdirectories(
            self._onecenter_searchdirs, ssc.ATOM_WORKDIR_PREFIX)
        resultdir = self._input.get_first_dir_with_matching_signature(
            previous_calc_dirs)
        if not resultdir:
            resultdir = self._do_calculation(eigenonly, eigenspinonly)
            if not (eigenonly or eigenspinonly):
                self._input.store_signature(resultdir)
        else:
            logger.info("Matching calculation found " + sc.log_path(resultdir))
        self._resultdir = resultdir


    def get_result(self):
        if self._resultdir is None:
            self.find_or_run_calculation()
        return SkgenAtomResult(self._resultdir)


    def get_result_directory(self):
        return self._resultdir


    def _do_calculation(self, eigenonly=False, eigenspinonly=False):
        workdir = ssc.create_onecenter_workdir(
            self._builddir, ssc.ATOM_WORKDIR_PREFIX, self._elem)
        calculation = SkgenAtomCalculation(self._input, workdir,
                                           self._onecenter_binary)
        calculation.run_and_convert_results(eigenonly, eigenspinonly)
        return workdir



class SkgenAtomInput(ssc.InputWithSignature):

    SIGNATURE_FILE = ssc.ATOM_SIGNATURE_FILE

    def __init__(self, skdefs, elem):
        self.elem = elem
        atomparams = skdefs.atomparameters[elem]
        self.atomconfig = atomparams.atomconfig
        self.xcf = skdefs.globals.xcf
        self.onecentpars = skdefs.onecenterparameters[elem]


    def get_signature(self):
        signature = {
            "atomconfig": self.atomconfig,
            "onecentpars": self.onecentpars,
            "xcf": self.xcf
        }
        return signature


class SkgenAtomCalculation:

    def __init__(self, myinput, workdir, binary):
        self._atomconfig = myinput.atomconfig
        self._delta_occ = myinput.onecentpars.deltafilling
        self._valence_shell_empty = self._get_valence_shell_empty(
            self._atomconfig, self._delta_occ)
        calculator = myinput.onecentpars.calculator
        self._oncenter_calculator = ssc.OnecenterCalculatorWrapper(calculator)
        self._xcf = myinput.xcf
        self._workdir = workdir
        self._binary = binary


    def run_and_convert_results(self, eigenonly, eigenspinonly):
        spin_needed = self._atomconfig.spinpolarized and not eigenonly
        result_spin_atom = None
        if eigenspinonly or spin_needed:
            result_spin_atom = self._calculate_spinpolarized_atom()
        result_spinavg_atom = None
        if eigenonly or not eigenspinonly:
            result_spinavg_atom = self._calculate_spinaveraged_atom()
        if eigenonly or eigenspinonly:
            return

        xcn = self._xcf.__class__.__name__.lower()
        if xcn in ('xclcpbe', 'xclcbnl'):
            hubbus = self._calculate_hubbus_corrected(result_spinavg_atom,
                                            replace_empty_with_homo=True)
        else:
            hubbus = self._calculate_hubbus(result_spinavg_atom,
                                            replace_empty_with_homo=True)

        self._log_substitutions(result_spinavg_atom)
        self._log_hubbus(hubbus)
        spinws = self._calculate_spinws(result_spinavg_atom,
                                        replace_empty_with_homo=True)
        self._log_spinws(spinws)
        self._convert_results(result_spinavg_atom, result_spin_atom, hubbus,
                              spinws)


    @staticmethod
    def _get_valence_shell_empty(atomconfig, delta_occ):
        valence_shell_empty = [
            atomconfig.occupations[ll][nn - ll - 1][0] < delta_occ
            for nn, ll in atomconfig.valenceshells]
        return valence_shell_empty


    def _calculate_spinpolarized_atom(self):
        workdir = os.path.join(self._workdir, "atom0_spin")
        logger.info("Calculating spin polarized atom " + sc.log_path(workdir))
        self._atomconfig.make_spinpolarized()
        result_spin = self._calculate_free_atom(workdir)
        return result_spin


    def _calculate_spinaveraged_atom(self):
        workdir = os.path.join(self._workdir, "atom0")
        logger.info("Calculating spin averaged atom " + sc.log_path(workdir))
        self._atomconfig.make_spinaveraged()
        result_spinavg = self._calculate_free_atom(workdir)
        return result_spinavg


    def _calculate_free_atom(self, workdir):
        output = self._oncenter_calculator.do_calculation(
            self._atomconfig, self._xcf, None, self._binary, workdir)
        result = self._collect_free_atom_result(output)
        self._log_free_atom_result(result)
        return result


    def _collect_free_atom_result(self, output):
        result = sc.ClassDict()
        result.etot = output.get_energy()
        eigvals0 = []
        occs0 = []
        for nn, ll in self._atomconfig.valenceshells:
            eigval = (output.get_eigenvalue(0, nn, ll),
                      output.get_eigenvalue(1, nn, ll))
            eigvals0.append(eigval)
            occ = (output.get_occupation(0, nn, ll),
                   output.get_occupation(1, nn, ll))
            occs0.append(occ)
        homo0 = output.get_homo_or_lowest_nl(0)
        homo1 = output.get_homo_or_lowest_nl(1)
        result.valence_eigvals = np.array(eigvals0, dtype=float)
        result.valence_occs = np.array(occs0, dtype=float)
        result.homo = (homo0, homo1)
        result.homo_eigval = (output.get_eigenvalue(0, homo0[0], homo0[1]),
                              output.get_eigenvalue(1, homo1[0], homo1[1]))
        return result


    def _calculate_hubbus(self, result_spinavg, replace_empty_with_homo):
        workdir = os.path.join(self._workdir, "hubbu")
        logger.info("Calculating Hubbard U values " + sc.log_path(workdir))
        shells, ihomo, energies = self._get_shells_and_energies_for_deriv_calc(
            result_spinavg, replace_empty_with_homo)
        sc.create_workdir(workdir)
        all_derivs = self._calc_deriv_matrix(workdir, shells, energies,
                                             spin_averaged=True)
        all_derivs = 0.5 * (all_derivs + np.transpose(all_derivs))
        valence_hubbus = self._get_valence_derivs(all_derivs, ihomo,
                                                  replace_empty_with_homo)
        return valence_hubbus


    def _calculate_hubbus_corrected(self, result_spinavg,
                                    replace_empty_with_homo):
        workdir = os.path.join(self._workdir, "hubbu")
        logger.info("Calculating Hubbard U values " + sc.log_path(workdir))
        shells, ihomo, energies = self._get_shells_and_energies_for_deriv_calc(
            result_spinavg, replace_empty_with_homo)
        sc.create_workdir(workdir)
        all_derivs = self._calc_deriv_matrix(workdir, shells, energies,
                                             spin_averaged=True)
        all_derivs = 0.5 * (all_derivs + np.transpose(all_derivs))
        valence_hubbus = self._get_valence_derivs(all_derivs, ihomo,
                                                  replace_empty_with_homo)

        logger.info("Applying the Hubbard correction algorithm")
        # needs: ll, omega, Hubbard
        # value of l: print(shells[ihomo][1])
        # valence Hubbard: print(valence_hubbus[ihomo][ihomo])
        # omega: print(self._xcf.omega)
        corr_hubb = self._correct_hubbard(shells[ihomo][1],
                                          valence_hubbus[ihomo][ihomo])
        # correction
        valence_hubbus[ihomo][ihomo] = corr_hubb

        return valence_hubbus


    def _correct_hubbard(self, ll, hubbu):
        '''Hubbard correction for range-separated hybrid functionals'''

        tau2hub = 3.2
        convtol = 1.0e-8
        maxiter = 50
        degeneracy = 0.5 / (2 * ll + 1)
        uu = hubbu * tau2hub / self._xcf.omega
        xxmin = uu
        xxmax = uu + tau2hub * degeneracy
        xx = 0.5 * (xxmin + xxmax)
        err = xxmax - xxmin
        iteration = 0
        while True:
            yymin = (xxmin * degeneracy * (1.0 - self.funcP(xxmin)))
            yymax = (xxmax * degeneracy * (1.0 - self.funcP(xxmax)))
            xxmin = uu + yymin
            xxmax = uu + yymax
            err = xxmax-xxmin
            xx = (xxmin + xxmax) * 0.5
            iteration = iteration + 1
            if (abs(err) < convtol) or iteration > maxiter:
                break
        tau = xx * self._xcf.omega

        return tau / tau2hub


    @staticmethod
    def funcP(xx):
        tmp = xx * xx * (xx * xx + 0.8 * xx + 0.2) / (xx + 1.0)**4
        return tmp


    def _get_shells_and_energies_for_deriv_calc(self, result_spinavg,
                                                replace_empty_with_homo):
        spin = 0
        homoshell = result_spinavg.homo[spin]
        if (homoshell not in self._atomconfig.valenceshells and
                replace_empty_with_homo):
            homoshell_n, homoshell_l = homoshell
            shells_to_calculate = [(homoshell_n, homoshell_l)]
            reference_energies = [result_spinavg.homo_eigval[spin]]
            ihomo = 0
        else:
            shells_to_calculate = []
            reference_energies = []
            ihomo = self._atomconfig.valenceshells.index(homoshell)
        shells_to_calculate += [(nn, ll)
                                for nn, ll in self._atomconfig.valenceshells]
        reference_energies += [eigval[spin]
                               for eigval in result_spinavg.valence_eigvals]
        return shells_to_calculate, ihomo, reference_energies


    def _calc_deriv_matrix(self, workdir, shells_to_calculate,
                           reference_energies, spin_averaged):
        ncalcshells = len(shells_to_calculate)
        tmp = np.zeros((ncalcshells, ncalcshells), dtype=float)
        if spin_averaged:
            deriv_matrix = tmp
        else:
            deriv_matrix = (tmp, np.array(tmp))
        for ishell, shell_to_variate in enumerate(shells_to_calculate):
            deriv = self._calc_de_shells_docc(
                workdir, shells_to_calculate, reference_energies,
                shell_to_variate, spin_averaged=spin_averaged)
            if spin_averaged:
                deriv_matrix[ishell] = deriv
            else:
                deriv_matrix[0][ishell] = deriv[0]
                deriv_matrix[1][ishell] = deriv[1]
        return deriv_matrix


    def _get_valence_derivs(self, all_hubbus, ihomo, replace_empty_with_homo):
        if not replace_empty_with_homo:
            return all_hubbus
        nvalshells = len(self._atomconfig.valenceshells)
        # noinspection PyNoneFunctionAssignment
        valence_hubbus = np.empty(( nvalshells, nvalshells ), dtype=float)
        hubbu_inds = [ihomo if self._valence_shell_empty[ii] else ii
                      for ii in range(nvalshells)]
        for ii, ii_hubbu in enumerate(hubbu_inds):
            for jj, jj_hubbu in enumerate(hubbu_inds):
                valence_hubbus[ii, jj] = all_hubbus[ii_hubbu, jj_hubbu]
        return valence_hubbus


    def _calculate_spinws(self, result_spinavg, replace_empty_with_homo):
        workdir = os.path.join(self._workdir, "spinw")
        logger.info("Calculating spinw values " + sc.log_path(workdir))
        shells, ihomo, energies = self._get_shells_and_energies_for_deriv_calc(
            result_spinavg, replace_empty_with_homo)
        sc.create_workdir(workdir)
        all_derivs_up, all_derivs_dn = self._calc_deriv_matrix(
            workdir, shells, energies, spin_averaged=False)
        spinws = 0.5 * (all_derivs_up - all_derivs_dn)
        spinws = 0.5 * (spinws + np.transpose(spinws))
        valence_spinws = self._get_valence_derivs(spinws, ihomo,
                                                  replace_empty_with_homo)
        return valence_spinws


    def _calc_de_shells_docc(self, workdir, derived_shells, reference_eigvals,
                             variated_shell, spin_averaged=False):
        atomconfig = self._atomconfig
        orig_occ = copy.deepcopy(atomconfig.occupations)
        nvar, lvar = variated_shell
        if spin_averaged:
            delta_occ = [self._delta_occ / 2.0, self._delta_occ / 2.0]
        else:
            delta_occ = [self._delta_occ, 0.0]

        # Decide, whether backwards, forward or central difference must be
        # used and set approriate deltas for occupation variation and
        # prefactors
        occ_varshell = orig_occ[lvar][nvar - lvar - 1][0]
        if occ_varshell < self._delta_occ:
            delta_occ_prefacs = [1.0, 2.0]
            finite_diff_coeffs_delta = np.array([2.0, -0.5])
            finite_diff_coeff0 = -1.5
        elif occ_varshell > 2 * lvar + 1 - self._delta_occ:
            delta_occ_prefacs = [-1.0, -2.0]
            finite_diff_coeffs_delta = np.array([-2.0, 0.5])
            finite_diff_coeff0 = 1.5
        else:
            delta_occ_prefacs = [-1.0, 1.0]
            finite_diff_coeffs_delta = np.array([-0.5, 0.5])
            finite_diff_coeff0 = 0.0
        finite_diff_coeffs_delta = finite_diff_coeffs_delta / self._delta_occ
        finite_diff_coeff0 = finite_diff_coeff0 / self._delta_occ

        # Calculate derivative via finite differences
        tmp = finite_diff_coeff0 * np.array(reference_eigvals)
        if spin_averaged:
            de_shells_docc = [ tmp, ]
        else:
            de_shells_docc = [ tmp, np.array(tmp) ]

        for ii in range(len(delta_occ_prefacs)):
            localname = "{:d}{:s}_{:d}".format(nvar, sc.ANGMOM_TO_SHELL[lvar],
                                               ii + 1)
            localworkdir = os.path.join(workdir, localname)
            occs = self._atomconfig.occupations[lvar][nvar - lvar - 1]
            new_occs = ( occs[0] + delta_occ_prefacs[ii] * delta_occ[0],
                         occs[1] + delta_occ_prefacs[ii] * delta_occ[1] )
            atomconfig.occupations[lvar][nvar - lvar - 1] = new_occs
            result = self._oncenter_calculator.do_calculation(
                atomconfig, self._xcf, None, self._binary, localworkdir)
            for ss in range(len(de_shells_docc)):
                e_shells = [ result.get_eigenvalue(ss, nn, ll)
                             for nn, ll in derived_shells]
                de_shells_docc[ss] += (finite_diff_coeffs_delta[ii]
                                       * np.array(e_shells, dtype=float))
            atomconfig.occupations = copy.deepcopy(orig_occ)

        if spin_averaged:
            return de_shells_docc[0]
        else:
            return de_shells_docc


    def _log_free_atom_result(self, result):
        logger.debug("Total energy: {:.5f}".format(result.etot))
        logger.debug("Eigenvalues of valence orbitals:")
        eigvals = result.valence_eigvals
        occs = result.valence_occs
        for ii, nl in enumerate(self._atomconfig.valenceshells):
            nn, ll = nl
            msg = "  {:d}{:s}: {:13.8f} ({:6.4f}) {:13.8f} ({:6.4f})".format(
                nn, sc.ANGMOM_TO_SHELL[ll], eigvals[ii][0], occs[ii][0],
                eigvals[ii][1], occs[ii][1])
            logger.debug(msg)


    def _log_substitutions(self, refcalc):
        nhomo, lhomo = refcalc.homo[0]
        if np.any(self._valence_shell_empty):
            logger.debug("Shell substitutions:")
            av = self._atomconfig.valenceshells
            out = ["{:d}{:s}".format(av[ii][0], sc.ANGMOM_TO_SHELL[av[ii][1]])
                   for ii in range(len(av)) if self._valence_shell_empty[ii]]
            out.append("<-- {:d}{:s}".format(nhomo, sc.ANGMOM_TO_SHELL[lhomo]))
            logger.debug(" ".join(out))
        else:
            logger.debug("Shell substitutions: None")


    @staticmethod
    def _log_hubbus(hubbus):
        logger.debug(str(hubbus))


    @staticmethod
    def _log_spinws(spinws):
        logger.debug(spinws)


    def _convert_results(self, res_spinavg, res_spin, hubbus, spinws):
        results = {
            "eigenvalues": res_spinavg.valence_eigvals[:, 0],
            "occupations": 2.0 * res_spinavg.valence_occs[:, 0],
            "homo": np.array(res_spinavg.homo, dtype=int),
        }
        if res_spin is not None:
            results["spinpol_energy"] = res_spin.etot - res_spinavg.etot
        else:
            results["spinpol_energy"] = 0.0
        results["hubbardu"] = hubbus
        results["spinw"] = spinws
        sc.store_as_shelf(os.path.join(self._workdir, ssc.ATOM_RESULT_FILE),
                          results)



class SkgenAtomResult:

    def __init__(self, workdir):
        self._workdir = workdir
        self._result_db = sc.retrive_from_shelf(
            os.path.join(workdir, ssc.ATOM_RESULT_FILE))


    def get_eigenvalues(self):
        return self._result_db["eigenvalues"]


    def get_occupations(self):
        return self._result_db["occupations"]


    def get_homo_nl(self):
        return self._result_db["homo"]


    def get_spinpolarization_energy(self):
        return self._result_db["spinpol_energy"]


    def get_hubbardus(self):
        return self._result_db["hubbardu"]


    def get_spinws(self):
        return self._result_db["spinw"]
