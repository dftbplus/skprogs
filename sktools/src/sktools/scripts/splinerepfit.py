#!/usr/bin/env python3

'''
Simple tool for fitting repulsives on data.
'''


import argparse
import numpy as np
import scipy.interpolate as spinter
from scipy.interpolate import CubicSpline
import scipy.optimize as spopt
from scipy import linalg

DESCRIPTION = '''Fits a spline repulsive to reference data. The script expects
and provides data in atomic units.'''

EPSILON = 1e-12


def get_cmdline_parser():
    """Set up the command-line argument parser for the splinerepfit tool.

    Returns:
        argparse.ArgumentParser: The configured argument parser with various
            options (rdamp, min-pow, max-pow, tail-length, etc.) and
            positional arguments (refdata, dftbdata, rmin, rcut).
    """

    parser = argparse.ArgumentParser(description=DESCRIPTION)

    msg = 'beginning of the damped region (damping goes from 1.0 to 0.0 from ' \
          'here to the cutoff)'
    parser.add_argument('--rdamp', type=float, default=None, help=msg)

    msg = 'minimal power for the polynomial fit (default: 3)'
    parser.add_argument('--min-pow', type=int, default=3, dest='minpow',
                        metavar='POW', help=msg)

    msg = 'maximal power for the polynomial fit (default: 11)'
    parser.add_argument('--max-pow', type=int, default=11, dest='maxpow',
                        metavar='POW', help=msg)

    msg = 'length of the repulsive tail fitted by a 5th order polynomial ' \
          '(default: 0.2)'
    parser.add_argument('--tail-length', type=float, default=0.2, metavar='LEN',
                        dest='taillen', help=msg)

    msg = 'grid separation for the splines (default: 0.02)'
    parser.add_argument('--spline-grid', type=float, default=0.02,
                        dest='splgrid', metavar='DR', help=msg)

    msg = 'grid separation for the internal interpolation (default: 0.001)'
    parser.add_argument('--interpol-grid', type=float, default=0.001,
                        dest='intergrid', metavar='DR', help=msg)

    msg = 'file with reference total energies'
    parser.add_argument('refdata', type=str, help=msg)

    msg = 'file with DFTB electronic energies'
    parser.add_argument('dftbdata', type=str, help=msg)

    msg = 'minimal distance to be considered in the fit'
    parser.add_argument('rmin', type=float, help=msg)

    msg = 'cutoff for the repulsive'
    parser.add_argument('rcut', type=float, help=msg)

    return parser


def get_data(fname):
    """Load x- and y-data columns from a text file.

    The file is expected to have at least two columns:
    the first column is taken as x-values, and the second as y-values.

    Args:
        fname (str): The path or filename from which to load data.

    Returns:
        tuple:
            - xx (numpy.ndarray): 1D array of x-values.
            - yy (numpy.ndarray): 1D array of y-values.
    """

    data = np.loadtxt(fname)
    xx = data[:, 0]
    yy = data[:, 1]

    return xx, yy


def polynomial(cc, minpow, rr, rcut):
    """
    Evaluate a polynomial of the form

        sum_{n=0}^{len(cc)-1} [ c_n * (rcut - r)^(n + minpow) ],

    at the distance(s) `rr`.

    This is done by first evaluating the polynomial

        P(x) = c_0 + c_1*x + c_2*x^2 + ... + c_{N-1}*x^{N-1}

    at x = (rcut - r), then multiplying by (rcut - r)^minpow.

    Args:
        cc (1D array): Polynomial coefficients, [c_0, c_1, ..., c_{N-1}].
        minpow (int) : Minimal power offset for the polynomial exponents.
        rr (array-like): Interatomic distance(s) at which to evaluate.
        rcut (float): The cutoff radius used in the polynomial expression.

    Returns:
        array-like or float:
            The polynomial evaluated at each of the distances in `rr`.
            If `rr` is a single value, a float is returned.
    """

    dr = rcut - rr
    res = cc[-1]

    for ii in range(len(cc) - 2, -1, -1):
        res = res * dr + cc[ii]

    for _ in range(minpow):
        res *= dr

    return res


def get_polyderiv(nderiv, cc, minpow):
    """
    Compute the n-th derivative of a polynomial of the form:

        sum_{n=0}^{len(cc)-1} [ c_n * (rcut - r)^(n + minpow) ].

    This function adjusts both the coefficients `cc` and the minimum power
    `minpow` to represent the polynomial after taking `nderiv` derivatives
    with respect to `r`.

    Args:
        nderiv (int): Order of the derivative to compute (e.g. 1 for first derivative).
        cc (array_like): 1D array of polynomial coefficients `[c_0, c_1, ..., c_{N-1}]`.
        minpow (int): The minimum power associated with `(rcut - r)^(n + minpow)`.

    Returns:
        tuple:
            - numpy.ndarray: The updated polynomial coefficients after taking
              `nderiv` derivatives.
            - int: The new minimum power corresponding to the derivative polynomial.
    """

    for _ in range(1, nderiv + 1):
        cc = np.array([-cc[ii] * (ii + minpow) for ii in range(len(cc))])
        minpow -= 1
        if minpow == -1:
            cc = cc[1:]
            minpow = 0

    return cc, minpow


def damping_cos(rdamp, rcut, rr):
    """
    Compute a cosine-based damping factor from `rdamp` to `rcut`.

    The damping factor is defined piecewise:
      - 1.0 for `rr < rdamp`
      - cos( (rr - rdamp) / (rcut - rdamp) * π/2 ) for `rr >= rdamp`

    Args:
        rdamp (float): The distance at which damping starts.
        rcut (float): The distance at which the damping is fully applied (cutoff).
        rr (float or array_like): Distance(s) at which to evaluate the damping.

    Returns:
        float or ndarray: The damping factor(s). If `rr` is scalar, returns float;
        otherwise returns an array of the same shape as `rr`.
    """

    return np.where(rr >= rdamp, np.cos((rr - rdamp) / (rcut - rdamp) * np.pi / 2.0), 1.0)


def get_spline_coeffs(xx, yy, deriv0, deriv1):
    """
    Compute cubic spline coefficients with specified first derivatives at both ends.

    Args:
        xx (array_like): 1D array of x-coordinates (must be strictly increasing).
        yy (array_like): 1D array of y-values, same length as `xx`.
        deriv0 (float): First derivative at the left endpoint (x = xx[0]).
        deriv1 (float): First derivative at the right endpoint (x = xx[-1]).

    Returns:
        CubicSpline: A SciPy CubicSpline object encapsulating the piecewise polynomial.
        numpy.ndarray: A transposed matrix of shape (len(xx) - 1, 4) containing
            the spline coefficients for each interval. Each row is
            `[c0, c1, c2, c3]` for the interval `[xx[i], xx[i+1]]`.
    """

    cs = CubicSpline(xx, yy, bc_type=((1, deriv0), (1, deriv1)))

    c0 = cs.c[3, :]
    c1 = cs.c[2, :]
    c2 = cs.c[1, :]
    c3 = cs.c[0, :]

    mtx = np.array([c0, c1, c2, c3])

    return cs, np.transpose(mtx)


def get_splineval012(splcoeffs, r0, rr):
    """
    Evaluate a cubic polynomial and its first two derivatives at shifted coordinates.

    The polynomial is assumed to have the form:
        f(dr) = alpha + beta*dr + gamma*dr^2 + delta*dr^3
    where dr = rr - r0.

    Args:
        splcoeffs (array_like): Length-4 array `[alpha, beta, gamma, delta]`.
        r0 (float): Reference point for the shift.
        rr (float or array_like): Value(s) at which to evaluate the polynomial.

    Returns:
        tuple:
            - f0 (float or ndarray): The function value(s) f(dr).
            - f1 (float or ndarray): The first derivative f'(dr).
            - f2 (float or ndarray): The second derivative f''(dr).
    """

    alpha, beta, gamma, delta = splcoeffs
    dr = rr - r0
    f0 = ((delta * dr + gamma) * dr + beta) * dr + alpha
    f1 = (3.0 * delta * dr + 2.0 * gamma) * dr + beta
    f2 = 6.0 * delta * dr + 2.0 * gamma

    return f0, f1, f2


def get_poly5coeffs(derivs, r0, rcut):
    """
    Solve for coefficients of a 5th-degree polynomial under certain boundary conditions.

    The polynomial is:
        P(dr) = alpha + beta*dr + gamma*dr^2 + delta*dr^3 + epsilon*dr^4 + phi*dr^5
    where dr = r - r0.

    Args:
        derivs (tuple of float): Usually `(aa, bb, cc)` representing constraints
            on the polynomial and its derivatives at r0, e.g., P(r0), P'(r0), P''(r0).
        r0 (float): Reference position.
        rcut (float): Position where P, P', P'' are set to zero (e.g. cutoff).

    Returns:
        tuple of float: The six polynomial coefficients
            `(alpha, beta, gamma, delta, epsilon, phi)`.
    """

    aa, bb, cc = derivs
    alpha = aa
    beta = bb
    gamma = cc / 2.0
    dr = rcut - r0

    mtx = np.array([[ dr**3, dr**4, dr**5 ],
                    [ 3.0 * dr**2, 4.0 * dr**3, 5.0 * dr**4 ],
                    [ 6.0 * dr, 12.0 * dr**2, 20.0 * dr**3 ]])
    rhs = np.array([ -(alpha + beta * dr + gamma * dr**2),
                     -(beta + 2 * gamma * dr),
                     -(2.0 * gamma) ])
    delta, epsilon, phi = linalg.solve(mtx, rhs)

    return alpha, beta, gamma, delta, epsilon, phi


def get_poly5_values(coeffs, r0, rr):
    """
    Evaluate a 5th-degree polynomial using Horner's method.

    The polynomial has the form:
        P(dr) = c0 + c1*dr + c2*dr^2 + c3*dr^3 + c4*dr^4 + c5*dr^5,
    where dr = rr - r0.

    Args:
        coeffs (array_like): Length-6 array of polynomial coefficients
            [c0, c1, c2, c3, c4, c5].
        r0 (float): Reference shift, used in dr = rr - r0.
        rr (float or array_like): Point(s) at which to evaluate the polynomial.

    Returns:
        float or ndarray: The polynomial value(s). If `rr` is scalar, returns float;
        otherwise returns an ndarray.
    """

    dr = rr - r0
    res = coeffs[-1]

    for ii in range(len(coeffs) - 2, -1, -1):
        res = res * dr + coeffs[ii]

    return res


def get_expcoeffs(derivs, r0):
    """
    Compute coefficients for an exponential function that satisfies given conditions at r0.

    The function is:
        f(r) = exp(-alpha*r + beta) + gamma.

    Args:
        derivs (tuple of float): Typically (aa, bb, cc) representing constraints
            such as f(r0), f'(r0), f''(r0).
        r0 (float): Reference point where conditions are specified.

    Returns:
        tuple: (alpha, beta, gamma), where
            - alpha (float): Decay rate
            - beta (float): Exponential offset in the exponent
            - gamma (float): Vertical shift
    """

    aa, bb, cc = derivs
    alpha = -cc / bb
    beta = alpha * r0 + np.log(cc / alpha**2)
    gamma = aa - np.exp(-alpha * r0 + beta)

    return alpha, beta, gamma


def get_exp_values(coeffs, rr):
    """
    Evaluate the exponential function f(r) = exp(-alpha*r + beta) + gamma.

    Args:
        coeffs (tuple of float): (alpha, beta, gamma)
        rr (float or array_like): Point(s) at which to evaluate the function.

    Returns:
        float or ndarray: Function value(s). Float if `rr` is scalar,
        otherwise an ndarray.
    """

    aa, bb, cc = coeffs

    return np.exp(-aa * rr + bb) + cc


def write_splinerep(fname, expcoeffs, splcoeffs, poly5coeffs, rr, rcut):
    """
    Write a spline representation of a repulsive potential to a file.

    The file format is:
    1) "Spline"
    2) Number of spline points and cutoff
    3) Exponential coefficients (3 floats)
    4) Spline coefficients for each interval in `rr`
    5) 5th-degree polynomial coefficients for the last segment

    Args:
        fname (str): Output filename.
        expcoeffs (tuple or list): (alpha, beta, gamma) for the exponential portion.
        splcoeffs (array_like): Spline coefficient data for each interval in `rr`.
            Typically shape (len(rr)-1, 4).
        poly5coeffs (tuple or list): 6 polynomial coefficients
            for the 5th-degree segment at the end.
        rr (array_like): The list of spline breakpoints.
        rcut (float): The cutoff radius.
    """

    fp = open(fname, 'w')
    fp.write('Spline\n')
    fp.write('{:d} {:.4f}\n'.format(len(rr), rcut))
    fp.write('{:15.8E} {:15.8E} {:15.8E}\n'.format(*expcoeffs))
    splcoeffs_format = ' '.join(['{:6.3f}'] * 2 + ['{:15.8E}'] * 4) + '\n'
    for ir in range(len(rr) - 1):
        rcur = rr[ir]
        rnext = rr[ir + 1]
        fp.write(splcoeffs_format.format(rcur, rnext, *splcoeffs[ir]))
    poly5coeffs_format = ' '.join(['{:6.3f}'] * 2 + ['{:15.8E}'] * 6) + '\n'
    fp.write(poly5coeffs_format.format(rr[-1], rcut, *poly5coeffs))
    fp.close()
    print_io_log(fname, 'Repulsive in spline format')


def write_as_nxy(fname, datadesc, vectors, column_names):
    """
    Write arrays in column-oriented format to a text file with descriptive headers.

    Args:
        fname (str): Output filename.
        datadesc (str): Short description of the data (for logging).
        vectors (list of array_like): Each element is a 1D array for one column.
        column_names (list of str): Names or labels for each column in `vectors`.
    """

    header_parts = []
    for ii, colname in enumerate(column_names):
        header_parts.append('Column {}: {}'.format(ii + 1, colname))
    header = '\n'.join(header_parts)
    data = np.array(vectors).transpose()
    np.savetxt(fname, data, header=header)
    print_io_log(fname, datadesc)


def print_io_log(fname, fcontent):
    """
    Print a short log message for file I/O operations.

    Args:
        fname (str): The filename or path used in the operation.
        fcontent (str): Description of the operation or file content.
    """

    print(f"{fcontent} -> '{fname}'")


def main():
    """
    Main routine of the splinerepfit tool.
    """

    parser = get_cmdline_parser()
    args = parser.parse_args()
    refdata = args.refdata
    dftbdata = args.dftbdata
    rmin = args.rmin
    rcut = args.rcut
    rdamp = args.rdamp
    interp_grid = args.intergrid
    spline_grid = args.splgrid
    minpow = args.minpow
    maxpow = args.maxpow
    poly5len = args.taillen

    # Read and interpolate energy curves
    rref, eref = get_data(refdata)
    rdftb, edftb = get_data(dftbdata)
    fref = spinter.interp1d(rref, eref)
    fdftb = spinter.interp1d(rdftb, edftb)

    # Shift energy curves together to match at cutoff
    rmin_raw = max(rref[0], rdftb[0])
    rmax_raw = min(rref[-1], rdftb[-1])
    rr = np.arange(rmin_raw + EPSILON, rmax_raw - EPSILON, interp_grid)
    data1 = fref(rr) - fref(rcut)
    data2 = fdftb(rr) - fdftb(rcut)
    write_as_nxy('shifted.dat', 'Shifted energy profiles', (rr, data1, data2),
                 ('rr', 'reference curve (shifted)', 'dftb curve (shifted)'))

    # Build energy differences with smoothing
    rr = np.arange(rmin, rcut + EPSILON, interp_grid)
    eref = fref(rr)
    edftb = fdftb(rr)
    diff = eref - edftb - (eref[-1] - edftb[-1])
    if rdamp:
        diff_smooth = diff * damping_cos(rdamp, rcut, rr)
    else:
        diff_smooth = diff
    fdesc = 'Target repulsive'
    write_as_nxy('targetrep.dat', fdesc, (rr, diff_smooth),
                 ('rr', 'target repulsive'))

    # Fit polynomial on energy difference
    fitfunc = lambda cc, rr: polynomial(cc, minpow, rr, rcut)
    errfunc = lambda cc, rr, ref: fitfunc(cc, rr) - ref
    c0 = np.zeros((maxpow - minpow + 1), dtype=float)
    cc, success = spopt.leastsq(errfunc, c0, args=(rr, diff_smooth))
    print('Least square fit succesfull:', bool(success))

    # Write fitted repulsive and reference data
    erep = polynomial(cc, minpow, rr, rcut)
    write_as_nxy('polyrepfit.dat', 'Polynomial fit on target repulsive',
                 (rr, erep), ('rr', 'fitted repulsive'))

    # Fit spline
    rr = np.arange(rmin, rcut - poly5len + EPSILON, interp_grid)
    rspline = np.arange(rmin, rcut - poly5len + EPSILON, spline_grid)
    erep = polynomial(cc, minpow, rr, rcut)
    erepc = polynomial(cc, minpow, rspline, rcut)
    ccderiv, minpowderiv = get_polyderiv(1, cc, minpow)
    deriv0 = polynomial(ccderiv, minpowderiv, rspline[0:1], rcut)[0]
    derivn = polynomial(ccderiv, minpowderiv, rspline[-1:], rcut)[0]
    cs, splcoeffs = get_spline_coeffs(rspline, erepc, deriv0, derivn)
    splval = cs(rr)
    write_as_nxy('splinefit.dat', 'Spline fitted on polynomial fit',
                 (rr, splval), ('rr', 'fitted spline'))

    # Fit exponential start to spline
    splderivs = get_splineval012(splcoeffs[0], rspline[0], rspline[0])
    expcoeffs = get_expcoeffs(splderivs, rspline[0])
    rexp = np.arange(rmin - 0.5, rmin + EPSILON, interp_grid)
    expvals = get_exp_values(expcoeffs, rexp)
    write_as_nxy('headfit.dat', 'Exponentail head', (rexp, expvals),
                 ('rr', 'exponential head'))

    # Fit 5th order polynomial tail
    derivs = get_splineval012(splcoeffs[-1], rspline[-2], rspline[-1])
    poly5coeffs = get_poly5coeffs(derivs, rspline[-1], rcut)
    rpoly5 = np.arange(rcut - poly5len, rcut + EPSILON, interp_grid)
    poly5vals = get_poly5_values(poly5coeffs, rspline[-1], rpoly5)
    write_as_nxy('tailfit.dat', '5th order spline tail', (rpoly5, poly5vals),
                 ('rr', '5th order spline'))

    # Write SK-compatible representation
    write_splinerep('repulsive.spl', expcoeffs, splcoeffs, poly5coeffs,
                    rspline, rcut)

    # Write shifted etot curves
    rmin0 = rref[0]
    rmax0 = rref[-1]
    rr = []
    ff = []
    rvals = np.arange(rmin0, rmin + EPSILON, interp_grid)
    if len(rvals):
        rr.append(rvals)
        ff.append(get_exp_values(expcoeffs, rvals))
    rvals = np.arange(rmin, min(rcut - poly5len + EPSILON, rmax0),
                        interp_grid)
    if len(rvals):
        rr.append(rvals)
        ff.append(cs(rvals))
    rvals = np.arange(rcut - poly5len, min(rmax0, rcut), interp_grid)
    if len(rvals):
        rr.append(rvals)
        ff.append(get_poly5_values(poly5coeffs, rspline[-1], rvals))
    rr = np.hstack(rr)
    etot_dftb = np.hstack(ff) + fdftb(rr)
    etot_ref = fref(rr)
    imin = np.argmin(etot_ref)
    etot_ref -= etot_ref[imin]
    imin = np.argmin(etot_dftb)
    etot_dftb -= etot_dftb[imin]
    write_as_nxy('etotcomp.dat', 'Reference vs. DFTB total energy',
                 (rr, etot_ref, etot_dftb), ('rr', 'reference total energy',
                                             'dftb total energy'))


if __name__ == '__main__':
    main()
