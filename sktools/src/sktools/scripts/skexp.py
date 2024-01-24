#!/usr/bin/env python3

'''
Prints STO exponents according to the usual geometric series.
(see for example DOI: 10.1063/1.1679962)
'''

import argparse
import numpy as np


USAGE = \
    '''
    Prints a set of Slater-type orbital exponents, following the usual geometric
    series for a given atomic number.
    '''


def main(cmdlineargs=None):
    '''Main driver routine.'''

    args = parseargs(cmdlineargs)

    exponents = get_exponents(
        args.atnum, args.nexp, a0=args.a0, extra=args.extra)

    print(((args.nexp + int(args.extra)) * ' {:.4f}').format(*exponents))


def parseargs(cmdlineargs):
    '''Parses command line arguments and returns the parser instance.'''

    parser = argparse.ArgumentParser(description=USAGE)

    msg = 'atomic number of element to parametrize'
    parser.add_argument('atnum', metavar='Z', type=int, help=msg)

    msg = 'number of exponents'
    parser.add_argument('nexp', metavar='N', type=int, help=msg)

    msg = 'smallest exponent (default: 0.5)'
    parser.add_argument('-s', '--start', metavar='a0', dest='a0',
                        default=0.5, type=float, help=msg)

    msg = 'adds one additional exponent beyond the atomic number, proved to' + \
        ' be useful for large atoms like e.g. Au (default: False)'
    parser.add_argument(
        '-e', '--extra', dest='extra', action='store_true', help=msg)

    args = parser.parse_args(args=cmdlineargs)

    # check input consistence
    if args.atnum < 1 or args.atnum > 118:
        raise ValueError('Invalid atomic number, choose 0 < Z < 119.')

    if args.nexp < 1:
        raise ValueError('Invalid number of STO exponents, choose N > 0.')

    if args.a0 <= 0.0:
        raise ValueError('Invalid smallest exponents, choose positive a0.')

    return args


def get_exponents(atnum, nexp, a0=0.5, extra=False):
    '''Generates Slater exponents according to geometric series.

        a0, a0*r, a0*r**2, a0*r**3, a0*r**4, ...

    Args:

        atnum (int): atomic number
        nexp (int): number of exponents (excluding extra one)
        a0 (float): smallest exponent, usually 0.5 or similar
        extra (bool): add one more exponent beyond atomic number,
            proved to useful for larger atoms like e.g. Au

    Returns:

        exponents (1darray): generated STO exponents

    '''

    rr = (a0 / float(atnum))**(1.0 / float((nexp - 1)))

    if extra:
        nn = nexp + 1
    else:
        nn = nexp

    exponents = np.empty(nn, dtype=float)

    if extra:
        for iexp in range(nn):
            exponents[nn - iexp - 1] = atnum * rr**float(iexp - 1)
    else:
        for iexp in range(nn):
            exponents[nn - iexp - 1] = atnum * rr**float(iexp)

    return exponents


if __name__ == '__main__':
    main()
