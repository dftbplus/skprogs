#!/usr/bin/env python3

'''Collects spin coupling constants.'''

import argparse
from sktools.skdef import Skdef
from sktools import skgen
import sktools.common as sc


USAGE = \
    '''
    Collects spin coupling constants by iterating over the elements defined
    in skdef.hsd. If the atomic calculation has been done already, it will be
    reused, otherwised it is done on the fly.
    '''

SCRIPTNAME = sc.get_script_name()
SPINW_FILE_NAME = 'spinw.txt'


def main(cmdlineargs=None):
    '''Main driver routine.'''

    args = parseargs(cmdlineargs)

    logger = sc.get_script_logger(args.loglevel, SCRIPTNAME)
    logger.info('Collecting spinw constants')

    skdef = Skdef.fromfile(args.configfile)
    searchdirs = [args.builddir]
    elems = skdef.atomparameters.keys()

    with open(SPINW_FILE_NAME, 'w', encoding='utf8') as spinw_fd:

        for elem in elems:
            calculator = skgen.run_atom(skdef, elem, args.builddir, searchdirs,
                                        args.onecnt_binary)
            spinw_fd.write(sc.capitalize_elem_name(elem) + ':\n')
            results = calculator.get_result()
            spinw = results.get_spinws()
            ndim = spinw.shape[0]
            formstr = '{:13.5f}' * ndim + '\n'
            for line in spinw:
                spinw_fd.write(formstr.format(*line))
                spinw_fd.write('\n')

    logger.info(f"File '{SPINW_FILE_NAME}' written.")


def parseargs(cmdlineargs):
    '''Parses command line arguments and return the parser instance.'''

    parser = argparse.ArgumentParser(description=USAGE)

    msg = 'build directory (default: _build)'
    parser.add_argument('-b', '--build-dir', default='_build', dest='builddir',
                        help=msg)

    msg = 'binary to use for the one-center calculations (default: depends ' + \
        'on the calculator specified in the input)'
    parser.add_argument('-o', '--onecenter-binary', dest='onecnt_binary',
                        default=None, help=msg)

    parser.add_argument(
        '-c', '--config-file', default='skdef.hsd', dest='configfile',
        metavar='CONFIGFILE',
        help='config file to be parsed (default: ./skdef.hsd)')

    msg = 'Logging level (default: info)'
    parser.add_argument('-l', '--log-level', dest='loglevel', default='info',
                        choices=['debug', 'info', 'warning', 'error'], help=msg)

    return parser.parse_args(args=cmdlineargs)


if __name__ == '__main__':
    main()
