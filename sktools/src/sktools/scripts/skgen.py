#!/usr/bin/env python3

'''
Module to generate homo- and hetero-nuclear Slater-Koster files.
'''


import sys
import argparse
import numpy as np

import sktools.common as sc
from sktools import skgen
from sktools.skdef import Skdef
from sktools import PACKAGE_VERSION


if sys.hexversion < 0x03020000:
    sys.exit('Program only works with Python 3.2 or greater')

if np.__version__.startswith('1.6.'):
    sys.exit('Program only works with Numpy 1.7.x or greater')


SCRIPTNAME = sc.get_script_name()

# Global script logger, will be overriden by the setup_logger() method in
# the respective subcommands depending on the command line loglevel options
logger = None


def main(args=None):
    '''Main driver routine.'''

    parser, subparsers = get_parser_and_subparser_container()
    setup_parser_main(parser)
    onecnt_common = get_onecnt_common_parser()
    setup_parser_atom(subparsers, onecnt_common, run_atom)
    setup_parser_denscomp(subparsers, onecnt_common, run_denscomp)
    setup_parser_wavecomp(subparsers, onecnt_common, run_wavecomp)
    twocnt_common = get_twocnt_common_parser()
    setup_parser_twocnt(subparsers, twocnt_common, run_twocnt)
    setup_parser_sktable(subparsers, twocnt_common, run_sktable)
    parse_command_line_and_run_subcommand(parser, args)


def run_atom(args):
    setup_logger(args.loglevel)
    logger.info('Subcommand atom started')
    elements = convert_argument_to_elements(args.element)
    skdefs = merge_skdefs(args.configfiles)
    searchdirs = [args.builddir,] + args.includedirs
    resultdirs = []
    for elem in elements:
        calculator = skgen.run_atom(
            skdefs, elem, args.builddir, searchdirs, args.onecnt_binary,
            args.eigenonly, args.eigenspinonly)
        resultdirs.append(calculator.get_result_directory())
    logger.info('Subcommand atom finished')
    logger.info('Atom results in {}'.format(' '.join(resultdirs)))


def run_denscomp(args):
    setup_logger(args.loglevel)
    logger.info('Subcommand denscomp started')
    elements = convert_argument_to_elements(args.element)
    skdefs = merge_skdefs(args.configfiles)
    searchdirs = [args.builddir,] + args.includedirs
    resultdirs = []
    for elem in elements:
        calculator = skgen.run_denscomp(
            skdefs, elem, args.builddir, searchdirs, args.onecnt_binary)
        resultdirs.append(calculator.get_result_directory())
    logger.info('Subcommand densecomp finished')
    logger.info('Denscomp results in {}'.format(' '.join(resultdirs)))


def run_wavecomp(args):
    setup_logger(args.loglevel)
    logger.info('Subcommand wavecomp started')
    elements = convert_argument_to_elements(args.element)
    skdefs = merge_skdefs(args.configfiles)
    searchdirs = [args.builddir,] + args.includedirs
    resultdirs = []
    for elem in elements:
        calculator = skgen.run_wavecomp(
            skdefs, elem, args.builddir, searchdirs, args.onecnt_binary)
        dirnames = ' '.join(calculator.get_result_directories())
        resultdirs.append(dirnames)
    logger.info('Subcommand wavecomp finished')
    logger.info('Wavecomp results in {}'.format(' '.join(resultdirs)))


def run_twocnt(args):
    setup_logger(args.loglevel)
    logger.info('Subcommand twocnt started')
    skdefs = merge_skdefs(args.configfiles)
    builddir = args.builddir
    searchdirs = [builddir,] + args.includedirs
    resultdirs = []
    element_pairs = convert_arguments_to_element_pairs(args.element1,
                                                       args.element2)
    for elem1, elem2 in element_pairs:
        calculator = skgen.run_twocnt(
            skdefs, elem1, elem2, builddir, searchdirs, args.onecnt_binary,
            args.twocnt_binary)
        resultdirs.append(calculator.get_result_directory())
    logger.info('Subcommand twocnt finished')
    logger.info('Twocnt results in {}'.format(' '.join(resultdirs)))


def run_sktable(args):
    setup_logger(args.loglevel)
    logger.info('Subcommand sktable started')
    skdefs = merge_skdefs(args.configfiles)
    builddir = args.builddir
    searchdirs = [builddir,] + args.includedirs
    workdir = args.outdir
    add_dummy_rep = args.dummyrep
    skfiles_written = []
    element_pairs = convert_arguments_to_element_pairs(args.element1,
                                                       args.element2)
    for elem1, elem2 in element_pairs:
        skfiles_written += skgen.run_sktable(
            skdefs, elem1, elem2, builddir, searchdirs, args.onecnt_binary,
            args.twocnt_binary, workdir, add_dummy_rep)
    logger.info('Directory with assembled SK-file(s): {}'.format(workdir))
    logger.info('SK-file(s) written: {}'.format(' '.join(skfiles_written)))


def get_parser_and_subparser_container():
    parser = argparse.ArgumentParser(
        description='General tool for generating Slater-Koster tables.')
    subparsers = parser.add_subparsers(title='available subcommands',
                                       help='')
    return parser, subparsers


def get_onecnt_common_parser():
    '''Common settings for all one-center calculations.'''
    onecnt_common = argparse.ArgumentParser(add_help=False)
    onecnt_common.add_argument(
        'element', help='element to process: either one element (e.g. N) or a '
        'comma separated list of element names *without* spaces in between '
        '(e.g. N,C,H)')
    return onecnt_common


def get_twocnt_common_parser():
    twocnt_common = argparse.ArgumentParser(add_help=False)
    twocnt_common.add_argument(
        'element1', help='first element of the element pair to process: '
        'either one element (e.g. N) or a comma separated list of element '
        'names *without* spaces in between (e.g. N,C,H)')
    twocnt_common.add_argument(
        'element2', help='second element of the element pair to process: '
        'either one element (e.g. N) or a comma separated list of element '
        'names *without* spaces in between (e.g. N,C,H)')
    return twocnt_common


def setup_parser_main(parser):
    parser.add_argument('--version', action='version',
                        version='sktools {}'.format(PACKAGE_VERSION))
    parser.add_argument(
        '-I', '--include-dir', action='append', default=[],
        dest='includedirs',
        help='directory to include in the search for calculation '
        '(default: build directory only)')
    parser.add_argument(
        '-c', '--config-file', default='skdef.hsd', dest='configfiles',
        metavar='CONFIGFILE',
        help='config file to be parsed (default: ./skdef.hsd)'
    )
    parser.add_argument(
        '-b', '--build-dir', default='_build', dest='builddir',
        help='build directory (default: _build)')
    parser.add_argument(
        '-o', '--onecenter-binary', dest='onecnt_binary', default=None,
        help='binary to use for the one-center calculations (default: depends '
        'on the calculator specified in the input)')
    parser.add_argument(
        '-t', '--twocenter-binary', dest='twocnt_binary', default=None,
        help='binary to use for the two-center calculationrs (default: depends '
        'on the calculator speciefied in the input)')
    parser.add_argument(
        '-l', '--log-level', dest='loglevel', default='info',
        choices=['debug', 'info', 'warning', 'error'],
        help='Logging level (default: info)')


def setup_parser_atom(subparsers, onecnt_common, target_function):
    parser_atom = subparsers.add_parser(
        'atom', parents=[onecnt_common],
        help='calculates the free atom to get eigenlevels, hubbard values, spin'
        ' couplings, etc.')
    parser_atom.add_argument(
        '-e', '--eigenlevels-only', dest='eigenonly', action='store_true',
        default=False, help='calculates only eigenlevels of the spin '
                            'unpolarized atom but no derivatives.')
    parser_atom.add_argument(
        '-s', '--spin-polarized', dest='eigenspinonly', action='store_true',
        default=False, help='calculates only the eigenlevels of the spin '
                            'polarized atom but no derivatives')
    parser_atom.set_defaults(func=target_function)


def setup_parser_denscomp(subparsers, onecnt_common, target_function):
    parser_denscomp = subparsers.add_parser(
        'denscomp', parents=[onecnt_common],
        help='calculates density compression')
    parser_denscomp.set_defaults(func=target_function)


def setup_parser_wavecomp(subparsers, onecnt_common, target_function):
    parser_wavecomp = subparsers.add_parser(
        'wavecomp', parents=[onecnt_common],
        help='calculates wave function compression')
    parser_wavecomp.set_defaults(func=target_function)


def setup_parser_twocnt(subparsers, twocnt_common, target_function):
    parser_twocnt = subparsers.add_parser(
        'twocnt', parents=[twocnt_common],
        help='calculates two center integrals')
    parser_twocnt.set_defaults(func=target_function)


def setup_parser_sktable(subparsers, twocnt_common, target_function):
    parser_sktable = subparsers.add_parser(
        'sktable', parents=[twocnt_common],
        help='creates an sktable for a given element pair')
    parser_sktable.add_argument(
        '-d', '--dummy-repulsive', action='store_true', dest='dummyrep',
        default=False, help='add dummy repulsive spline to the sk tables')
    parser_sktable.add_argument(
        '-o', '--output-dir', dest='outdir', default='.',
        help='directory where the skfiles should be written to (default: .)')
    parser_sktable.set_defaults(func=target_function)


def parse_command_line_and_run_subcommand(parser, cmdlineargs=None):
    args = parser.parse_args(args=cmdlineargs)
    args.configfiles = [args.configfiles,]
    args.func(args)


def setup_logger(loglevel):
    global logger
    logger = sc.get_script_logger(loglevel, SCRIPTNAME)


def merge_skdefs(filenames):
    '''Returns a merged skdefs object using all specified skdef files.'''

    skdef = Skdef.fromfile(filenames[0])
    for filename in filenames[1:]:
        skdef2 = Skdef.fromfile(filename)
        skdef.update(skdef2)
    return skdef


def convert_argument_to_elements(argument):
    return argument.split(',')


def convert_arguments_to_element_pairs(argument1, argument2):
    elements1 = convert_argument_to_elements(argument1)
    elements2 = convert_argument_to_elements(argument2)
    processed = set()
    element_pairs = []
    for elem1 in elements1:
        elem1low = elem1.lower()
        for elem2 in elements2:
            elem2low = elem2.lower()
            already_processed = ((elem1low, elem2low) in processed
                                 or (elem2low, elem1low) in processed)
            if not already_processed:
                element_pairs.append((elem1, elem2))
                processed.add((elem1low, elem2low))
    return element_pairs


if __name__ == '__main__':
    try:
        main()
    except sc.SkgenException as ex:
        sc.fatalerror(str(ex))
