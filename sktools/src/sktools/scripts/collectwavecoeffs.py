#!/usr/bin/env python3

'''
Module to extract or calculate wavefunction coefficients for Waveplot.
'''


import os
import argparse

import sktools.common as sc
from sktools.skgen import compression
from sktools.scripts.skgen import get_onecnt_common_parser, \
    setup_parser_wavecomp, \
    parse_command_line_and_run_subcommand, \
    convert_argument_to_elements, merge_skdefs
from sktools import PACKAGE_VERSION
from sktools.oldskfile import OldSKFile
from sktools.taggedfile import TaggedFile


SCRIPTNAME = sc.get_script_name()

USAGE = \
    '''Collects coefficient information for Waveplot. It iterates over the
    elements defined in skdef.hsd and collects the wavefunction coefficients
    and other information necessary for Waveplot. The homonuclear SK-files as
    well as atomic calculations with compressed wavefunction must have been
    generated in advance.
    '''

# Global script logger, will be overriden by the setup_logger() method in
# the respective subcommands depending on the command line loglevel options
logger = None


def main(args=None):
    '''Main driver routine.'''

    parser, subparsers = get_parser_and_subparser_container()
    setup_parser_main(parser)
    onecnt_common = get_onecnt_common_parser()
    setup_parser_wavecomp(subparsers, onecnt_common, run_wavecomp)
    parse_command_line_and_run_subcommand(parser, args)


def writecoeffs(fd_wc, elem, atomconfig, homoskname, wavecompdir):
    '''Writes element-specific input, processed by Waveplot.

    Args:

        fd_wc (file object): file object to write to
        elem (str): element name to fetch information for
        atomconfig (AtomConfig): represents the configuration of a free atom
        homoskname (str): pathname of homonuclear Slater-Koster file
        wavecompdir (str): path to calculation of the compressed atom

    '''

    homosk = OldSKFile.fromfile(homoskname, True)
    cutoff = homosk.nr * homosk.dr / 2.0
    fd_wc.write('{} {{\n'.format(elem))
    fd_wc.write('  AtomicNumber = {:d}\n'.format(int(atomconfig.atomicnumber)))
    for nn, ll in atomconfig.valenceshells:
        coeffsname = 'coeffs_{:02d}{:1s}.tag'.format(nn, sc.ANGMOM_TO_SHELL[ll])
        coeffs = TaggedFile.fromfile(os.path.join(wavecompdir, coeffsname),
                                     transpose=True)
        fd_wc.write('  Orbital {\n')
        fd_wc.write('    AngularMomentum = {:d}\n'.format(ll))
        fd_wc.write('    Occupation = {:.1f}\n'.format(coeffs['occupation']))
        fd_wc.write('    Cutoff = {:5.2f}\n'.format(cutoff))
        fd_wc.write('    Exponents {\n')
        sc.writefloats(fd_wc, coeffs['exponents'], indent=6, numperline=3,
                       formstr='{:21.12E}')
        fd_wc.write('    }\n')
        fd_wc.write('    Coefficients {\n')
        sc.writefloats(fd_wc, coeffs['coefficients'], indent=3, numperline=3,
                       formstr='{:21.12E}')
        fd_wc.write('    }\n')
        fd_wc.write('  }\n')
    fd_wc.write('}\n')


def run_wavecomp(args):
    setup_logger(args.loglevel)
    logger.info('Checking if required SK-files are present')
    elements = convert_argument_to_elements(args.element)
    elements_lower = [elem.lower() for elem in elements]

    for ielem, elem in enumerate(elements_lower):
        homoskname = '{elem}-{elem}.skf'.format(elem=elem.capitalize())
        if not os.path.exists(homoskname):
            sc.fatalerror('Error while processing element ' + elem.capitalize()
                          + '. SK-file ' + homoskname + ' absent.')

    skdef = merge_skdefs(args.configfiles)
    searchdirs = [args.builddir,] + args.includedirs
    resultdirs = []

    logger.info('Subcommand wavecomp started')
    for elem in elements:
        calculator = compression.search_wavecomp(
            skdef, elem, args.builddir, searchdirs)
        dirnames = ' '.join(calculator.get_result_directories())
        resultdirs.append(dirnames)

    logger.info('Subcommand wavecomp finished')
    logger.info('Wavecomp results in {}'.format(' '.join(resultdirs)))

    atomparameters = skdef.atomparameters

    with open('wfc.hsd', 'w', encoding='utf8') as fd_wc:

        for ielem, elem in enumerate(elements_lower):
            homoskname = '{elem}-{elem}.skf'.format(elem=elem.capitalize())
            wavecompdir = os.path.join(os.getcwd(), resultdirs[ielem])
            if not os.path.exists(homoskname):
                sc.fatalerror('Error while processing element '
                              + elem.capitalize() + '. SK-file ' + homoskname
                              + ' absent.')

            logger.info('Writing element ' + elem.capitalize()
                        + ' to wfc.hsd file.')
            atomconfig = atomparameters[elem].atomconfig
            writecoeffs(fd_wc, elem.capitalize(), atomconfig, homoskname,
                        wavecompdir)


def setup_parser_main(parser):
    parser.add_argument('--version', action='version',
                        version=f'sktools {PACKAGE_VERSION}')

    parser.add_argument(
        '-I', '--include-dir', action='append', default=[],
        dest='includedirs',
        help='directory to include in the search for calculation '
        '(default: build directory only)')

    parser.add_argument(
        '-c', '--config-file', default='skdef.hsd', dest='configfiles',
        metavar='CONFIGFILE',
        help='config file to be parsed (default: ./skdef.hsd)')

    parser.add_argument(
        '-b', '--build-dir', default='_build', dest='builddir',
        help='build directory (default: _build)')

    parser.add_argument(
        '-l', '--log-level', dest='loglevel', default='info',
        choices=['debug', 'info', 'warning', 'error'],
        help='Logging level (default: info)')


def get_parser_and_subparser_container():
    parser = argparse.ArgumentParser(description=USAGE)
    subparsers = parser.add_subparsers(
        title='available subcommands', help='')
    return parser, subparsers


def setup_logger(loglevel):
    global logger
    logger = sc.get_script_logger(loglevel, SCRIPTNAME)


if __name__ == '__main__':
    try:
        main()
    except sc.CollectwavecoeffsException as ex:
        sc.fatalerror(str(ex))
