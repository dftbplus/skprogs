#!/usr/bin/env python3

'''General tool for manipulating SK-tables.'''


import sys
import argparse
import re
import xml.etree.ElementTree as etree
from sktools import PACKAGE_VERSION
import sktools.common as sc
from sktools.oldskfile import OldSKFile

SCRIPTNAME = sc.get_script_name()


FNAME_PATTERN = re.compile('(?P<elem1>\w+)-(?P<elem2>\w+)\.skf')


def main(cmdlineargs=None):
    '''Main driver routine.'''

    parser, subparsers = get_parser_and_subparser_container()
    setup_parser_main(parser)
    common = get_common_parser()
    setup_parser_getdoc(subparsers, common, run_getdoc)
    setup_parser_setdoc(subparsers, common, run_setdoc)
    parse_command_line_and_run_subcommand(parser, cmdlineargs)


def run_getdoc(args):
    '''Extracts documentation out of SK-file.'''

    skfile = args.skfile
    if args.sktype == 'auto':
        homo = is_homo_file(skfile)
    else:
        homo = (args.sktype == 'homo')
    sk = OldSKFile.fromfile(skfile, homo)
    doc = sk.documentation
    fobj = sys.stdout if args.file == '-' else args.file
    fp, tobeclosed = sc.openfile(fobj, 'w')
    fp.write(etree.tostring(doc, encoding='UTF-8').decode('UTF-8'))
    if tobeclosed:
        fp.close()


def run_setdoc(args):
    '''Sets documentation in SK-file.'''

    skfile = args.skfile
    if args.sktype == 'auto':
        homo = is_homo_file(skfile)
    else:
        homo = (args.sktype == 'homo')
    sk = OldSKFile.fromfile(skfile, homo)
    fobj = sys.stdin if args.file == '-' else args.file
    fp, tobeclosed = sc.openfile(fobj, 'r')
    xml = fp.read()
    if tobeclosed:
        fp.close()
    doc = etree.fromstring(xml)
    sk.documentation = doc
    sk.tofile(skfile)


def is_homo_file(filename):
    '''Infers if given SK-file is homonuclear.'''

    match = FNAME_PATTERN.match(filename)
    if match:
        homo = (match.group('elem1') == match.group('elem2'))
    else:
        homo = False

    return homo


def get_parser_and_subparser_container():
    '''Returns parser and subparser of skmanip.'''

    parser = argparse.ArgumentParser(
        description='General tool for manipulating SK-tables.')
    subparsers = parser.add_subparsers(title='available subcommands',
                                       help='')
    return parser, subparsers


def get_common_parser():
    '''Common settings for all one-center calculations.'''

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('skfile', help='skfile to process.')
    common.add_argument(
        '-t', '--type', dest='sktype', choices=['homo', 'hetero', 'auto'],
        default='auto', help='Type of skfile (default: auto)')
    common.add_argument(
        '-f', '--file', default='-',
        help='Reads/writes from/into file instead using stdin/stderr')

    return common


def setup_parser_main(parser):
    '''Add skmanip version to main parser.'''

    parser.add_argument('--version', action='version',
                        version='skmanip {}'.format(PACKAGE_VERSION))


def setup_parser_getdoc(subparsers, common, target_function):
    '''setup_parser_getdoc'''

    parser = subparsers.add_parser(
        'get_documentation', parents=[common],
        help='Extracts the documentation into a file')
    parser.set_defaults(func=target_function)


def setup_parser_setdoc(subparsers, common, target_function):
    '''setup_parser_setdoc'''

    parser = subparsers.add_parser(
        'set_documentation', parents=[common],
        help='Replaces the documentation in an SK-file')
    parser.set_defaults(func=target_function)


def parse_command_line_and_run_subcommand(parser, cmdlineargs):
    '''Parses command line arguments and runs skmanip subcommand.'''

    args = parser.parse_args(args=cmdlineargs)
    args.func(args)


if __name__ == '__main__':
    try:
        main()
    except sc.SkgenException as ex:
        sc.fatalerror(str(ex))
