'''
Basic sktable Regression Test Class

Initialized with the path to a working directory this class figures out the
rudimentary specifications of an sktable run and compares the existing
reference file with the output of the current version of the program.
'''


import os
import warnings
from collections import Counter

from sktools.oldskfile import OldSKFile


SKDEF = 'skdef.hsd'

ATOL = 1e-10
RTOL = 1e-09


class TestWithWorkDir:
    '''Basic sktable Regression Test Class.'''


    def __init__(self, wdir):
        '''Initializes a TestWithWorkDir object.

        Args:

            wdir (str): path to working directory to perform regression testing

        '''

        self._wdir = wdir

        # check for skdef file
        if not os.path.isfile(os.path.join(self._wdir, SKDEF)):
            raise TestWithWorkDirError('Skdef file absent.')

        # check for reference SK-files and additional junk
        allfiles = os.listdir(self._wdir)
        reffiles = [fname for fname in os.listdir(self._wdir)
                    if fname.startswith('_') and fname.endswith('.skf')]
        if len(reffiles) < 1:
            raise TestWithWorkDirError('No reference SK-files found.')

        # check if there is any additional junk
        if Counter(reffiles) != Counter(allfiles):
            exceptions = ['config', '_build', SKDEF]
            curfiles = [fname.strip('_') for fname in reffiles]
            junkfiles = [fname for fname in allfiles
                         if fname not in reffiles
                         and fname not in curfiles
                         and fname not in exceptions
                         and not fname.endswith('~')]
            if len(junkfiles) > 0:
                warnings.warn('Found additional junk in the directory: {}'
                              .format(junkfiles))

        # store list of SK-files to compare (without '_' at the beginning)
        self._skfiles = [fname.strip('_') for fname in reffiles]


    def test(self):
        '''Performs regression testing by comparing all relevant files.

        Returns:

            passed (bool): True, if all tests have passed, otherwise False

        '''

        skpassed = True

        for skfile in self._skfiles:
            ishomo = infer_homonuclear(skfile)
            if ishomo:
                print('Comparing homonuclear files ({}, {}):'
                      .format(skfile, '_' + skfile))
            else:
                print('Comparing heteronuclear files ({}, {}):'
                      .format(skfile, '_' + skfile))
            curskfile = OldSKFile.fromfile(
                os.path.join(self._wdir, skfile), ishomo)
            refskfile = OldSKFile.fromfile(
                os.path.join(self._wdir, '_' + skfile), ishomo)
            skpassed = skpassed and curskfile.equals(refskfile, atol=ATOL,
                                                     rtol=RTOL)
            if skpassed:
                print('Passed.')
            else:
                print('Failed.')

        return skpassed


def infer_homonuclear(fname):
    '''Tries to infer whether the filename of an SK-file
       corresponds to a homo- or heteronuclear configuration.

    Args:

        fname (str): pathname of SK-file

    Returns:

        ishomo (bool): true, if SK-file seems to be homonuclear

    '''

    fname = fname.strip('_')
    fname = fname.strip('.skf')
    elements = fname.split('-')

    ishomo = elements[0] == elements[1]

    return ishomo


class TestWithWorkDirError(Exception):
    '''Exception thrown by the TestWithWorkDir class.'''


class SkfileError(Exception):
    '''Exception thrown by the Skfile class.'''
