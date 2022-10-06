#!/usr/bin/env python3

'''Collects coefficient information for waveplot.'''


import os.path

from sktools.common import ANGMOM_TO_SHELL, writefloats
from sktools.taggedfile import TaggedFile
from sktools.skdef import Skdef
from sktools.oldskfile import OldSKFile


USAGE = \
    '''Collects coefficient information for waveplot. It iterates over the
    elements defined in skdefs.py and collects the wave function coefficients
    and other information necessary for waveplot. The homonuclear SK-files for
    those elements must have been created already. If it is missing, the given
    element will be ignored.
    '''


def writecoeffs(fp, elem, atomconfig, homoskname, wavecompdir):
    '''Writes element-specific input, processed by Waveplot.

    Args:

        fp (file object): file object to write to
        elem (str): element name to fetch information for
        atomconfig (AtomConfig): represents the configuration of a free atom
        homoskname (str): pathname of homonuclear Slater-Koster file
        wavecompdir (str): path to calculation of the compressed atom

    '''

    homosk = OldSKFile.fromfile(homoskname, True)
    cutoff = homosk.nr * homosk.dr / 2.0
    fp.write('{} {{\n'.format(elem))
    fp.write('  AtomicNumber = {:d}\n'.format(atomconfig.znuc))
    for nn, ll in atomconfig.valenceorbs:
        coeffsname = 'coeffs_{:02d}{:1s}.tag'.format(nn, ANGMOM_TO_SHELL[ll])
        coeffs = TaggedFile.fromfile(os.path.join(wavecompdir, coeffsname),
                                     transpose=True)
        fp.write('  Orbital {\n')
        fp.write('    AngularMomentum = {:d}\n'.format(ll))
        fp.write('    Occupation = {:.1f}\n'.format(coeffs['occupation']))
        fp.write('    Cutoff = {:5.2f}\n'.format(cutoff))
        fp.write('    Exponents {\n')
        writefloats(fp, coeffs['exponents'], indent=6, numperline=3,
                    formstr='{:21.12E}')
        fp.write('    }\n')
        fp.write('    Coefficients {\n')
        writefloats(fp, coeffs['coefficients'], indent=3, numperline=3,
                    formstr='{:21.12E}')
        fp.write('    }\n')
        fp.write('  }\n')
    fp.write('}\n')


def main():
    '''Main driver routine.'''

    skdefs = Skdef.fromfile('skdefs.py')
    atomconfigs = skdefs.atomconfigs
    elems = atomconfigs.keys()

    with open('wfc.hsd', 'w', encoding='utf8') as fp:

        for elem in elems:
            homoskname = '{elem}-{elem}.skf'.format(elem=elem)
            wavecompdir = os.path.join(elem, 'wavecomp')
            filespresent = (os.path.exists(homoskname)
                            and os.path.exists(wavecompdir))
            if not filespresent:
                print('*** Skipping: ', elem)
                continue

            print('*** Processing: ', elem)
            atomconfig = atomconfigs[elem]
            writecoeffs(fp, elem, atomconfig, homoskname, wavecompdir)


if __name__ == '__main__':
    main()
