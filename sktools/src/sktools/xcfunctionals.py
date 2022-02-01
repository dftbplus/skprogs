'''Defines supported xc-functionals.'''


import sktools.hsd as hsd
import sktools.hsd.converter as conv
import sktools.common as sc


class XCLCBNL(sc.ClassDict):
    '''Long-range corrected BNL xc-functional.

    Attributes
    ----------
    omega (float): range-separation parameter
    '''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        omega, child = query.getvalue(root, 'omega', conv.float0,
                                      returnchild=True)
        if omega <= 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid rs-parameter {:f}'.format(omega),
                node=child)

        myself = cls()
        myself.type = 'lc-bnl'
        myself.omega = omega
        return myself


class XCLCPBE(sc.ClassDict):
    '''Long-range corrected PBE xc-functional.

    Attributes
    ----------
    omega (float): range-separation parameter
    '''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        omega, child = query.getvalue(root, 'omega', conv.float0,
                                      returnchild=True)
        if omega <= 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid rs-parameter {:f}'.format(omega),
                node=child)

        myself = cls()
        myself.type = 'lc-pbe'
        myself.omega = omega
        return myself


class XCLocal(sc.ClassDict):
    '''Local xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = 'local'
        return myself


class XCPBE(sc.ClassDict):
    '''Semi-local PBE xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = 'pbe'
        return myself


class XCBLYP(sc.ClassDict):
    '''Semi-local BLYP xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = 'blyp'
        return myself


class XCLDA(sc.ClassDict):
    '''Local LDA xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = 'lda'
        return myself


# Registered xc-functionals with corresponding HSD name as key:
XCFUNCTIONALS = {
    'lc-bnl': XCLCBNL,
    'lc-pbe': XCLCPBE,
    'local': XCLocal,
    'pbe': XCPBE,
    'lda': XCLDA,
    'blyp': XCBLYP,
}
