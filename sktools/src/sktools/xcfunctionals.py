'''Defines supported xc-functionals.'''


import sktools.hsd as hsd
import sktools.hsd.converter as conv
import sktools.common as sc


class XCPBE0(sc.ClassDict):
    '''Globald PBE0 hybrid xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''

        myself = cls()
        myself.type = 'pbe0'
        return myself


class XCB3LYP(sc.ClassDict):
    '''Globald B3LYP hybrid xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''

        myself = cls()
        myself.type = 'b3lyp'
        return myself


class XCCAMB3LYP(sc.ClassDict):
    '''Range-separated CAMY-B3LYP xc-functional.

    Attributes
    ----------
    omega (float): range-separation parameter
    alpha (float): fraction of the global exact HF exchange
    beta (float): determines (alpha + beta) fraction of long-range HF exchange
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

        alpha, child = query.getvalue(root, 'alpha', conv.float0,
                                      returnchild=True)
        if not 0.0 <= alpha <= 1.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid alpha CAM-parameter {:f}'.format(alpha),
                node=child)

        beta, child = query.getvalue(root, 'beta', conv.float0,
                                     returnchild=True)
        if not 0.0 <= beta <= 1.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid beta CAM-parameter {:f}'.format(beta),
                node=child)

        if not 0.0 <= alpha + beta <= 1.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid CAM-parameter combination alpha={:f}, beta={:f}!\n'
                .format(alpha, beta) +
                'Should satisfy 0.0 <= alpha + beta <= 1.0', node=child)

        myself = cls()
        myself.type = 'cam-b3lyp'
        myself.omega = omega
        myself.alpha = alpha
        myself.beta = beta
        return myself


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
    'pbe0': XCPBE0,
    'b3lyp': XCB3LYP,
    'cam-b3lyp': XCCAMB3LYP
}
