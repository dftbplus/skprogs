'''Defines supported xc-functionals.'''


import sktools.hsd as hsd
import sktools.hsd.converter as conv
import sktools.common as sc


class XCPBE0(sc.ClassDict):
    '''Globald PBE0 hybrid xc-functional.

    Attributes
    ----------
    alpha (float): fraction of the global exact HF exchange
    '''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''

        alpha, child = query.getvalue(root, 'alpha', conv.float0,
                                      returnchild=True)
        if not 0.0 <= alpha <= 1.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid alpha CAM-parameter {:f}'.format(alpha),
                node=child)

        if not 0.0 <= alpha <= 1.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid global HFX portion alpha={:f}!\n'
                .format(alpha) +
                'Should satisfy 0.0 <= alpha <= 1.0', node=child)

        myself = cls()
        myself.type = 'pbe0'
        # dummy omega
        myself.omega = 1.0
        myself.alpha = alpha
        # dummy beta
        myself.beta = 0.0
        return myself


class XCB3LYP(sc.ClassDict):
    '''Globald B3LYP hybrid xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''

        myself = cls()
        myself.type = 'b3lyp'
        # dummy omega
        myself.omega = 1.0
        myself.alpha = 0.2
        myself.beta = 0.0
        return myself


class XCCAMYB3LYP(sc.ClassDict):
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

        beta, child = query.getvalue(root, 'beta', conv.float0,
                                     returnchild=True)

        if not alpha + beta > 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid CAM-parameter combination alpha={:f}, beta={:f}!\n'
                .format(alpha, beta) +
                'Should satisfy alpha + beta > 0.0', node=child)

        myself = cls()
        myself.type = 'camy-b3lyp'
        myself.omega = omega
        myself.alpha = alpha
        myself.beta = beta
        return myself


class XCCAMYPBEH(sc.ClassDict):
    '''Range-separated CAMY-PBEh xc-functional.

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

        beta, child = query.getvalue(root, 'beta', conv.float0,
                                     returnchild=True)

        if not alpha + beta > 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid CAM-parameter combination alpha={:f}, beta={:f}!\n'
                .format(alpha, beta) +
                'Should satisfy alpha + beta > 0.0', node=child)

        myself = cls()
        myself.type = 'camy-pbeh'
        myself.omega = omega
        myself.alpha = alpha
        myself.beta = beta
        return myself


class XCLCYBNL(sc.ClassDict):
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
        myself.type = 'lcy-bnl'
        myself.omega = omega
        myself.alpha = 0.0
        myself.beta = 1.0
        return myself


class XCLCYPBE(sc.ClassDict):
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
        myself.type = 'lcy-pbe'
        myself.omega = omega
        myself.alpha = 0.0
        myself.beta = 1.0
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
    'lcy-bnl': XCLCYBNL,
    'lcy-pbe': XCLCYPBE,
    'local': XCLocal,
    'pbe': XCPBE,
    'lda': XCLDA,
    'blyp': XCBLYP,
    'pbe0': XCPBE0,
    'b3lyp': XCB3LYP,
    'camy-b3lyp': XCCAMYB3LYP,
    'camy-pbeh': XCCAMYPBEH
}
