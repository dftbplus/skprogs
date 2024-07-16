'''Contains various compression types.'''


import sktools.hsd as hsd
import sktools.hsd.converter as conv
import sktools.common as sc


#######################################################################
# Compressions
#######################################################################


class PowerCompression(sc.ClassDict):
    '''Compression by a power function (r/r0)^n.

    Attributes
    ----------
    power : float
        Power of the compression function (n).
    radius : float
        Radius of the compression (r0)
    '''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''

        power, child = query.getvalue(root, 'power', conv.float0,
                                      returnchild=True)
        if power < 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid compression power {:f}'.format(power), node=child)
        radius, child = query.getvalue(root, 'radius', conv.float0,
                                       returnchild=True)
        if radius <= 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid compression radius {:f}'.format(radius),
                node=child)

        myself = cls()
        myself.compid = SUPPORTED_COMPRESSIONS[
            myself.__class__.__name__.lower()]
        myself.power = power
        myself.radius = radius

        return myself


    def tohsd(self, root, query, parentname=None):
        ''''''

        if parentname is None:
            mynode = root
        else:
            mynode = query.setchild(root, 'PowerCompression')
        query.setchildvalue(mynode, 'power', conv.float0, self.power)
        query.setchildvalue(mynode, 'radius', conv.float0, self.radius)


    def __eq__(self, other):
        power_ok = abs(self.power - other.power) < 1e-8
        radius_ok = abs(self.radius - other.radius) < 1e-8
        return power_ok and radius_ok


class WoodsSaxonCompression(sc.ClassDict):
    '''Compression by the Woods Saxon potential.

    Attributes
    ----------
    ww : float
        Height (W) of the compression
    aa : float
        Slope (a) of the compression
    r0 : float
        Half-height radius of the compression
    '''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''

        ww, child = query.getvalue(
            root, 'w', conv.float0, defvalue=100.0, returnchild=True)
        if ww < 0.0:
            msg = 'Invalid potential height (W) {:f}'.format(ww)
            raise hsd.HSDInvalidTagValueException(msg=msg, node=child)

        aa, child = query.getvalue(root, 'a', conv.float0, returnchild=True)
        if aa <= 0.0:
            msg = 'Invalid potential slope (a) {:f}'.format(aa)
            raise hsd.HSDInvalidTagValueException(msg=msg, node=child)

        r0, child = query.getvalue(root, 'r0', conv.float0, returnchild=True)
        if r0 < 0.0:
            msg = 'Invalid potential half-height radius {:f}'.format(r0)
            raise hsd.HSDInvalidTagValueException(msg=msg, node=child)

        myself = cls()
        myself.compid = SUPPORTED_COMPRESSIONS[
            myself.__class__.__name__.lower()]
        myself.ww = ww
        myself.aa = aa
        myself.r0 = r0

        return myself


    def tohsd(self, root, query, parentname=None):
        ''''''

        if parentname is None:
            mynode = root
        else:
            mynode = query.setchild(root, 'WoodsSaxonCompression')

        query.setchildvalue(mynode, 'w', conv.float0, self.ww)
        query.setchildvalue(mynode, 'a', conv.float0, self.aa)
        query.setchildvalue(mynode, 'r0', conv.float0, self.r0)


    def __eq__(self, other):
        ww_ok = abs(self.ww - other.ww) < 1e-8
        aa_ok = abs(self.aa - other.aa) < 1e-8
        r0_ok = abs(self.r0 - other.r0) < 1e-8
        return ww_ok and aa_ok and r0_ok


# Registered compressions with corresponding hsd name as key
COMPRESSIONS = {
    'powercompression': PowerCompression,
    'woodssaxoncompression': WoodsSaxonCompression,
}
SUPPORTED_COMPRESSIONS = {
    'nocompression': 0,
    'powercompression': 1,
    'woodssaxoncompression': 2,
}


#######################################################################
# Compression containers
#######################################################################


class SingleAtomCompressions(sc.ClassDict):
    '''Compression container for cases where all compressed wavefunctions are
    determined from one single atomic calculation.

    Attributes
    ----------
    0,1,2.. : compression object
        Compression type for the given object.
    '''

    def getatomcompressions(self, atomconfig):
        '''Returns compressions for one or more atomic calculations.

        Parameters
        ----------
        atomconfig : AtomConfig
            Configuration of the atom, for which the compression container
            had been specified.

        Returns
        -------
        atomcompressions : list
            List of ( compressions, valenceshells ) tuples. Compressions
            is a list of compression objects with one compression for every
            angular momentum of the atom, representing a complete compression
            for an atomic calculation. Valencshells is a list of (nn, ll)
            tuples containing principal quantum number and angular momentum of
            the valenceshells, for which the wave function should be taken
            from that compressed calculation.
        '''
        compressions = []
        for ll in range(atomconfig.maxang + 1):
            if ll not in self:
                msg = 'Missing wave compression for shell {:s}'.format(
                    sc.ANGMOM_TO_SHELL[ll])
                raise sc.SkgenException(msg)
            compressions.append(self[ll])
        atomcompressions = [(compressions, atomconfig.valenceshells)]
        return atomcompressions


    @classmethod
    def fromhsd(cls, root, query):
        myself = cls()
        for ll, shellname in enumerate(sc.ANGMOM_TO_SHELL):
            child = query.findchild(root, shellname, optional=True)
            if child is None:
                break
            compr = sc.hsd_node_factory(
                'wavefunction compression', COMPRESSIONS,
                query.getvaluenode(child, '.'), query)
            myself[ll] = compr
        return myself


class MultipleAtomCompressions(sc.ClassDict):

    def getatomcompressions(self, atomconfig):
        '''Returns compressions for one or more atomic calculations.

        Parameters
        ----------
        atomconfig : AtomConfig
            Configuration of the atom, for which the compression container
            had been specified.

        Returns
        -------
        atomcompressions : list
            List of ( compressions, valenceshells ) tuples. Compressions
            is a list of compression objects with one compression for every
            angular momentum of the atom, representing a complete compression
            for an atomic calculation. Valencshells is a list of (nn, ll)
            tuples containing principal quantum number and angular momentum of
            the valenceshells, for which the wave function should be taken
            from that compressed calculation.
        '''
        atomcompressions = []
        for nn, ll in atomconfig.valenceshells:
            if (nn, ll) not in self:
                msg = 'Missing compression for shell {:d}{:s}'.format(
                    nn, sc.ANGMOM_TO_SHELL[ll])
                raise sc.SkgenException(msg)
            comprs = [self[(nn, ll)],] * (atomconfig.maxang + 1)
            atomcompressions.append((comprs, [(nn, ll), ]))
        return atomcompressions



    @classmethod
    def fromhsd(cls, root, query):
        myself = cls()
        for shellnode in root:
            try:
                nn, ll = sc.shell_name_to_ind(shellnode.tag)
            except ValueError:
                raise hsd.HSDInvalidTagException(
                    "Invalid shell name '{}'".format(shellnode.tag), shellnode)
            wavecompr = sc.hsd_node_factory(
                'wavefunction compression', COMPRESSIONS,
                query.getvaluenode(shellnode, '.'), query)
            myself[(nn, ll)] = wavecompr
        return myself


# Registered compression containers with corresponing hsd name as key
COMPRESSION_CONTAINERS = {
    'singleatomcompressions': SingleAtomCompressions,
    'multipleatomcompressions': MultipleAtomCompressions,
}
