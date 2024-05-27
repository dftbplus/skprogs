"""Parser for the skdefs.hsd file."""

import re
import copy
import numpy as np
import sktools.hsd as hsd
import sktools.hsd.converter as conv
from sktools.hsd.treebuilder import HSDTreeBuilder, VariableTreeBuilder
from sktools.hsd.query import HSDQuery
from sktools.hsd.parser import HSDParser
from . import common as sc
from . import compressions
from . import twocenter_grids
from . import calculators
from . import xcfunctionals


CURRENT_SKDEF_VERSION = 1
ENABLED_SKDEF_VERSIONS = frozenset([CURRENT_SKDEF_VERSION])


class Skdef(sc.ClassDict):
    """Represents the full input file 'skdef.hsd'.

    Attributes
    ----------
    globals : Globals
        Global settings
    atomparameters : AtomParameters
        Various atomic parameters
    oncenterparameters : OnecenterParameters
        Parameters influencing the technical details of the one-center
        calculation.
    twocenterparameters : TwocenterParameters
        Parameters influencing the technical details of the two-center
        calculation.
    """
    @classmethod
    def fromhsd(cls, root, query):
        myself = cls()
        version = query.getvalue(root, "skdefversion", conv.int0)
        cls._check_version(version)
        node = query.getchild(root, "globals")
        myself.globals = Globals.fromhsd(node, query)
        node = query.getchild(root, "atomparameters")
        myself.atomparameters = AtomParameters.fromhsd(node, query)
        node = query.getchild(root, "onecenterparameters")
        myself.onecenterparameters = OnecenterParameters.fromhsd(node, query)
        node = query.getchild(root, "twocenterparameters")
        myself.twocenterparameters = TwocenterParameters.fromhsd(node, query)
        return myself

    @classmethod
    def fromfile(cls, fileobj):
        parser = HSDParser()
        builder = VariableTreeBuilder()
        treebuilder = HSDTreeBuilder(parser=parser, builder=builder)
        openclose = isinstance(fileobj, str)
        if openclose:
            fp = open(fileobj, "r")
        else:
            fp = fileobj
        tree = treebuilder.build(fp)
        if openclose:
            fp.close()
        query = HSDQuery(checkuniqueness=True, markprocessed=True)
        return cls.fromhsd(tree, query)


    def update(self, other):
        """Extends data with the data in an other skdefs.

        Parameters
        ----------
        other : Skdef
            Data to use for extending.
        """
        if other.globals != self.globals:
            raise sc.SkgenException(
                "Incompatible globals, skdefs can not be merged.")
        self.atomparameters.update(other.atomparameters)


    @staticmethod
    def _check_version(version):
        if version not in ENABLED_SKDEF_VERSIONS:
            msg = "Invalid skdef version {:d}".format(version)
            raise sc.SkgenException(msg)



class Globals(sc.ClassDict):
    """Global settings.

    Attributes
    ----------
    functional : int
        DFT functional (sktools.common.FUNCTIONAL_{LDA,PBE}).
    radius : int
        Superposition type (sktools.common.SUPERPOSITION_{POTENTIAL,DENSITY})
    """
    @classmethod
    def fromhsd(cls, root, query):
        """Creates instance from a HSD-node and with given query object."""

        superpos, child = query.getvalue(root, "superposition", conv.str0,
                                         returnchild=True)
        if superpos not in sc.SUPERPOSITION_TYPES:
            raise hsd.HSDInvalidTagValueException(
                "Invalid superposition type '{}'".format(superpos), child)

        # read the functional
        xcf = sc.hsd_node_factory('xc', xcfunctionals.XCFUNCTIONALS,
                                  query.getvaluenode(root, 'xcfunctional'),
                                  query)
        if xcf.__class__ not in xcfunctionals.XCFUNCTIONALS.values():
            raise hsd.HSDInvalidTagValueException(
                "Invalid functional type '{}'".format(xcf), child)

        myself = cls()
        myself.superposition = sc.SUPERPOSITION_TYPES[superpos]
        myself.xcf = xcf

        return myself


class AtomParameters(sc.ClassDict):
    """Atomic parameters

    Attributes
    ----------
    "elementname" : ClassDict
        ClassDict with fields `atomconfig` (type `AtomConfig`) and `dftbatom`
        (type `DftbAtom`) fields containing those settings for the given atom.
    """
    @classmethod
    def fromhsd(cls, root, query):
        myself = cls()
        for elemnode in query.findchildren(root, "*"):
            try:
                atomparam = sc.ClassDict()
                node = query.getchild(elemnode, "atomconfig")
                atomparam.atomconfig = AtomConfig.fromhsd(node, query)
                node = query.getchild(elemnode, "dftbatom")
                atomparam.dftbatom = DftbAtom.fromhsd(node, query)
            except sc.SkgenException as ex:
                msg = "AtomParameters/{}:\n{}".format(
                    elemnode.tag, ex)
                raise sc.SkgenException(msg)
            myself[elemnode.tag] = atomparam
        return myself


class AtomConfig(sc.ClassDict):
    """Represents the configuration of a free atom.

    Attributes
    ----------
    znuc : float
        Nuclear charge.
    mass : float
        Mass of the atom.
    occupations : list
        Either spin polarized (by default) or spin averaged occupation.
        It can be changed via the `make_spinpolarized` `make_spinaveraged`
        methods.
    occupations_spinpol : list
        List of (nup, ndown) tuples for each shell (e.g.
        [[ (1.0, 1.0), (1.0, 1.0) ], [ (3.0, 2.0), ]] for N)
    occupations_spinavg : list
        Same as occupations but averaged out for spin up and spin down.
    valenceshells : list
        List of (n, l) tuples representing the valence shells.
    relativistics : int
        Type of relativistics. (None, "zora")
    maxang : int
        Maximal angular momentum.
    nelec : float
        Number of electrons.
    spinpolarized : bool
        Whether atom is spinpolarized.
    charge : float
        Charge of the atom.
    charged : bool
        Whether atom has a net charge.
    """

    # Tolerance for treating electron populations being equal
    _ELECTRON_TOL = 1e-8


    def __init__(self, atomicnumber, mass, occupations, valenceshells,
                 occshells, relativistics, charge=0.0):
        super().__init__()
        self.atomicnumber = atomicnumber
        self.mass = mass

        # Sort valenceshells (and occupations) by ascending nn and ll
        tmp = [nn * (sc.MAX_ANGMOM + 1) + ll for nn, ll in valenceshells]
        self.valenceshells = [valenceshells[ii] for ii in np.argsort(tmp)]

        # check for uniqueness and continuity of angular quantum numbers
        angmom = sorted([tpl[1] for tpl in self.valenceshells])
        unique_and_continuous = all(ii + 1 == jj
                                    for ii, jj in zip(angmom, angmom[1:]))
        if not unique_and_continuous:
            shell_str = ' '.join([sc.shell_ind_to_name(nn, ll)
                                  for nn, ll in self.valenceshells])
            raise sc.SkgenException(
                "Invalid valence shell configuration '" + shell_str \
                + "' found:\nDuplicate angular momenta and/or omitting " + \
                "intermediate shells is not supported by the SK-file format.")

        # Sort occshells by ascending nn and ll
        tmp = [qn[0] * (sc.MAX_ANGMOM + 1) + qn[1] for qn, occ in occshells]
        self.occshells = [occshells[ii] for ii in np.argsort(tmp)]

        self.occupations_spinpol = occupations
        self.occupations = self.occupations_spinpol

        self.relativistics = sc.RELATIVISTICS_TYPES.get(relativistics, None)
        if self.relativistics is None:
            raise sc.SkgenException(
                "Invalid relativistics type '{}'".format(relativistics))

        # If any valenceshell has higher n or l as occupations are listed for,
        # fill up occupations with zeros accordingly
        maxl = 0
        maxn = [0,] * (sc.MAX_ANGMOM + 1)
        for nn, ll in valenceshells:
            maxl = max(ll, maxl)
            maxn[ll] = max(nn, maxn[ll])
        if maxl > len(self.occupations) - 1:
            self.occupations += [[],] * (maxl - len(self.occupations) + 1)
        for ll, occ_l in enumerate(occupations):
            # At least one occupation for each angular momentum up to lmax.
            if not len(occ_l):
                occ_l.append((0.0, 0.0))
            # Extend occupations up to highest principal quantum number in
            # valence shells
            if maxn[ll] - ll > len(occ_l):
                occ_l.extend([(0.0, 0.0)] * (maxn[ll] - ll - len(occ_l)))
        self.maxang = len(self.occupations) - 1

        self.nelec = 0.0
        self.spinpolarized = False
        self.occupations_spinavg = []
        for shellocc in self.occupations:
            occ_l = []
            for nup, ndown in shellocc:
                nn = nup + ndown
                self.nelec += nn
                occ_l.append(( nn / 2.0, nn / 2.0))
                self.spinpolarized = (self.spinpolarized
                                      or abs(nup - ndown) > self._ELECTRON_TOL)
            self.occupations_spinavg.append(occ_l)
        self.charge = self.atomicnumber - self.nelec
        if abs(self.charge - charge) > self._ELECTRON_TOL:
            msg = "Mismatch between specified total charge and occupations " \
                "({:.8f} vs. {:.8f})".format(charge, self.charge)
            raise sc.SkgenException(msg)
        self.charged = abs(self.charge > self._ELECTRON_TOL)


    def make_spinpolarized(self):
        """Make sure `occupation` attribute represents spin polarized state."""
        self.occupations = copy.deepcopy(self.occupations_spinpol)


    def make_spinaveraged(self):
        """Make sure `occupation` attribute represents spin averaged state."""
        self.occupations = copy.deepcopy(self.occupations_spinavg)


    @classmethod
    def fromhsd(cls, root, query):
        """Initializes an AtomProperties object from a HSD-tree.

        Parameters
        ----------
        root : HSDTree instance
            Root of the node containing the information.
        query : HSDQuery instance
            Object used for querying the tree.

        Returns
        -------
        atomconfig : AtomProperties
            Initialized Atomconfig instance.
        """
        znuc, child = query.getvalue(root, "atomicnumber", conv.float0,
                                     returnchild=True)
        if znuc < 0.0 or znuc > 95.0:
            raise hsd.HSDInvalidTagValueException(
                msg="Invalid nuclear charge {:f}".format(znuc),
                node=child)
        mass, child = query.getvalue(root, "mass", conv.float0,
                                     returnchild=True)
        if mass < 0 or mass > 250.0:
            raise hsd.HSDInvalidTagValueException(
                msg="Invalid atomic mass {:f}".format(mass), node=child)

        occshellnames = []
        occupations = []
        occnode = query.findchild(root, "occupations")
        for ll, shellname in enumerate(sc.ANGMOM_TO_SHELL):
            occ_l = []
            for nn in range(ll + 1, sc.MAX_PRINCIPAL_QN + 1):
                txt = "{:d}{:s}".format(nn, shellname)
                shelloccnode = query.findchild(occnode, txt, optional=True)
                if shelloccnode is None:
                    break
                tmp = query.getvalue(shelloccnode, ".", conv.float1)
                if len(tmp) != 2:
                    raise hsd.HSDInvalidTagValueException(
                        msg="Invalid number of occupation numbers",
                        node=shelloccnode)
                occ_l.append((tmp[0], tmp[1]))
                occshellnames.append((txt, (tmp[0] + tmp[1])))
            if len(occ_l):
                occupations.append(occ_l)

        valshellnames, child = query.getvalue(root, "valenceshells",
                                              conv.str1, returnchild=True)
        valshells = []
        for valshellname in valshellnames:
            try:
                valshell = sc.shell_name_to_ind(valshellname)
                valshells.append(valshell)
            except ValueError:
                raise hsd.HSDInvalidTagValueException(
                    msg="Invalid shell name '{}'".format(valshellname),
                    node=child)

        occshells = []
        for occshellname, occ in occshellnames:
            occshell = sc.shell_name_to_ind(occshellname)
            occshells.append((occshell, occ))

        relattype, child = query.getvalue(root, "relativistics", conv.str0,
                                          "none", returnchild=True)
        relattype = relattype.lower()
        if relattype not in sc.RELATIVISTICS_TYPES:
            raise hsd.HSDInvalidTagValueException(
                msg="Invalid relativistics type '{}'".format(relattype))

        return cls(znuc, mass, occupations, valshells, occshells, relattype)


    def __eq__(self, other):
        if not isinstance(other, AtomConfig):
            return False
        if (abs(self.atomicnumber - other.atomicnumber)
                > sc.INPUT_FLOAT_TOLERANCE):
            return False
        if abs(self.mass - other.mass) > sc.INPUT_FLOAT_TOLERANCE:
            return False
        if len(self.occupations_spinpol) != len(other.occupations_spinpol):
            return False
        for occ_l1, occ_l2 in zip(self.occupations_spinpol,
                                  other.occupations_spinpol):
            if len(occ_l1) != len(occ_l2):
                return False
            occ1 = np.array(occ_l1)
            occ2 = np.array(occ_l2)
            if np.any(np.abs(occ1 - occ2) > sc.INPUT_FLOAT_TOLERANCE):
                return False
        if self.valenceshells != other.valenceshells:
            return False
        if self.relativistics != other.relativistics:
            return False
        return True



class DftbAtom(sc.ClassDict):
    """Contains settings related to atoms in DFTB.

    Attributes
    ----------
        shellresolved :  bool
            Whether shell resolved Hubbard U values should be used.
        customizedonsites : dict
            (n, l) indexed dictionary of onsite values which should be
            overriden.
        customizedhubbards : dict
            (n, l) indexed dictionary with override values for the
            Hubbard parameter.
        customizedoccupations: dict
            (n, l) indexed dictionary with override values for the
            occupations.
        densitycompression : compression object
            Contains the details how density should be compressed.
        wavecompressions : compression objects
            Contains the type of compressions for the wavefunction.
    """

    @classmethod
    def fromhsd(cls, root, query):
        """Creates instance from a HSD-node and with given query object."""

        shellresolved = query.getvalue(root, "shellresolved", conv.bool0)

        customonsites_node = query.getchild(root, "customizedonsites",
                                            optional=True)
        customonsites = sc.get_shellvalues(customonsites_node, query)

        customhubbards_node = query.getchild(root, "customizedhubbards",
                                             optional=True)
        customhubbards = sc.get_shellvalues(customhubbards_node, query)

        customoccupations_node = query.getchild(root, "customizedoccupations",
                                                optional=True)
        customoccupations = sc.get_shellvalues(customoccupations_node, query)

        denscompr = sc.hsd_node_factory(
            "density compression", compressions.COMPRESSIONS,
            query.getvaluenode(root, "densitycompression"), query)
        wavecomprs = sc.hsd_node_factory(
            "wave compression container", compressions.COMPRESSION_CONTAINERS,
            query.getvaluenode(root, "wavecompressions"), query)

        myself = cls()
        myself.shellresolved = shellresolved
        myself.densitycompression = denscompr
        myself.wavecompressions = wavecomprs
        myself.customizedonsites = customonsites
        myself.customizedhubbards = customhubbards
        myself.customizedoccupations = customoccupations

        return myself


class OnecenterParameters(sc.ClassDict):
    """One center parameters with defaults.

    Attributes
    ----------
    elementname : ClassDict
        Contains one center settings in the fields `deltafilling` and
        `calculator`.
    """

    @classmethod
    def fromhsd(cls, root, query):
        """Returns one center parameters with substituted defaults."""
        myself = cls()

        # Parse all other nodes
        for node in query.findchildren(root, "*"):
            name = node.tag
            try:
                myself[name] = OnecenterParameter.fromhsd(node, query)
            except sc.SkgenException as ex:
                msg = "onecenterparameters/{}:\n{}".format(name, ex)
        return myself


class OnecenterParameter(sc.ClassDict):

    @classmethod
    def fromhsd(cls, root, query):
        myself = cls()
        myself.deltafilling = query.getvalue(root, "deltafilling", conv.float0)
        myself.calculator = sc.hsd_node_factory(
            "one-center calculator",
            calculators.ONECENTER_CALCULATOR_SETTINGS,
            query.getvaluenode(root, "calculator"), query)
        return myself

    def __eq__(self, other):
        if not isinstance(other, OnecenterParameter):
            return False
        if (abs(self.deltafilling - other.deltafilling)
                > sc.INPUT_FLOAT_TOLERANCE):
            return False
        if self.calculator != other.calculator:
            return False
        return True



class TwocenterParameters(sc.ClassDict):
    """Two center parameters with defaults.

    Attributes
    ----------
    elementname : ClassDict
        Contains two center settings in fields `grid` and
        `calculator`.
    """

    _PATTERN_DEFAULT = re.compile(
        r"^([a-z][a-z0-9_]*)-([a-z][a-z0-9_]*)$", re.IGNORECASE)

    @classmethod
    def fromhsd(cls, root, query):
        """Returns two center parameters with substituted defaults."""
        myself = cls()

        # Parse all other nodes
        for node in query.findchildren(root, "*"):
            name = node.tag
            match = cls._PATTERN_DEFAULT.match(name)
            if not match:
                msg = "Invalid two center interaction '{}'".format(name)
                raise sc.SkgenException(msg)
            name1, name2 = match.groups()
            key = min(name1, name2), max(name1, name2)
            try:
                myself[key] = TwocenterParameter.fromhsd(node, query)
            except sc.SkgenException as ex:
                msg = "twocenterparameters/{}-{}:\n{}".format(name1, name2, ex)
                raise sc.SkgenException(msg)
        return myself


class TwocenterParameter(sc.ClassDict):

    @classmethod
    def fromhsd(cls, root, query):
        myself = cls()
        myself.grid = sc.hsd_node_factory(
            "two-center grid", twocenter_grids.TWOCENTER_GRIDS,
            query.getvaluenode(root, "grid"), query)
        myself.calculator = sc.hsd_node_factory(
            "two-center calculator",
            calculators.TWOCENTER_CALCULATOR_SETTINGS,
            query.getvaluenode(root, "calculator"), query)
        return myself


def _test_module():
    from sktools.hsd.treebuilder import HSDTreeBuilder
    from sktools.hsd.query import HSDQuery
    from sktools.hsd.parser import HSDParser

    parser = HSDParser(lowertagnames=True)
    treebuilder = HSDTreeBuilder(parser=parser)
    fp = open("skdefs.hsd", "r")
    tree = treebuilder.build(fp)
    fp.close()
    query = HSDQuery(checkuniqueness=True, markprocessed=True)
    skdefs = Skdef.fromhsd(tree, query)
    print(skdefs.onecenterparameters["n"].calculator.exponents)
    print(skdefs.onecenterparameters["n"].deltafilling)

    unprocessed_list = query.findunprocessednodes(tree)
    for unprocessed in unprocessed_list:
        print("Unprocessed element '{}' at line {:d}!".format(
            unprocessed.hsdattrib[hsd.HSDATTR_TAG],
            unprocessed.hsdattrib[hsd.HSDATTR_LINE] + 1))


if __name__ == "__main__":
    _test_module()
