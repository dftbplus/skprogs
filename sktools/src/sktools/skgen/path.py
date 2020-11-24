import sktools.common as sc
import re
from sktools.hsd.tree import Element, SubElement
from sktools.hsd.query import HSDQuery
from sktools.hsd.parser import HSDParser
from sktools.hsd.treebuilder import VariableTreeBuilder, HSDTreeBuilder


_PATTERN_ONECENTER_TAG = re.compile(r"^([a-z:]+)$", re.IGNORECASE)
_PATTERN_TWOCENTER_TAG = re.compile(r"^([a-z:]+)-([a-z:]+)$", re.IGNORECASE)


class SkgenPaths:
    """Stores working paths used by skgen."""

    def __init__(self, root=None, query=None):
        """Initializes an SkgenPaths instance.

        Args:
            root: Root of the hsd tree storing the paths (default=None).
            query: Query object to use to query the path tree (default=None).
        """
        if root is None:
            self._root = Element("hsd")
        else:
            self._root = root
        if query is None:
            self._query = HSDQuery(checkuniqueness=True)
        else:
            self._query = query
        self._onecenter_nodes = {}
        self._twocenter_nodes = {}
        self._store_nodes()


    @classmethod
    def fromhsd(cls, root, query):
        """Initializes an SkgenPaths instance from an existing HSD tree.

        Args:
            root: Root of the hsd tree storing the paths.
            query: Query object to use to query the path tree.

        Returns:
            Initialized instance.
        """
        myself = cls(root, query)
        return myself


    @classmethod
    def fromfile(cls, fileobj):
        """Initializes an SkgenPaths instance from a file.

        Args:
            fileobj: File name or file like object containing the text
                representation of a path tree.

        Returns:
            Initialized instance.
        """
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
        query = HSDQuery(checkuniqueness=True)
        myself = cls(tree, query)
        return myself


    def get_onecenter_workdir(self, elem, calctype, default):
        """Delivers working directory for a onecenter calculation.

        Args:
            elem: Element to process (must be lower case!)
            calctype: Type of the calculation (e.g. 'atom', 'potcomp', ...)
            default: Directory to return (and store), if no directory was
                found in the path tree yet.

        Returns:
            Working directory for the given calculation type.
        """
        elemnode = self._onecenter_nodes.get(elem)
        if elemnode is None:
            elemnode = SubElement(self._root, elem)
            self._onecenter_nodes[elem] = elemnode
        workdir = self._query.getvalue(elemnode, calctype, defvalue=default)
        return workdir


    def get_twocenter_workdir(self, elem1, elem2, calctype, default):
        """Delivers working directory for a two-center calculation.

        Args:
            elem1: First element (must be lower case!)
            elem2: Second element (must be lower case!)
            calctype: Type of the calculation (e.g. 'atom', 'potcomp', ...)
            default: Directory to return (and store), if no directory was
                found in the path tree yet.

        Returns:
            Working directory for the given calculation type.
        """
        elem1, elem2 = min(elem1, elem2), max(elem1, elem2)
        elemnode = self._twocenter_nodes.get(( elem1, elem2 ))
        if elemnode is None:
            name = elem1 + "-" + elem2
            elemnode = SubElement(self._root, name)
            self._twocenter_nodes[elem1, elem2] = elemnode
        workdir = self._query.getvalue(elemnode, calctype, defvalue=default)
        return workdir


    def get_paths(self):
        """Returns an hsd-tree with the stored paths.
        """
        return self._root


    def _store_nodes(self):
        """Sort out nodes into one-center and two-center ones.
        """
        for child in self._query.findchildren(self._root, "*"):
            name = child.tag
            match = _PATTERN_TWOCENTER_TAG.match(name)
            if match:
                elem1, elem2 = match.groups()
                elem1, elem2 = ( min(elem1, elem2), max(elem1, elem2) )
                if ( elem1, elem2 ) in self._twocenter_nodes:
                    msg = "Multiple two-center defintions for {}-{}".format(
                        elem1, elem2)
                    raise sc.SkgenException(msg)
                self._twocenter_nodes[elem1, elem2] = child
            else:
                match = _PATTERN_ONECENTER_TAG.match(name)
                if match:
                    elem = match.groups(0)
                    if elem in self._onecenter_nodes:
                        msg = "Multiple one-center defintions for {}".format(
                            elem)
                        raise sc.SkgenException(msg)
                    self._onecenter_nodes[elem] = child
                else:
                    msg = "Invalid node name '{}'".format(name)
                    raise sc.SkgenException(msg)
