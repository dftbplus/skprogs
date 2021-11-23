"""Contains the object needed to query a HSD-tree in a customized way.
"""

import sktools.hsd as hsd
from sktools.hsd.tree import Element

__all__ = ["HSDQueryError", "HSDMissingTagException", "HSDInvalidTagException",
           "HSDInvalidTagValueException", "HSDMissingAttributeException",
           "HSDInvalidAttributeException", "HSDInvalidAttributeValueException",
           "HSDQuery"]


class HSDQuery:
    """Class providing methods for querying a HSD-tree.

    Parameters
    ----------
    checkuniqueness : bool, optional
        Whether all query methods except `findchildren()` should check
        for the uniqueness of the child found.
    markprocessed : bool, optional
        Whether the nodes which have been queried should be marked as
        processed.
    """

    def __init__(self, checkuniqueness=False, markprocessed=False):
        self.chkunique = checkuniqueness
        self.mark = markprocessed


    def findchild(self, node, name, optional=False):
        """Finds a child of a node with a given name.

        Parameters
        ----------
        node : Element
            Parent node.
        name : string
            Name of the child to look for.
        optional : bool, optional
            Whether the child is optional only.

        Returns
        -------
         child : Element or None
            A hsd node if child has been found or None.

        Raises
        -------
        HSDMissingTagException
            If child was not found and the optional flag was False.
        HSDInvalidTagException
            Iff there are duplicates of the child and the query object was
            initialized with `check_uniqueness=True`.
        """
        if self.chkunique:
            children = node.findall("./" + name)
            if len(children) > 1:
                raise hsd.HSDInvalidTagException(
                    node=children[1],
                    msg="Double occurance of unique tag '{}'.".format(name))
            child = children[0] if children else None
        else:
            child = node.find("./" + name)
        if child is None and not optional:
            raise hsd.HSDMissingTagException(
                msg="Required tag '{}' not found.".format(name), node=node)
        self.markprocessed(child)
        return child


    def findchildren(self, node, name, optional=False):
        """Finds children of a node with given name.

        Parameters
        ----------
        node : Element
            Parent node.
        name : string
            Name of the children to look for.
        optional : bool, optional
            Whether the presence of at least one child is optional.

        Returns
        -------
        childlist : list
            List of child nodes or empty list.

        Raises
        ------
        HSDMissingTagException
            if no children were not found and the optional flag was False.
        """
        children = node.findall("./" + name)
        if not children and not optional:
            raise hsd.HSDMissingTagException(
                node=node,
                msg="No occurrence of required tag '{}' found.".format(
                    name))
        self.markprocessed(*children)
        return children


    def getchild(self, node, name, optional=False, defattribs=None):
        """Returns child with a given name.

        Parameters
        ----------
        node : Element
            Parent node.
        name : string
            Name of the child to look for.
        optional : bool, optional
            If set to `True`, an empty child node will be createad if a child
            with the given name does not exist.
        defattribs: dict, optional
            Default attribute dictionary for the child if it is created.
            Only makes sense if keyword argument `optional` was set to `True`.

        Returns
        -------
        child : Element
            The child with the given name. Either from the original HSD-tree
            or the one, which had been created. In latter case, the appropriate
            child is inserted into the tree.

        Raises
        ------
        HSDMissingTagException
            if the child was not found and keyword argument `optional`
            was not set to `True`.
        """
        child = self.findchild(node, name, optional)
        # findchild only returns if child has been found or optional is True.
        if child is None:
            child = self.setchild(node, name, defattribs)
        return child


    def getvalue(self, node, name, converter=None, defvalue=None,
                 attribs=None, defattribs=None, hsdblock=False,
                 returnchild=False):
        """Returns the converted value of the data stored in a child with a
        given name.

        Parameters
        ----------
        node : Element
            Parent node.
        name : string
            Name of the child to look for.
        converter : converter object
            Object with methods fromhsd() and tohsd() which can
            convert between the hsd element and the desired type. See
            converters in hsd.converter for examples.
        defvalue : arbitrary, optional
            Default value used if child has not been found. It will be
            converted to text by the tohsd() method of the converter.
        attribs : frozen set
            Set of attributes the node is allowed to have.
        defattribs: dict, optional
            Default attribute dictionary used if child has not been found.
        hsdblock : bool, optional
            Whether the given value should be added in hsd block notation
            (enclosed in curly braces) instead of an assignment.
        returnchild : bool, optional
            Whether not only the value but also the child node should be
            returned.

        Returns
        -------
        value : arbitrary
            The converted value of the child node's text or the default value
            if the child had not been found. In latter case, an appropriate
            node with the appropriate text representation of the default
            value is inserted into the tree.
        child : Element, optional
            Child node. Only returned, if `returnchild=True` was set.

        Raises
        ------
        HSDMissingTagException
            If child was not found and no default value had been specified.
        HSDInvalidTagValueException
            If conversion from tag values was unsuccessful.
        HSDInvalidAttributeException
            If node posses an attribute which is not allowed.

        Notes
        -----
        This method may store a reference to the converter object.
        Make sure you pass something which does not change afterwards.
        """
        optional = defvalue is not None
        child = self.findchild(node, name, optional)
        if child is not None:
            if len(child):
                raise hsd.HSDInvalidTagException("Unexpected children")
            self._checkattribs(child, attribs)
            if converter:
                try:
                    value = converter.fromhsd(child.text)
                except ValueError as ex:
                    raise hsd.HSDInvalidTagValueException(
                        msg="Conversion error: " + str(ex), node=child)
            else:
                value = child.text
            return (value, child) if returnchild else value
        else:
            child = self.setvalue(node, name, converter, defvalue, defattribs,
                                  hsdblock)
            return (defvalue, child) if returnchild else defvalue


    def getvaluenode(self, node, name, defvalue=None, defvaluename=None,
                     defattribs=None, hsdblock=False, allowtextvalue=False,
                     returnchild=False):
        """Returns the value node stored in a child with a given
        name. The child should contain either no nodes at all or only this one.

        Parameters
        ----------
        node : Element
            Parent node.
        name : string
            Name of the child to look for.
        defvalue : Element, optional
            If child is not found, it will be created and contain the specified
            node as subnode.
        defvaluename : string, optional
            If child is not found, it will be created and contain a subnode
            with the given name. If name is "", the child node will not
            contain a subnode. It is ignored, if defvalue had been specified.
        defattribs : dict, optional
            Default attribute dictionary used if child has not been found.
        hsdblock : bool, optional
            Whether the given value should be added in hsd block notation
            (enclosed in curly braces) instead of an assignment.
        allowtextvalue : bool, optional
            If set to yes, the child (if it exists) is allowed to have
            no subnode, but a text value instead. In that case the text
            value will be deleted and converted into an empty node.
        returnchild : bool, optional

        Returns
        -------
        node : Element or None
            The child node's first child or a node with the specified default
            name if child had not been found. In latter case, an
            appropriate child node with this as subnode is inserted into the
            tree.
        child : Element, optional
            The child not itself. Only returned if `returnchild=Yes` was set.

        Raises
        ------
        HSDMissingTagException
            If child was not found and no default value had been specified.
        HSDInvalidTagException
            If child has more than one child.

        Notes
        -----
        This routine should be used, if the child is not a leaf but
        contains a further child.
        """
        optional = defvalue is not None or defvaluename is not None
        child = self.findchild(node, name, optional)
        if child is not None:
            if len(child) > 1:
                raise hsd.HSDInvalidTagException(
                    "Tag '{}' is not allowed to have"
                    " more than one child".format(child.tag), node=child)
            self.markprocessed(child)
            if len(child):
                self.markprocessed(child[0])
                return (child[0], child) if returnchild else child[0]
            elif allowtextvalue and child.text:
                valuenode = Element(child.text)
                child.text = ""
                child.append(valuenode)
                self.markprocessed(valuenode)
                return (valuenode, child) if returnchild else valuenode
            else:
                return (None, child) if returnchild else None
        else:
            if defvalue is None:
                defvalue = Element(defvaluename)
            child, valuenode = self.setvaluenode(node, name, defvalue,
                                                 defattribs, hsdblock)
            return (valuenode, child) if returnchild else valuenode


    def setchild(self, node, name, attribs=None):
        """Creates an empty child with given name and attributes.

        Parameters
        ----------
        node : Element
            Parent node
        name : str
            Name of the child to create.
        attribs : dict
            Dictionary of attributes for the child.

        Returns
        -------
        child : Element
            Node which had been created and added to the tree.
        """
        child = Element(name, attribs or {})
        self.markprocessed(child)
        node.append(child)
        return child


    def setvalue(self, node, name, converter, value, attribs=None,
                 hsdblock=False):
        """Creates a child with the given value.

        Parameters
        ----------
        node : Element
            Parent node.
        name : str
            Name of the child node to create.
        converter : converter object
            Object with methods fromhsd() and tohsd() which can
            convert between the hsd element and the desired type. See
            converters in hsd.converter for examples.
        value : arbitrary
            Value which should be converted to text by the converter.
        attribs : dict
            Dictionary with attributes for the child.
        hsdblock : bool, optional
            Whether the given value should be added in hsd block notation
            (enclosed in curly braces) instead of an assignment.

        Returns
        -------
        child : Element
            The child node which had been created and added.
        """

        child = Element(name, attribs or {})
        if converter:
            child.text = converter.tohsd(value)
        else:
            child.text = value
        self.markprocessed(child)
        if not hsdblock:
            child.sethsd(hsd.HSDATTR_EQUAL, True)
        node.append(child)
        return child


    def setvaluenode(self, node, name, value=None, attribs=None,
                     hsdblock=False):
        """Creates a child with a node with a given name as only child.

        Parameters
        ----------
        node : Element
            Parent node.
        name : str
            Name of the child to create.
        value : str, optional
            Name of the node to create as child of the child node. If not
            specified, no subchild node is created.
        attribs : dict, optional
            Attributes of the child created.
        hsdblock : bool, optional
            Whether the given value should be added in hsd block notation
            (enclosed in curly braces) instead of an assignment.
        """
        child = Element(name, attribs or {})
        self.markprocessed(child)
        if value is not None:
            valuenode = Element(value)
            self.markprocessed(valuenode)
            child.append(valuenode)
        else:
            valuenode = None
        if not hsdblock:
            child.sethsd(hsd.HSDATTR_EQUAL, True)
        node.append(child)
        return child, valuenode





    def markprocessed(self, *nodes):
        """Marks nodes as having been processed, if the query object had been
        initialized with the appropriate option.

        Parameters
        ----------
        *nodes : list
            List of nodes to mark as processed.
        """
        if self.mark:
            for node in nodes:
                if node is not None:
                    node.sethsd(hsd.HSDATTR_PROC, True)


    def findunprocessednodes(self, node, allnodes=False):
        """Returns list of all nodes which had been not marked as processed.

        Parameters
        ----------
        node : Element
            Parent node.
        allnodes : bool, optional
            By default, only highest unprocessed nodes are returned, but
            not their children (which should be also unprocessed then). Setting
            `allnodes` to True, retuns all nodes.

        Returns
        -------
        nodelist : list of Elements
            List of all nodes, which have not been queried by a HSDQuery
            instance.
        """
        unprocessed = []
        for child in node:
            if child.gethsd(hsd.HSDATTR_PROC, None) is None:
                unprocessed.append(child)
                if not allnodes:
                    continue
            unprocessed += self.findunprocessednodes(child, allnodes)
        return unprocessed


    @staticmethod
    def _checkattribs(node, attribs):
        """Checks whether the node has only the allowed attributes

        Parameters
        ----------
        node : Element
            Node to investigate.
        attribs : frozen set
            Set of allowed attributes

        Raises
        ------
        HSDInvalidAttributeException
            If an invalid attribute is found.
        """
        nodekeys = frozenset(node.keys())
        if not nodekeys:
            return
        if not attribs:
            raise hsd.HSDInvalidAttributeException(
                node=node, msg="No attributes allowed.")
        if not attribs >= nodekeys:
            tmp = "', '".join(list(nodekeys - attribs))
            raise hsd.HSDInvalidAttributeException(
                node=node,
                msg="Tag '{}' contains invalid attribute(s) '{}'.".format(
                    node.tag, tmp))


def _test_module():
    """Testing module capabilities."""
    from io import StringIO
    from sktools.hsd.treebuilder import HSDTreeBuilder
    from sktools.hsd.parser import HSDParser
    from sktools.hsd.tree import HSDTree
    import sktools.hsd.converter as conv

    unit_attr = "unit"
    unit_only = frozenset([unit_attr])
    parser = HSDParser(defattrib=unit_attr)
    builder = HSDTreeBuilder(parser=parser)

    # Defining force type (scalar, list)
    force_units = {"ev/aa": 0.0194469050555}

    # Trivial unit conversion routine.
    def multiply_unit(child, value, unitattr, units):
        unit = child.get(unitattr)
        convfact = units.get(unit.lower(), None)
        if convfact is None:
            hsd.HSDInvalidAttributeValueException(
                node=child, msg="Invalid unit '{}'".format(unit))
        return value * convfact

    stream = StringIO("""
# Various driver possibilities
#                                  # No driver specification
#Driver {}                         # Use the default driver (whatever it is)
#Driver = None {}                  # Use the driver None {}
#Driver = None
Driver = ConjugateGradient {
    MaxForceComponent [eV/AA] = 1e-2
}

Hamiltonian = DFTB {
  # SCC = True
  # SCCTolerance = 1e-4
  # MaxSCCIterations = 100
  MaxAngularMomentum {
    O = "p"
    H = "s"
  }
  Mixer = Broyden
  #Mixer = Broyden {
  #  MixingParameter = 0.3
  #}
  #ReadInitialCharges = No
  KPointsAndWeights {
     0.0   0.0  0.0   0.25
     0.25 0.25 0.25   0.75
  }
}

Options {
  WriteAutotestTag = Yes
  UnknownOption = No
}

#ParserOptions {
#  ParserVersion = 4
#}
""")
    root = builder.build(stream)
    qy = HSDQuery(markprocessed=True)
    # A complex case: If driver was not specified, it defaults to None {}
    # If it was specified but nothing was assinged to it (no child)
    # it defaults to ConjugateGradient {}.
    dtype, driver = qy.getvaluenode(root, "Driver", "None",
                                    allowtextvalue=True, returnchild=True)
    # Since the in the previous getvaluenode() call a default had been specified
    # dtype can only be None, if "Driver" was in the input, but had no
    # child (e.g. 'Driver {}' or 'Driver = ;'). In this case we set
    # it to ConjugateGradient
    if dtype is None:
        dtype = qy.getchild(driver, "ConjugateGradient", optional=True)
    if dtype.tag == "None":
        pass
    elif dtype.tag == "ConjugateGradient":
        forcetol, child = qy.getvalue(
            dtype, "MaxForceComponent", conv.float0, 1e-4, returnchild=True,
            attribs=unit_only)
        multiply_unit(child, forcetol, unit_attr, force_units)
    elif dtype.tag == "SteepestDescent":
        forcetol, child = qy.getvalue(
            dtype, "MaxForceComponent", conv.float0, 1e-4, returnchild=True,
            attribs=unit_only)
        multiply_unit(child, forcetol, unit_attr, force_units)
        stepsize = qy.getvalue(dtype, "StepSize", conv.float0, 40.0)
        pass
    else:
        raise hsd.HSDInvalidTagException(
            node=dtype, msg="Unknown driver type '{}'".format(dtype.tag))

    ham = qy.getchild(root, "Hamiltonian")
    dftb = qy.getchild(ham, "DFTB")
    scc = qy.getvalue(dftb, "SCC", conv.bool0, True)
    scctol = qy.getvalue(dftb, "SCCTolerance", conv.float0, 1e-4)
    scciter = qy.getvalue(dftb, "MaxSCCIterations", conv.int0, 100)
    mangmom = qy.getchild(dftb, "MaxAngularMomentum")
    maxangs = [qy.getvalue(mangmom, species, conv.str0)
               for species in ["O", "H"]]
    mixer = qy.getvaluenode(dftb, "Mixer", "Broyden", allowtextvalue=True)
    if mixer.tag == "Broyden":
        mixparam = qy.getvalue(mixer, "MixingParameter", conv.float0, 0.2)
    else:
        raise hsd.HSDInvalidTagException(node=mixer,
                                        msg="Unknown mixer type '{}'.".format(
                                            mixer.tag))
    readcharges = qy.getvalue(dftb, "ReadInitalCharges", conv.bool0, False)
    kpoints = qy.getvalue(dftb, "KPointsAndWeights", conv.float1)
    if len(kpoints) % 4:
        raise hsd.HSDInvalidTagValueException(node=kpoints,
                                             msg="Incorrect number of floats")
    options = qy.getchild(root, "Options", optional=True)
    autotest = qy.getvalue(options, "WriteAutotestTag", conv.bool0, False)
    parseroptions = qy.getchild(root, "ParserOptions", optional=True)
    parserversion = qy.getvalue(parseroptions, "ParserVersion", conv.int0, 4)
    tree = HSDTree(root)
    tree.writehsd()
    print("\nUnprocessed: ", qy.findunprocessednodes(root))


if __name__ == "__main__":
    _test_module()
