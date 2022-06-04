import sktools.hsd as hsd
import sktools.hsd.parser as hsdparser
import sktools.hsd.tree as hsdtree


HSD_ROOT_TAG_NAME = "hsd"


class TreeBuilder:
    """Event driven HSD tree builder.

    Attributes:
        path: List of elements representing the path to the current element.
    """

    def __init__(self, element_factory=None, root=None):
        """Initializes the builder.

        Args:
            element_factory: Element factory to be used. If not specified
                hsdtree.Element is used. The interface of the element factory
                must be compatible with hsdtree.Element.
            root: Element representing the root of the tree. If not specified
                the builder automatically creates a root element.
        """
        if element_factory is None:
            self._element_factory = hsdtree.Element
        else:
            self._element_factory = element_factory
        if root is None:
            root = self._element_factory(HSD_ROOT_TAG_NAME)
        self.path = [ root, ]


    def start(self, tag, attribs, hsdattribs):
        """Opens a new element.

        Args:
            tag: Name of element.
            attribs: Dictionary with element attributes.
            hsdattribs: Dictionary with HSD attributes.

        Returns:
            The opened element.
        """
        if self.path[-1].text is not None:
            msg = "No children allowed, as element already contains data"
            raise hsd.HSDException(msg)
        elem = self._element_factory(tag, attribs, hsdattribs)
        self.path.append(elem)
        return elem


    def end(self, tag):
        """Closes the current element.

        Args:
            tag: Name of element.

        Returns:
            The closed element.
        """
        elem = self.path[-1]
        if elem.tag != tag:
            raise hsd.HSDException("Invalid closing tag name")
        del self.path[-1]
        self.path[-1].append(elem)
        return elem


    def data(self, data):
        """Adds text to the current element.

        Args:
            data: Text to be added.
        """
        elem = self.path[-1]
        if len(elem):
            msg = "No data allowed, as element already contains children"
            raise hsd.HSDException(msg)
        elem.text = data


    def close(self):
        """Closes the builder and returns the created tree.

        The builder should not be used after this call any more!

        Raises:
            HSDException: If there were some open elements.
        """
        if len(self.path) != 1:
            raise hsd.HSDException("Open elements remained")
        return self.path[0]



class VariableTreeBuilder:
    """Event driven HSD tree builder with variables."""

    SCOPE_TAG_NAME = "_scope_"
    VARDEF_TAG_PREFIX = "_var_"
    VARIABLE_NAME_PREFIX = "$"

    def __init__(self, element_factory=None, root=None):
        """Initializes the builder.

        Args:
            element_factory: Element factory to be used. If not specified
                hsdtree.Element is used. The interface of the element factory
                must be compatible with hsdtree.Element.
            root: Element representing the root of the tree. If not specified
                the builder automatically creates a top element.
        """
        builder = TreeBuilder(element_factory=element_factory, root=root)
        self._builders = [ builder, ]
        if element_factory is None:
            self._element_factory = hsdtree.Element
        else:
            self._element_factory = element_factory
        self._scopes = [ self._element_factory(self.SCOPE_TAG_NAME) ]


    def start(self, tag, attribs, hsdattribs):
        """Opens a new element.

        Args:
            tag: Name of element.
            attribs: Dictionary with element attributes.
            hsdattribs: Dictionary with HSD attributes.

        Returns:
            The opened element.
        """
        if self._is_variable(tag):
            if attribs:
                msg = "Variable defintion can not have any attributes"
                raise hsd.HSDInvalidAttributeException(msg)
            hsdattribs[hsd.HSDATTR_TAG] = tag
            elemname = self._convert_variable_to_tag_name(tag)
            variable = self._element_factory(elemname, {}, hsdattribs)
            variable_builder = TreeBuilder(root=variable)
            self._builders.append(variable_builder)
        else:
            builder = self._builders[-1]
            builder.start(tag, attribs, hsdattribs)
        self._scopes.append(self._element_factory(self.SCOPE_TAG_NAME))


    def end(self, tag):
        """Closes the current element.

        Args:
            tag: Name of element.

        Returns:
            The closed element.
        """
        del self._scopes[-1]
        if self._is_variable(tag):
            variable_builder = self._builders[-1]
            variable = variable_builder.close()
            del self._builders[-1]
            last_scope = self._scopes[-1]
            last_scope.append(variable)
        else:
            builder = self._builders[-1]
            builder.end(tag)


    def data(self, data):
        """Adds text to the current element.

        Args:
            data: Text to be added.
        """
        builder = self._builders[-1]
        elem = builder.path[-1]
        if self._is_variable(data):
            varname = self._convert_variable_to_tag_name(data.strip())
            vardef = self._lookup_variable(varname)
            if vardef is None:
                msg = "Undefined variable '" + data.strip() + "'"
                raise hsd.HSDException(msg)
            elem.text = vardef.text
            for child in vardef:
                elem.append(child)
        else:
            elem.text = data


    def close(self):
        """Closes the builder and returns the created tree.

        The builder should not be used after this call any more!

        Raises:
            HSDException: If there were some open elements.
        """
        if len(self._builders) > 1:
            raise hsd.HSDException("Unclosed variable defintion")
        builder = self._builders[0]
        elem = builder.close()
        return elem


    @property
    def path(self):
        """List of elements representing the path to the current element."""
        path = []
        for builder in self._builders:
            path += builder.path
        return path


    def _is_variable(self, name):
        """Checks whether given name is a variable."""
        return name.startswith(self.VARIABLE_NAME_PREFIX)


    def _lookup_variable(self, varname):
        """Looks up a variable in the currently opened scopes.

        Args:
            varname: Name of the variable.

        Returns:
            Variable element or None, if not found.
        """
        variable = None
        path_expression = "./" + varname
        for ind in range(len(self._scopes) - 1, -1, -1):
            scope = self._scopes[ind]
            variable = scope.find(path_expression)
            if variable is not None:
                break
        return variable


    def _convert_variable_to_tag_name(self, varname):
        """Converts a variable to a valid XML-tag"""
        tagname = self.VARDEF_TAG_PREFIX + varname[1:].lower()
        return tagname



class HSDTreeBuilder:
    """Builds HSD-tree by connecting parser with builder."""

    def __init__(self, parser=None, builder=None):
        """Initializes a HSDTreeBuilder instance.

        Args:
            parser: Event-driven HSD-parser (default: HSDParser)
            builder: Event-driven tree builder (default: TreeBuilder)
        """
        if parser:
            self.parser = parser
        else:
            self.parser = hsdparser.HSDParser()
        if builder:
            self.builder = builder
        else:
            self.builder = TreeBuilder()
        self.parser.start_handler = self.builder.start
        self.parser.close_handler = self.builder.end
        self.parser.text_handler = self.builder.data


    def build(self, fileobj):
        """Builds a HSD-tree from a file-like object.

        Args:
            fileobj: File like object containing an HSD in text form.

        Returns:
            HSD-tree
        """
        self.parser.feed(fileobj)
        return self.builder.close()



if __name__ == "__main__":
    from io import StringIO
    import sys
    stream = StringIO("""Geometry = GenFormat {
2  S
Ga As
1    1    0.00000000000E+00   0.00000000000E+00   0.00000000000E+00
2    2    0.13567730000E+01   0.13567730000E+01   0.13567730000E+01
0.00000000000E+00   0.00000000000E+00   0.00000000000E+00
0.27135460000E+01   0.27135460000E+01   0.00000000000E+00
0.00000000000E+00   0.27135460000E+01   0.27135460000E+01
0.27135460000E+01   0.00000000000E+00   0.27135460000E+01
}
Test[unit=1,
    dim=3]{}

$SlakoDef {
  SlaterKosterFiles [format=old] {
    $SpecValue = "As-As.skf"
    Ga-Ga = "./Ga-Ga.skf"
    Ga-As = "./Ga-As.skf"
    As-Ga = "./As-Ga.skf"
    As-As = $SpecValue
  }
}

Hamiltonian = DFTB {
  $MyTemp = 1000
  $Filling = Fermi {
    Temperature [Kelvin] = $mytemp
  }

  SCC = Yes
  SCCTolerance = 1.0E-007
  MaxSCCIterations = 1000
  Mixer = Broyden {}
  MaxAngularMomentum {
    Ga = "d"
    As = "p"
  }
  Filling = $Filling
  $SlakoDef
  KPointsAndWeights {
    0.0 0.0 0.0   1.0
  }
}

Options {
  AtomResolvedEnergies = No
  RestartFrequency = 20
  RandomSeed = 0
  WriteHS = No
}""")
    mybuilder = HSDTreeBuilder(parser=hsd.parser.HSDParser(),
                               builder=VariableTreeBuilder())
    hsdnodes = mybuilder.build(stream)
    tree = hsdtree.HSDTree(hsdnodes)
    tree.write(sys.stdout, encoding="unicode")
    tree.writehsd(sys.stdout)
