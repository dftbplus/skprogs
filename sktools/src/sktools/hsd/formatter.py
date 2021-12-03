"""Formatting utilities for HSD content.
"""
import sys
import sktools.hsd as hsd

__all__ = [ "HSDFormatter", "HSDStreamFormatter" ]


class HSDFormatter:
    """Event controlled formatter producing HSD output."""

    def __init__(self, indentstring="  ", closecomments=False, defattrib=None):
        """Initializes HSDFormatter instance.

        Args:
            indentstring: String used for indenting (default: "  ").
            closecomments: Whether comments after tag closing should indicate
                which tag was closed (default: False).
            defattrib: When specified, attribute with that name is handled as
                default. When it is the only attribute, the name is not printed
                just the value. (default: None)

        Note:
            Per default, formatter writes to stdout. You can override this
            by calling its set_output() method.
        """
        self.output = sys.stdout
        self._closecomments = closecomments
        self._indent = indentstring
        self._defattrib = defattrib
        self._firsttag = True
        self._curindent = ""
        self._indentlist = []
        self._equalsigns = [ False, ]
        self._last2 = self._last = 0


    def set_output(self, output):
        """Sets the output for the formatter.

        Args:
            output: Open file or file object, to write the formatted output in.
        """
        self.output = output


    def start_tag(self, tagname, options, hsdoptions):
        """Starts an HSD tag.

        Args:
            tagname: Name of the tag to be started.
            options: Dictionary of the tag options.
        """
        tagname = hsdoptions.get(hsd.HSDATTR_TAG, tagname)
        equalsign = hsdoptions.get(hsd.HSDATTR_EQUAL, False)  # opens with '='?
        if options:
            if (self._defattrib and len(options) == 1
                    and self._defattrib in options):
                optstr = " [" + options[self._defattrib] + "]"
            else:
                optlist = [ key + "=" + value
                            for key, value in options.items() ]
                optstr = " [" + ",".join(optlist) + "]"
        else:
            optstr = ""
        if self._firsttag:
            indent = self._curindent
            self._firsttag = False
        else:
            indent = "" if self._equalsigns[-1] else "\n" + self._curindent
        trailing = " = " if equalsign else " {"
        self.output.write(indent + tagname + optstr + trailing)
        self._equalsigns.append(equalsign)
        self._increaseindentation()
        self._last2, self._last = self._last, 1


    def close_tag(self, tagname):
        """Closes an HSD tag.

        Args:
            tagname: Name of the tag to be closed.
        """
        self._decreaseindentation()
        if not self._equalsigns[-1]:
            if self._last == 1:
                self.output.write("}")
            else:
                self.output.write("\n" + self._curindent + "}")
                if self._closecomments:
                    self.output.write(" # " + tagname)
        elif self._closecomments and self._last == 2 and self._last2 != 1:
            self.output.write(", " + tagname)
        del self._equalsigns[-1]
        self._last2, self._last = self._last, 2

    def text(self, text):
        """Adds text between tag opening and closing.

        Args:
            text: Text to be added.
        """
        if self._last == 1 and not self._equalsigns[-1]:
            self.output.write("\n")
        self.output.write(text)
        self._last2, self._last = self._last, 3

    def _increaseindentation(self):
        """Increases indentation level and adjusts indentation string."""
        self._indentlist.append(self._curindent)
        if not self._equalsigns[-1]:
            self._curindent += self._indent

    def _decreaseindentation(self):
        """Decreases indentation level and adjusts indentation string."""
        self._curindent = self._indentlist.pop()


class HSDStreamFormatter:
    """Reads a HSD feed and writes it on the fly formatted into a stream."""

    def __init__(self, parser, formatter):
        """Intializes HSDFeedPrinter instance.

        Args:
            parser: Event controled parser to be used.
            formatter: Formatter to be used.
        """
        self._parser = parser
        self._formatter = formatter
        self._parser.start_handler = self._formatter.start_tag
        self._parser.close_handler = self._formatter.close_tag
        self._parser.text_handler = self._formatter.text

    def feed(self, fileobj):
        """Feeds the printer with content.

        The contant in fileobj is passed to the parser, and output is generated
        depending on the events.

        Args:
            fileobj: File with HSD-content.
        """
        self._parser.feed(fileobj)


if __name__ == "__main__":
    import io
    from sktools.hsd.parser import HSDParser

    fp = io.StringIO("""
Geometry = GenFormat {
2  S
Ga As
1    1    0.00000000000E+00   0.00000000000E+00   0.00000000000E+00
2    2    0.13567730000E+01   0.13567730000E+01   0.13567730000E+01
0.00000000000E+00   0.00000000000E+00   0.00000000000E+00
0.27135460000E+01   0.27135460000E+01   0.00000000000E+00
0.00000000000E+00   0.27135460000E+01   0.27135460000E+01
0.27135460000E+01   0.00000000000E+00   0.27135460000E+01
}

Hamiltonian = DFTB {
  SCC [unit=None,dim=0] = Yes
  SCCTolerance [toosmall=sure] = 1.0E-007
  MaxSCCIterations = 1000
  Mixer = Broyden {
    MixingParameter = 0.200000000000000
    CachedIterations = -1
  }
  MaxAngularMomentum {
    Ga = "d"
    As = "p"
  }
  Filling = Fermi {
    Temperature = 1.0E-006
    IndependentKFilling = No
  }
  SlaterKosterFiles {
    Ga-Ga = "./Ga-Ga.skf"
    Ga-As = "./Ga-As.skf"
    As-Ga = "./As-Ga.skf"
    As-As = "./As-As.skf"
  }
  KPointsAndWeights {
0.000000000000000E+000 0.000000000000000E+000 0.000000000000000E+000 1.00000000000000
  }
  Charge = 0.000000000000000E+000
  ReadInitialCharges = No
  DampXH = No
  EwaldParameter = 0.000000000000000E+000
  Eigensolver = DivideAndConquer {}
  ThirdOrder = No
}

Options {
  RandomSeed = 0
  WriteHS = No
  ShowFoldedCoords = No
}
""")
    streamformatter = HSDStreamFormatter(HSDParser(),
                                         HSDFormatter(closecomments=True))
    streamformatter.feed(fp)
