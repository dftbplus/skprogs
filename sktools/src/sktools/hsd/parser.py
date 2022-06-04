from collections import OrderedDict
import sktools.hsd as hsd


__all__ = [ "HSDParser",
            "SYNTAX_ERROR", "UNCLOSED_TAG_ERROR", "QUOTATION_ERROR",
            "BRACKET_ERROR" ]

SYNTAX_ERROR = 1
UNCLOSED_TAG_ERROR = 2
UNCLOSED_OPTION_ERROR = 3
UNCLOSED_QUOTATION_ERROR = 4
ORPHAN_TEXT_ERROR = 5

GENERAL_SPECIALS = "{}[]<=\"'#;"
OPTION_SPECIALS = ",]=\"'#{};"



class HSDParser:
    """Event based parser for the HSD format.

    The methods `start_handler()`, `close_handler()`, `text_handler()`
    and `error_handler()` should be overridden by the actual application.
    """

    def __init__(self, defattrib="default"):
        """Initializes the parser.

        Args:
            defattrib: Name of the default attribute (default: 'default')
        """
        self._fname = ""                   # Name of file being processed
        self._defattrib = defattrib.lower()  # def. attribute name
        self._checkstr = GENERAL_SPECIALS  # special characters to look for
        self._oldcheckstr = ""             # buffer fo checkstr
        self._currenttags = []             # info about opened tags
        self._buffer = []                  # buffering plain text between lines
        self._options = OrderedDict()      # options for current tag
        self._hsdoptions = OrderedDict()   # hsd-options for current tag
        self._key = ""                     # current option name
        self._currline = 0                 # nr. of current line in file
        self._flag_equalsign = False       # last tag was opened with equal sign
        self._flag_option = False          # parser inside option specification
        self._flag_quote = False           # parser inside quotation
        self._flag_haschild = False
        self._oldbefore = ""


    def feed(self, fileobj):
        """Feeds the parser with data.

        Args:
            fileobj: File like object or name of a file containing the data.
        """
        isfilename = isinstance(fileobj, str)
        if isfilename:
            fp = open(fileobj, "r")
            self._fname = fileobj
        else:
            fp = fileobj
        for line in fp.readlines():
            self._parse(line)
            self._currline += 1
        if isfilename:
            fp.close()

        # Check for errors
        if self._currenttags:
            line0 = self._currenttags[-1][1]
        else:
            line0 = 0
        if self._flag_quote:
            self._error(UNCLOSED_QUOTATION_ERROR, (line0, self._currline))
        elif self._flag_option:
            self._error(UNCLOSED_OPTION_ERROR, (line0, self._currline))
        elif self._currenttags:
            self._error(UNCLOSED_TAG_ERROR, (line0, line0))
        elif ("".join(self._buffer)).strip():
            self._error(ORPHAN_TEXT_ERROR, (line0, self._currline))


    def start_handler(self, tagname, options, hsdoptions):
        """Handler which is called when a tag is opened.

        It should be overriden in the application to handle the event in a
        customized way.

        Args:
            tagname: Name of the tag which had been opened.
            options: Dictionary of the options (attributes) of the tag.
            hsdoptions: Dictionary of the options created during the processing
                in the hsd-parser.
        """
        pass


    def close_handler(self, tagname):
        """Handler which is called when a tag is closed.

        It should be overriden in the application to handle the event in a
        customized way.

        Args:
            tagname: Name of the tag which had been closed.
        """
        pass


    def text_handler(self, text):
        """Handler which is called with the text found inside a tag.

        It should be overriden in the application to handle the event in a
        customized way.

        Args:
           text: Text in the current tag.
        """
        pass


    def error_handler(self, error_code, file, lines):
        """Handler which is called if an error was detected during parsing.

        The default implementation throws a HSDException or a descendant of it.

        Args:
            error_code: Code for signalizing the type of the error.
            file: Current file name (empty string if not known).
            lines: Lines between the error occurred.
        """
        error_msg = (
            "Parsing error ({}) between lines {} - {} in file '{}'.".format(
            error_code, lines[0] + 1, lines[1] + 1, file))
        raise hsd.HSDParserError(error_msg)


    def interrupt_handler_hsd(self, command):
        """Handles hsd type interrupt.

        The base class implements following handling: Command is interpreted as
        a file name (quotes eventually removed). A parser is opened with the
        same handlers as the current one, and the given file is feeded in it.

        Args:
            command: Unstripped string as specified in the HSD input after
                the interrupt sign.
        """
        fname = hsd.unquote(command.strip())
        parser = HSDParser(defattrib=self._defattrib)
        parser.start_handler = self.start_handler
        parser.close_handler = self.close_handler
        parser.text_handler = self.text_handler
        parser.feed(fname)


    def interrupt_handler_txt(self, command):
        """Handles text type interrupt.

        The base class implements following handling: Command is interpreted as
        a file name (quotes eventually removed). The file is opened and its
        content is read (without parsing) and added as text.

        Args:
            command: Unstripped string as specified in the HSD input after
                the interrupt sign.

        Returns:
            Unparsed text to be added to the HSD input.
        """
        fname = hsd.unquote(command.strip())
        fp = open(fname, "r")
        txt = fp.read()
        fp.close()
        return txt


    def _parse(self, line):
        """Parses a given line."""

        while True:
            sign, before, after = _splitbycharset(line, self._checkstr)

            # End of line
            if not sign:
                if self._flag_quote:
                    self._buffer.append(before)
                elif self._flag_equalsign:
                    self._text("".join(self._buffer) + before.strip())
                    self._closetag()
                    self._flag_equalsign = False
                elif not self._flag_option:
                    self._buffer.append(before)
                elif before.strip():
                    self._error(SYNTAX_ERROR, (self._currline, self._currline))
                break

            # Special character is escaped
            elif before.endswith("\\") and not before.endswith("\\\\"):
                self._buffer.append(before + sign)

            # Equal sign outside option specification
            elif sign == "=" and not self._flag_option:
                # Ignore if followed by "{" (DFTB+ compatibility)
                if after.lstrip().startswith("{"):
                    self._oldbefore = before
                else:
                    self._flag_haschild = True
                    self._hsdoptions[hsd.HSDATTR_EQUAL] = True
                    self._starttag(before, False)
                    self._flag_equalsign = True

            # Equal sign inside option specification
            elif sign == "=":
                self._key = before.strip()
                self._buffer = []

            # Opening tag by curly brace
            elif sign == "{" and not self._flag_option:
                self._flag_haschild = True
                self._starttag(before, self._flag_equalsign)
                self._buffer = []
                self._flag_equalsign = False

            # Closing tag by curly brace
            elif sign == "}" and not self._flag_option:
                self._text("".join(self._buffer) + before)
                self._buffer = []
                # If 'test { a = 12 }' occurs, curly brace closes two tags
                if self._flag_equalsign:
                    self._flag_equalsign = False
                    self._closetag()
                self._closetag()

            # Closing tag by semicolon
            elif sign == ";" and self._flag_equalsign and not self._flag_option:
                self._flag_equalsign = False
                self._text(before)
                self._closetag()

            # Comment line
            elif sign == "#":
                self._buffer.append(before)
                after = ""

            # Opening option specification
            elif sign == "[" and not self._flag_option:
                if "".join(self._buffer).strip():
                    self._error(SYNTAX_ERROR, (self._currline, self._currline))
                self._oldbefore = before
                self._buffer = []
                self._flag_option = True
                self._key = ""
                self._currenttags.append(("[", self._currline, None))
                self._checkstr = OPTION_SPECIALS

            # Closing option specification
            elif sign == "]" and self._flag_option:
                value = "".join(self._buffer) + before
                key = self._key.lower() if self._key else self._defattrib
                self._options[key] = value.strip()
                self._flag_option = False
                self._buffer = []
                self._currenttags.pop()
                self._checkstr = GENERAL_SPECIALS

            # Quoting strings
            elif sign == "'" or sign == '"':
                if self._flag_quote:
                    self._checkstr = self._oldcheckstr
                    self._flag_quote = False
                    self._buffer.append(before + sign)
                    self._currenttags.pop()
                else:
                    self._oldcheckstr = self._checkstr
                    self._checkstr = sign
                    self._flag_quote = True
                    self._buffer.append(sign)
                    self._currenttags.append(('"', self._currline, None))

            # Closing attribute specification
            elif sign == "," and self._flag_option:
                value = "".join(self._buffer) + before
                key = self._key.lower() if self._key else self._defattrib
                self._options[key] = value.strip()

            # Interrupt
            elif (sign == "<" and not self._flag_option
                  and not self._flag_equalsign):
                txtint = after.startswith("<<")
                hsdint = after.startswith("<!")
                if txtint:
                    self._text("".join(self._buffer) + before)
                    self._buffer = []
                    self.text_handler(self.interrupt_handler_txt(after[2:]))
                    break
                elif hsdint:
                    self.interrupt_handler_hsd(after[2:])
                    break
                else:
                    self._buffer.append(before + sign)

            else:
                self._error(SYNTAX_ERROR, (self._currline, self._currline))

            line = after


    def _text(self, text):
        stripped = text.strip()
        if stripped:
            self.text_handler(stripped)


    def _starttag(self, tagname, closeprev):
        txt = "".join(self._buffer)
        if txt:
            self._text(txt)
        tagname_stripped = tagname.strip()
        if self._oldbefore:
            if tagname_stripped:
                self._error(SYNTAX_ERROR, ( self._currline, self._currline ))
            else:
                tagname_stripped = self._oldbefore.strip()
        if len(tagname_stripped.split()) > 1:
            self._error(SYNTAX_ERROR, (self._currline, self._currline))
        self._hsdoptions[hsd.HSDATTR_LINE] = self._currline
        self._hsdoptions[hsd.HSDATTR_TAG] = tagname_stripped
        tagname_stripped = tagname_stripped.lower()
        self.start_handler(tagname_stripped, self._options, self._hsdoptions)
        self._currenttags.append(
            ( tagname_stripped, self._currline, closeprev, self._flag_haschild))
        self._buffer = []
        self._oldbefore = ""
        self._flag_haschild = False
        self._options = OrderedDict()
        self._hsdoptions = OrderedDict()


    def _closetag(self):
        if not self._currenttags:
            self._error(SYNTAX_ERROR, (0, self._currline))
        self._buffer = []
        tag, line, closeprev, self._flag_haschild = self._currenttags.pop()
        self.close_handler(tag)
        if closeprev:
            self._closetag()

    def _error(self, code, lines):
        self.error_handler(code, self._fname, lines)



def _splitbycharset(txt, charset):
    """Splits a string at the first occurrence of a character in a set.

    Args:
        txt: Text to split.
        chars: Chars to look for.

    Returns:
        Tuple (char, before, after). Char is the character which had been found
        (or empty string if nothing was found). Before is the substring before
        the splitting character (or the entire string). After is the substring
        after the splitting character (or empty string).
    """
    for firstpos, char in enumerate(txt):
        if char in charset:
            break
    else:
        return '', txt, ''
    return txt[firstpos], txt[:firstpos], txt[firstpos + 1:]



def _test_module():
    from io import StringIO
    from sktools.hsd.formatter import HSDStreamFormatter, HSDFormatter
    formatter = HSDFormatter(closecomments=True)
    parser = HSDParser(defattrib="unit")
    streamformatter = HSDStreamFormatter(parser, formatter)
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
Hamiltonian = DFTB {
  SCC = Yes
  SCCTolerance = 1.0E-007
  MaxSCCIterations = 1000
  $MyVariable = 12
  Mixer = Broyden {}
  MaxAngularMomentum = {
    Ga = "d"
    As = "p"
  }
  Filling = Fermi {
    Temperature [Kelvin] = 1.0E-006
  }
  SlaterKosterFiles [format=old] {
    Ga-Ga = "./Ga-Ga.skf"
    Ga-As = "./Ga-As.skf"
    As-Ga = "./As-Ga.skf"
    As-As = "./As-As.skf"
  }
  KPointsAndWeights {
    0.0 0.0 0.0   1.0
  }
}
Options {
  AtomResolvedEnergies = No
  RestartFrequency = 20
  RandomSeed = 0
  WriteHS = No
}
""")
    streamformatter.feed(stream)


if __name__ == "__main__":
    _test_module()
