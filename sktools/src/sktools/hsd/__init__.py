"""Implements various functionalities for creating and querying the
HSD (Human readable Structured Data) format.
"""


class HSDException(Exception):
    """Base class for exceptions in the HSD package."""
    pass


class HSDQueryError(HSDException):
    """Base class for errors detected by the HSDQuery object.


    Attributes:
        filename: Name of the file where error occured (or empty string).
        line: Line where the error occurred (or -1).
        tag: Name of the tag with the error (or empty string).
    """

    def __init__(self, msg="", node=None):
        """Initializes the exception.

        Args:
            msg: Error message
            node: HSD element where error occured (optional).
        """
        super().__init__(msg)
        if node is not None:
            self.tag = node.gethsd(HSDATTR_TAG, node.tag)
            self.file = node.gethsd(HSDATTR_FILE, -1)
            self.line = node.gethsd(HSDATTR_LINE, None)
        else:
            self.tag = ""
            self.file = -1
            self.line = None


class HSDMissingTagException(HSDQueryError): pass
class HSDInvalidTagException(HSDQueryError): pass
class HSDInvalidTagValueException(HSDQueryError): pass
class HSDMissingAttributeException(HSDQueryError): pass
class HSDInvalidAttributeException(HSDQueryError): pass
class HSDInvalidAttributeValueException(HSDQueryError): pass


class HSDParserError(HSDException):
    """Base class for parser related errors."""
    pass


def unquote(txt):
    """Giving string without quotes if enclosed in those."""
    if len(txt) >= 2 and (txt[0] in "\"'") and txt[-1] == txt[0]:
        return txt[1:-1]
    else:
        return txt


HSDATTR_PROC = "processed"
HSDATTR_EQUAL = "equal"
HSDATTR_FILE = "file"
HSDATTR_LINE = "line"
HSDATTR_TAG = "tag"