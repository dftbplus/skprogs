"""Contains various converters for the query module.
"""


class _HSDConvInt0:

    @staticmethod
    def fromhsd(txt):
        return int(txt)

    @staticmethod
    def tohsd(value):
        return str(value)

#: Converts an 
int0 = _HSDConvInt0


class _HSDConvFloat0:
    formstr = "{:.12E}"

    @staticmethod
    def fromhsd(txt):
        return float(txt)

    @staticmethod
    def tohsd(value):
        return _HSDConvFloat0.formstr.format(value)

float0 = _HSDConvFloat0


class _HSDConvInt1:
    @staticmethod
    def fromhsd(txt):
        words = txt.split()
        return [ _HSDConvInt0.fromhsd(word) for word in words ]

    @staticmethod
    def tohsd(values):
        return " ".join([ _HSDConvInt0.tohsd(val) for val in values ])

int1 = _HSDConvInt1


class _HSDConvFloat1:

    @staticmethod
    def fromhsd(txt):
        words = txt.split()
        return [ _HSDConvFloat0.fromhsd(word) for word in words ]

    @staticmethod
    def tohsd(values):
        return " ".join([ _HSDConvFloat0.tohsd(val) for val in values ])

float1 = _HSDConvFloat1


class _HSDConvStr0:

    @staticmethod
    def fromhsd(txt):
        return txt

    @staticmethod
    def tohsd(value):
        return value

str0 = _HSDConvStr0


class _HSDConvStr1:

    @staticmethod
    def fromhsd(txt):
        return txt.split()

    @staticmethod
    def tohsd(values):
        return " ".join(values)    

str1 = _HSDConvStr1


class _HSDConvBool0:
    truewords = frozenset(("true", "yes", "on"))
    falsewords = frozenset(("false", "no", "off"))
    default_true = "Yes"
    default_false = "No"

    @staticmethod
    def fromhsd(txt):
        lowtxt = txt.lower()
        if lowtxt in _HSDConvBool0.truewords:
            return True
        elif lowtxt in _HSDConvBool0.falsewords:
            return False
        else:
            raise ValueError("Unknown boolean value '{}'".format(txt))

    @staticmethod
    def tohsd(value):
        if value:
            return _HSDConvBool0.default_true
        else:
            return _HSDConvBool0.default_false

bool0 = _HSDConvBool0


class _HSDConvBool1:

    @staticmethod
    def fromhsd(txt):
        words = txt.split()
        return [ _HSDConvBool0.fromhsd(word) for word in words ]

    @staticmethod
    def tohsd(values):
        return " ".join([ _HSDConvBool0.tohsd(value) for value in values ])

bool1 = _HSDConvBool1
