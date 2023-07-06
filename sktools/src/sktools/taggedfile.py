from collections import OrderedDict
import numpy as np


class TaggedFile(OrderedDict):

    CONVERTER = { "real": np.float64,
                  "integer": np.int64,
                  "logical": lambda x: x.lower() == "t"
                  }

    DTYPE_NAMES = {
        np.dtype("int64"): "integer",
        np.dtype("float64"): "real",
        np.dtype("bool"): "logical",
        int: "integer",
        float: "real",
        bool: "logical",
        str: "string",
    }

    DTYPE_FORMATS = {
        np.dtype("int64"): ( 3, "{:20d}"),
        np.dtype("float64"): (3, "{:23.15E}"),
        np.dtype("bool"): (40, "{:2s}"),
        str: (72, "  {:s}"),
        int: "  {:20d}",
        float: "  {:23.15E}",
        bool: "  {:s}",
    }

    def __init__(self, initvalues=None):
        if initvalues:
            super().__init__(initvalues)
        else:
            super().__init__()

    def tofile(self, fp):
        for key, value in self.items():
            if isinstance(value, np.ndarray):
                dtype = value.dtype
                fp.write(
                    "@{:s}:{:s}:{:d}:{:s}\n".format(
                        key,
                        self.DTYPE_NAMES[dtype],
                        len(value.shape),
                        ",".join([str(dd) for dd in value.shape]),
                    )
                )
                nitem, formstr = self.DTYPE_FORMATS[dtype]
                if dtype == np.dtype("bool"):
                    value = np.where(value, "T", "F")
                self._lineformattedwrite(fp, nitem, formstr, value.ravel(order="F"))
            elif isinstance(value, str):
                dtype = str
                nn = len(value)
                fp.write(
                    "@{:s}:{:s}:{:d}:{:d}\n".format(key, self.DTYPE_NAMES[dtype], 1, nn)
                )
                nitem, formstr = self.DTYPE_FORMATS[dtype]
                for ii in range(nitem, nn + 1, nitem):
                    fp.write(formstr.format(value[ii - nitem : ii]) + "\n")
                remaining = nn % nitem
                if remaining:
                    fp.write(formstr.format(value[nn - remaining : nn]) + "\n")
            else:
                dtype = type(value)
                fp.write("@{:s}:{:s}:{:d}:\n".format(key, self.DTYPE_NAMES[dtype], 0))
                if isinstance(value, bool):
                    value = "T" if value else "F"
                fp.write(self.DTYPE_FORMATS[dtype].format(value))
                fp.write("\n")

    def _lineformattedwrite(self, fp, nitem, formstr, valuelist):
        nn = len(valuelist)
        lineformstr = (
            " ".join(
                nitem
                * [
                    formstr,
                ]
            )
            + "\n"
        )
        for ii in range(nitem, nn + 1, nitem):
            fp.write(lineformstr.format(*valuelist[ii - nitem : ii]))
        remaining = nn % nitem
        if remaining:
            lineformstr = (
                " ".join(
                    remaining
                    * [
                        formstr,
                    ]
                )
                + "\n"
            )
            fp.write(lineformstr.format(*valuelist[nn - remaining : nn]))

    @classmethod
    def fromfile(cls, fp, transpose=False):
        fname = isinstance(fp, str)
        if fname:
            fp = open(fp, "r")
        tagvalues = []
        line = fp.readline()
        tmp = []
        while line:
            tagline = line
            tmp = []
            line = fp.readline()
            while line and line[0] != "@":
                tmp += line.split()
                line = fp.readline()
            words = tagline.split(":")
            tag = words[0][1:]
            dtype = words[1]
            dim = int(words[2])
            if dim:
                shape = [int(dd) for dd in words[3].split(",")]
                if dtype == "string":
                    value = "".join(tmp)
                else:
                    elems = [cls.CONVERTER[dtype](ss) for ss in tmp]
                    value = np.array(elems).reshape(shape, order="F")
                    if transpose:
                        value = value.transpose()
            else:
                value = cls.CONVERTER[dtype](tmp[0])
            tagvalues.append((tag, value))
        if fname:
            fp.close()
        return cls(tagvalues)
