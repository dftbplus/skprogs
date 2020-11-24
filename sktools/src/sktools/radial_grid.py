import numpy as np
import sktools.common as sc


class RadialGrid:

    FLOAT_TOLERANCE = 1e-8

    def __init__(self, rr, weights):
        if len(rr) != len(weights):
            raise ValueError("Length of radial grid and weights not compatible")
        self.nr = len(rr)
        self.rr = rr
        self.weights = weights

    def __eq__(self, other):
        if self.nr != other.nr:
            return False
        if np.max(np.abs(self.rr - other.rr)) > self.FLOAT_TOLERANCE:
            return False
        if np.max(np.abs(self.weights - other.weights)) > self.FLOAT_TOLERANCE:
            return False
        return True

    def dot(self, f1, f2):
        return np.sum(self.rr * self.rr * self.weights * f1 * f2)


class GridData:

    FLOAT_FORMAT = "{:21.12E}"

    def __init__(self, grid, data):
        if grid.nr != len(data):
            raise ValueError("Incompatible grids")
        self.grid = grid
        self.data = np.reshape(data, (grid.nr, -1))


    def tofile(self, fobj):
        with sc.FileFromStringOrHandler(fobj, "w") as fp:
            fp.write("{:d}\n".format(self.grid.nr))
            ndata = len(self.data[0])
            formstr = self.FLOAT_FORMAT * (ndata + 2) + "\n"
            for ii in range(self.grid.nr):
                fp.write(formstr.format(self.grid.rr[ii], self.grid.weights[ii],
                                        *self.data[ii]))


VNUC = 0
VHARTREE = 1
VXCUP = 2
VXCDOWN = 3