import sktools.hsd.converter as conv
import sktools.common as sc


class EquidistantGrid(sc.ClassDict):
    """Equidistant grid.

    Attributes
    ----------
    gridstart : float
        Starting point of the grid.
    gridseparation : float
        Distance between grid points.
    maxdistance : float
        Maximal grid distance.
    tolerance:
        Stopping criterion for grid (when represented value on the grid is
        below tolerance).
    """

    @classmethod
    def fromhsd(cls, node, query):
        myself = cls()
        myself.gridstart = query.getvalue(node, "gridstart", conv.float0)
        myself.gridseparation = query.getvalue(node, "gridseparation",
                                               conv.float0)
        myself.tolerance = query.getvalue(node, "tolerance", conv.float0)
        myself.maxdistance = query.getvalue(node, "maxdistance", conv.float0)
        return myself

    def __eq__(self, other):
        if abs(self.gridstart - other.gridstart) > sc.INPUT_FLOAT_TOLERANCE:
            return False
        if (abs(self.gridseparation - other.gridseparation)
                > sc.INPUT_FLOAT_TOLERANCE):
            return False
        if abs(self.tolerance - other.tolerance) > sc.INPUT_FLOAT_TOLERANCE:
            return False
        if abs(self.maxdistance - other.maxdistance) > sc.INPUT_FLOAT_TOLERANCE:
            return False
        return True


TWOCENTER_GRIDS = {
    "equidistantgrid": EquidistantGrid
}
