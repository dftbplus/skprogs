import xml.etree.ElementTree as ETree
import sktools.hsd as hsd
from sktools.hsd.formatter import HSDFormatter


HSD_ATTRIB_PREFIX = "_hsd_"
HSD_ATTRIB_PREFIX_LEN = len(HSD_ATTRIB_PREFIX)



class HSDTree(ETree.ElementTree):
    """Wrapper around an entire tree."""  

    def writehsd(self, file, formatter=None):
        """Writes the tree in HSD format.
        
        Args:
            file: File to write the tree into. It can be file name or a
                file like object.
        """
        if formatter is None:
            formatter = HSDFormatter()
        isfilename = isinstance(file, str)
        fp = open(file, "w") if  isfilename else file
        formatter.set_output(fp)
        self.processtree(formatter)
        if isfilename:
            fp.close()


    def processtree(self, processor):
        self._processtree(self.getroot(), processor)


    def _processtree(self, parent, processor):
        """Private helper routine for writehsd."""
        if parent.text:
            processor.text(parent.text)
        else:
            for child in parent:
                attribs = dict(child.items())
                hsdattribs = dict(child.hsditems())
                processor.start_tag(child.tag, attribs, hsdattribs)
                self._processtree(child, processor)
                processor.close_tag(child.tag)



class Element(ETree.Element):
    """Element Interface containing extra dictionary with hsd attributes."""

    def __init__(self, tag, attribs=None, hsdattribs=None):
        if attribs is None:
            attribs = {}
        super().__init__(tag, attribs)
        self.sethsdfromdict(hsdattribs)

    def sethsdfromdict(self, hsdattribs):
        if hsdattribs is not None:
            for key, value in hsdattribs.items():
                self.sethsd(key, value)

    def gethsd(self, attribname, defaultvalue=None):
        return super().get(HSD_ATTRIB_PREFIX + attribname, defaultvalue)

    def sethsd(self, attribname, attribvalue):
        super().set(HSD_ATTRIB_PREFIX + attribname, attribvalue)

    def hsdkeys(self):
        return [ key[HSD_ATTRIB_PREFIX_LEN:]
                 for key in super().keys()
                 if key.startswith(HSD_ATTRIB_PREFIX) ]

    def hsditems(self):
        return [ ( key, self.gethsd(key) ) for key in self.hsdkeys() ]

    def get(self, key, default=None):
        if key.startswith(HSD_ATTRIB_PREFIX):
            msg = "Invalid attribute name '{}'".format(key)
            raise hsd.HSDException(msg)
        return super().get(key, default)

    def set(self, key, value):
        if key.startswith(HSD_ATTRIB_PREFIX):
            msg = "Invalid attribute name '{}'".format(key)
            raise hsd.HSDException(msg)
        super().set(key, value)

    def keys(self):
        return [ key for key in super().keys()
                 if not key.startswith(HSD_ATTRIB_PREFIX) ]

    def items(self):
        return [ ( key, self.get(key) ) for key in self.keys() ]

    def makeelement(tag, attrib, hsdattrib=None):
        return Element(tag, attrib, hsdattrib)
    


def SubElement(parent, tag, attribs=None, hsdattribs=None):
    """Subelement factory with extra hsd attributes."""
    element = parent.makeelement(tag, attribs, hsdattribs)
    parent.append(element)
    return element
