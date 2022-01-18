'''Common functionality used by the project.'''


import sys
import re
import os.path
import shelve
import dbm
import tempfile
import shutil
import logging

import sktools.hsd as hsd
import sktools.hsd.converter as conv


LOGGER = logging.getLogger('common')


# Maximal angular momentum
MAX_ANGMOM = 4

# Translate between angular momentum and shell name
ANGMOM_TO_SHELL = ['s', 'p', 'd', 'f', 'g']

# Translate between shell name and angular momentum
SHELL_TO_ANGMOM = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}

# Name of the spin channels
SPIN_NAMES = ['u', 'd']

# Max. principal quantum number
MAX_PRINCIPAL_QN = 7

RELATIVISTICS_NONE = 0
RELATIVISTICS_ZORA = 1
RELATIVISTICS_TYPES = {'none': RELATIVISTICS_NONE,
                       'zora': RELATIVISTICS_ZORA}

# XC_FUNCTIONAL_LDA = 2
# XC_FUNCTIONAL_PBE = 3
# XC_FUNCTIONAL_BLYP = 4
# XC_FUNCTIONAL_LCPBE = 5
# XC_FUNCTIONAL_LCBNL = 6
# XC_FUNCTIONAL_TYPES = {'lda': XC_FUNCTIONAL_LDA,
#                        'pbe': XC_FUNCTIONAL_PBE,
#                        'blyp': XC_FUNCTIONAL_BLYP,
#                        'lcpbe': XC_FUNCTIONAL_LCPBE,
#                        'lcbnl': XC_FUNCTIONAL_LCBNL}

XC_FUNCTIONAL_LDA = 0
XC_FUNCTIONAL_PBE = 1
XC_FUNCTIONAL_RSH = 2
XC_FUNCTIONAL_TYPES = {"lda": XC_FUNCTIONAL_LDA,
                       "pbe": XC_FUNCTIONAL_PBE,
                       "rsh": XC_FUNCTIONAL_RSH}

SUPERPOSITION_POTENTIAL = 0
SUPERPOSITION_DENSITY = 1
SUPERPOSITION_TYPES = {'potential': SUPERPOSITION_POTENTIAL,
                       'density': SUPERPOSITION_DENSITY}

WAVEFUNC_FILE_NAME_FORMAT = 'wave_{:02d}{:s}.dat'
POTENTIAL_FILE_NAME = 'pot.dat'
DENSITY_FILE_NAME = 'dens.dat'


# Tolerance for float numbers in user input
INPUT_FLOAT_TOLERANCE = 1E-8


class SkgenException(Exception):
    '''Custom exception of the skgen script.'''


def openfile(fobj, mode):
    '''Opens a file or passes a file object.

    Args:

        fobj (file object): file object
        mode (str): mode to open file in

    Returns:

        fp (file object): file object
        isfname (bool): true, if file object got opened from file name

    '''

    isfname = isinstance(fobj, str)

    if isfname:
        fp = open(fobj, mode)
    else:
        fp = fobj

    return fp, isfname


def writefloats(fp, nums, indent=0, indentstr=None, numperline=4,
                formstr='{:23.15E}'):
    '''Writes (nested) data array to formatted file.

    Args:

        fp (file object): file object
        nums (ndarray): data
        indent (int): number of space indentations while writing data
        indentstr (str): if none, indentation string build from indent
        numperline (int): number of values to write per line
        formstr (str): string formatter

    '''

    if indentstr is None:
        indentstr = ' ' * indent

    lineform = indentstr + formstr * numperline + '\n'
    nums1d = nums.flat
    nnumber = len(nums1d)
    nline = nnumber // numperline

    for ii in range(nline):
        fp.write(
            lineform.format(*nums1d[ii * numperline:(ii + 1) * numperline]))

    res = nnumber % numperline
    if res:
        lineform = indentstr + formstr * res + '\n'
        fp.write(lineform.format(*nums1d[nnumber - res:nnumber]))


# Fortran float pattern with possibility for reccurance
PAT_FORTRAN_FLOAT = re.compile(
    r'^(?:(?P<occurance>[0-9]+)\*)?(?P<value>[+-]?\d*\.?\d*(?:[eE][+-]?\d+)?)$')
PAT_FORTRAN_SEPARATOR = re.compile(r'[,]?\s+')


def split_fortran_fields(sep, maxsplit=0):
    ''''Splits a line containing Fortran (numeric) fields.

    Args:

        sep (str): separator to use when splitting the string
        maxsplit (int): maximum number of splits allowed

    '''

    return [field for field in
            PAT_FORTRAN_SEPARATOR.split(sep, maxsplit=maxsplit)
            if len(field) > 1]


def convert_fortran_floats(txt):
    '''Converts floats in fortran notation to intrinsic floats.'''

    result = []
    words = split_fortran_fields(txt)
    for word in words:
        match = PAT_FORTRAN_FLOAT.match(word)
        if not match:
            result.append(None)
            continue
        occ = match.group('occurance')
        if occ is not None:
            occ = int(occ)
        else:
            occ = 1
        val = float(match.group('value'))
        result += [val,] * occ
    return result


# Shell name pattern
PAT_SHELLNAME = re.compile(r'^(?P<n>[0-9])(?P<shell>[spdfg])$')


def shell_name_to_ind(txt):
    '''Converts a named shell (e.g. '1s') into (n, l) tuple (e.g. (1, 0)).

    Parameters
    ----------
    txt : str
        Text to parse.

    Returns
    -------
    n : int
        Principal quantum number
    l : int
        Angular momentum

    Raises
    ------
    ValueError
        If conversion was not successfull.

    '''

    match = PAT_SHELLNAME.match(txt)
    if not match:
        raise ValueError("Invalid shell name '{}'".format(txt))

    return int(match.group('n')), SHELL_TO_ANGMOM[match.group('shell')]


def shell_ind_to_name(nn, ll):
    '''Converts the shell index, i.e. angular momentum, to the shell string.

    Args:

        nn (int): principal quantum number, i.e. 1, 2, 3, ...
        ll (int): angular momentum quantum number, i.e. ll = 0, ..., nn - 1

    Returns:

        shell string

    '''

    return '{:d}{}'.format(nn, ANGMOM_TO_SHELL[ll])


class FileFromStringOrHandler:
    '''Class that handles file I/O based on a handler or filename.'''

    def __init__(self, fname_or_handler, mode):
        '''Initializes a FileFromStringOrHandler object.'''

        if isinstance(fname_or_handler, str):
            self._fp = open(fname_or_handler, mode)
            self._tobeclosed = True
        else:
            self._fp = fname_or_handler
            self._tobeclosed = False

    def __enter__(self):
        '''Overload __enter__ function.'''
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        '''Overload __exit__ function.'''
        if self._tobeclosed:
            self._fp.close()

    def write(self, *args, **kwargs):
        '''Overload write function.'''
        return self._fp.write(*args, **kwargs)

    def writelines(self, *args, **kwargs):
        '''Overload writelines function.'''
        return self._fp.writelines(*args, **kwargs)

    def read(self, *args, **kwargs):
        '''Overload read function.'''
        return self._fp.read(*args, **kwargs)

    def readline(self, *args, **kwargs):
        '''Overload readline function.'''
        return self._fp.readline(*args, **kwargs)

    def readlines(self, *args, **kwargs):
        '''Overload readlines function.'''
        return self._fp.readlines(*args, **kwargs)


class ClassDict:
    '''Dictionary like object accessible in class notation.'''

    def __init__(self, initdata=None):
        '''Initializes a ClassDict object.'''

        self._dict = {}
        if initdata is not None:
            self._dict.update(initdata)

    def __setattr__(self, key, value):
        if key.startswith('_'):
            super().__setattr__(key, value)
        else:
            self[key] = value

    def __getattr__(self, item):
        if item.startswith('_'):
            return super().__getattribute__(item)
        return self[item]

    def __contains__(self, item):
        return item in self._dict

    def __setitem__(self, key, value):
        self._dict[key] = value

    def __getitem__(self, item):
        try:
            return self._dict[item]
        except KeyError:
            pass
        msg = '{} instance has no key/attribute "{}"'.format(
            self.__class__.__name__, item)
        raise KeyError(msg)

    def __iter__(self):
        return iter(self._dict)

    def __len__(self):
        return len(self._dict)

    def __eq__(self, other):
        if isinstance(other, ClassDict):
            return self._dict == other._dict
        return self._dict == other

    def update(self, other):
        '''Adds other iterable to the dictionary.'''
        self._dict.update(other._dict)


    def get(self, key, default=None):
        '''Returns the value of the item with the specified key.'''
        return self._dict.get(key, default)


    def keys(self):
        '''Returns view that contains the keys of the dictionary.'''
        return self._dict.keys()


def fatalerror(msg, errorcode=-1):
    '''Issue error message and exit.

    Args:

        msg (str): error message
        errorcode (int): error code to raise

    '''

    LOGGER.critical(msg)
    sys.exit(errorcode)


def get_shellvalues(node, query):
    '''Returns dictionary with the values assigned to individual shells.

    Args:

        node (Element): parent node
        query (HSDQuery): queries an HSD-tree

    Returns:

        values (dict): dictionary with the values assigned to individual shells

    '''

    values = {}

    for child in node:

        try:
            shell = shell_name_to_ind(child.tag)
        except ValueError:
            raise hsd.HSDInvalidTagException(
                msg="Invalid shell name '{}'".format(child.tag), node=child)

        value = query.getvalue(child, '.', conv.float0)
        values[shell] = value

    return values


def get_shellvalues_list(node, query, converter):
    '''Returns a list of converted shell values.

    Args:

        node (Element): parent node
        query (HSDQuery): queries an HSD-tree
        converter (converter object): object with methods fromhsd() and tohsd()
            which can convert between the hsd element and the desired type

    Returns:

        values (list): list of values with their type depending on the converter

    '''

    values = []

    for shellname in ANGMOM_TO_SHELL:
        shellnode = query.findchild(node, shellname, optional=True)
        if shellnode is None:
            break
        value = query.getvalue(shellnode, '.', converter)
        values.append(value)

    return values


def hsd_node_factory(classtype, classes, node, query):
    '''Creates an object depending on the node and a class dictionary.

    Parameters
    ----------
    classtype : str
        Textual name of the class to create for error messages
        (e.g. 'twocenter calculator')
    classes : dict
        Contains classes (not instances!) with their corresponding hsd-name.
        Each must have a `fromhsd(node, query)` class method.
    node : Element
        HSD representation of the node.
    query : query object
        Query object to use.

    Returns
    -------
    node : Element or None
        Returns the element created using the hsd input in the node or None
        if the node passed was None.
    '''

    if node is None:
        return None
    myclass = classes.get(node.tag)
    if myclass is None:
        raise hsd.HSDInvalidTagException(
            "Unknown {} '{}'".format(classtype, node.tag))

    return myclass.fromhsd(node, query)


def store_as_shelf(fname, shelfdict=None, **kwargs):
    '''Stores the given keyword arguments in a shelf.

    Parameters
    ----------
    fname : str
        Name of the file which will contain the shelf content.
    shelfdict : dict, optional
        Dictionary with values to be stored in the shelf file.
    **kwargs : arbitrary, optional
        Keyword value pairs to be stored in the shelf file.
    '''

    db = shelve.open(fname, 'n')
    if shelfdict is not None:
        for key, value in shelfdict.items():
            db[key] = value
    for key in kwargs:
        db[key] = kwargs[key]
    db.close()


def retrive_from_shelf(fname):
    '''Open dictionary from shelf.'''

    db = shelve.open(fname, 'r')
    resdict = dict(db)
    db.close()

    return resdict


def create_unique_workdir(workroot, subdirprefix):
    '''Create uniquely named directory.

    Args:

        workroot (str): root directory where to create temporary directory
        subdirprefix (str): file name will begin with this prefix

    Returns:

        workdir (str): created temporary directory

    '''

    workdir = tempfile.mkdtemp(prefix=subdirprefix, dir=workroot)
    LOGGER.debug('Created working directory %s', workdir)

    return workdir


def create_workdir(workdir, reuse_existing=False):
    '''Creates a working directory.

    Parameters
    ----------
    workdir : str
        Working directory to create. If directory already exists, it will
        be deleted, unless reuse_existing is set to True.
    resuse_existing : bool, optional
        Reuse if working directory already exists.
    '''

    if os.path.exists(workdir):
        if reuse_existing:
            return
        LOGGER.debug('Removing existing working directory %s', workdir)
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    LOGGER.debug('Created working directory %s', workdir)


def find_dir_with_matching_shelf(search_dirs, shelf_file, **kwargs):
    '''Returns the directory containing a shelve with given content.

    Paramters
    ---------
    search_dirs : directories to scan
        Directories to scan.
    shelf_file : str
        Name of the file containing the shelve.
    **kwargs : arbitrary
        Name and content of the items the shelve should contain.

    Returns
    -------
    directory : str
        The directory, where a shelve file containing at least the given
        content exist. If no such directory was found, None is returned.
    '''

    for directory in search_dirs:
        if is_shelf_file_matching(os.path.join(directory, shelf_file), kwargs):
            return directory

    return None


def is_shelf_file_matching(shelf_file, mydict):
    '''Returns true, if the dictionary in shelf file matches reference.'''

    try:
        db = shelve.open(shelf_file, 'r')
    except dbm.error:
        return False
    match = True
    for key in mydict:
        match = key in db and db[key] == mydict[key]
        if not match:
            return False
    return True


def get_dirs_with_matching_shelf(search_dirs, shelf_file, **kwargs):
    '''Searches multiple directories for given shelf and returns those with
    matching entries.'''

    matching_dirs = []
    for directory in search_dirs:
        shelf_path = os.path.join(directory, shelf_file)
        if is_shelf_file_matching(shelf_path, kwargs):
            matching_dirs.append(directory)

    return matching_dirs


def shelf_exists(shelf_name):
    '''Infers whether given dictionary-like object exists in shelve.

    Args:

        shelf_name (dict): dictionary-like object

    Returns:

        result (bool): true, if dictionary-like object exists in shelve

    '''

    try:
        db = shelve.open(shelf_name, 'r')
    except dbm.error:
        result = False
    else:
        db.close()
        result = True

    return result


def capitalize_elem_name(elem):
    '''Converts element name into a capitalized one.

    Args:

        elem (str): element string to convert

    Returns:

        proper, capitalized element name

    '''

    return elem[0].upper() + elem[1:].lower()


class ScriptLogFormatter(logging.Formatter):
    '''Defines the general log formatting.'''

    log_formats = {
        logging.CRITICAL: '!!! [{logrecord.name}] {logrecord.message}',
        logging.ERROR: '!!! [{logrecord.name}] {logrecord.message}',
        logging.WARNING: '! [{logrecord.name}] {logrecord.message}',
        logging.INFO: '[{logrecord.name}] {logrecord.message}',
        logging.DEBUG: '[{logrecord.name}] {logrecord.message}'
    }
    default_log_format = '{logrecord.levelno}: {logrecord.message}'


    def __init__(self):
        super().__init__('{message}', style='{')


    def format(self, logrecord):
        # Make sure, message attribute of logrecord is generated
        super().format(logrecord)
        formatstr = self.log_formats.get(logrecord.levelno,
                                         self.default_log_format)
        result = formatstr.format(logrecord=logrecord)
        return result


def log_path(path):
    '''Generate path shown in logging messages.

    Args:

        path (str): path to build message string from

    Returns:

        pathname (str): modified pathname of logging message

    '''

    cwd = os.path.curdir
    pathname_abs = os.path.abspath(path)
    pathname_rel = os.path.relpath(path, cwd)
    if len(pathname_abs) < len(pathname_rel):
        pathname = pathname_abs
    else:
        pathname = pathname_rel

    pathname = '(' + pathname + ')'

    return pathname


def get_script_logger(loglevel, scriptname):
    '''Generate script logger with proper loglevel.

    Args:

        loglevel (str): logging level, i.e. debug, info, warning, error
        scriptname (str): name of the current script

    Returns:

        logger (logger): script logger

    '''

    loghandler = logging.StreamHandler()
    myformatter = ScriptLogFormatter()
    loghandler.setFormatter(myformatter)
    logging.root.addHandler(loghandler)
    numeric_level = getattr(logging, loglevel.upper(), None)
    logging.root.setLevel(numeric_level)
    logger = logging.getLogger(name=scriptname)

    return logger


def get_script_name():
    '''Returns the name of the invoked script.'''
    return os.path.basename(sys.argv[0])
