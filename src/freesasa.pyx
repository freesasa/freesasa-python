# -*- mode: python; python-indent-offset: 4 -*-
#
# The cython directives will fail if we don't have a few lines of comments above them. (Why?)
#
# cython: c_string_type=str, c_string_encoding=ascii

"""
The :py:mod:`freesasa` python module wraps the FreeSASA `C API`_
"""

from libc.stdio cimport FILE, fopen, fclose
from libc.stdlib cimport free, realloc, malloc
from libc.string cimport memcpy
from cfreesasa cimport *

include "parameters.pyx"
include "result.pyx"
include "classifier.pyx"
include "structure.pyx"

## Used for classification
polar = 'Polar'

## Used for classification
apolar = 'Apolar'

## int: Suppress all warnings and errors (used by setVerbosity())
silent = FREESASA_V_SILENT

## int: Suppress all warnings but not errors (used by setVerbosity())
nowarnings = FREESASA_V_NOWARNINGS

## int: Normal verbosity (used by setVerbosity())
normal = FREESASA_V_NORMAL

## int: Print debug messages (used by setVerbosity())
debug = FREESASA_V_DEBUG


def calc(structure,parameters=None):
    """
    Calculate SASA of Structure

    Args:
        structure: :py:class:`.Structure` to be used
        parameters: :py:class:`.Parameters` to use (if not specified defaults are used)

    Returns:
        :py:class:`.Result`: The results

    Raises:
        Exception: something went wrong in calculation (see C library error messages)
    """
    cdef const freesasa_parameters *p = NULL
    cdef const freesasa_structure *s = NULL
    if parameters is not None:  parameters._get_address(<size_t>&p)
    structure._get_address(<size_t>&s)
    result = Result()
    result._c_result = <freesasa_result*> freesasa_calc_structure(s,p)
    result._c_root_node = <freesasa_node*> freesasa_tree_init(result._c_result,
                                                              s, "Structure")
    if result._c_result is NULL:
        raise Exception("Error calculating SASA.")

    return result

def calcCoord(coord, radii, parameters=None):
    """
    Calculate SASA for a set of coordinates and radii

    Args:
        coord (list): array of size 3*N with atomic coordinates
           `(x1, y1, z1,  x2, y2, z2, ..., x_N, y_N, z_N)`.
        radii (list): array of size N with atomic radii `(r_1, r_2, ..., r_N)`.
        parameters: :py:class:`.Parameters` to use (if not specified, defaults are used)
    Raises:
        AssertionError: mismatched array-sizes
        Exception: Out of memory
        Exception: something went wrong in calculation (see C library error messages)
    """
    assert(len(coord) == 3*len(radii))

    cdef const freesasa_parameters *p = NULL
    cdef double *c = <double*> malloc(len(coord)*sizeof(double))
    cdef double *r = <double*> malloc(len(radii)*sizeof(double))
    if c is NULL or r is NULL:
        raise Exception("Memory allocation error")

    for i in xrange(len(coord)):
        c[i] = coord[i]
    for i in xrange(len(radii)):
        r[i] = radii[i]

    if parameters is not None: parameters._get_address(<size_t>&p)

    result = Result()
    result._c_result = <freesasa_result*> freesasa_calc_coord(c, r, len(radii), p)

    if result._c_result is NULL:
        raise Exception("Error calculating SASA.")

    free(c)
    free(r)

    return result

def classifyResults(result,structure,classifier=None):
    """
    Break SASA result down into classes.

    Args:
        result: :py:class:`.Result` from SASA calculation.
        structure: :py:class:`Structure` used in calculation.
        classifier: :py:class:`.Classifier` to use (if not specified default is used).

    Returns:
        dict: Dictionary with names of classes as keys and their SASA values as values.

    Raises:
        Exception: Problems with classification, see C library error messages
            (or Python exceptions if run with derived classifier).
    """
    if classifier is None:
        classifier = Classifier()
    ret = dict()
    for i in range(0,structure.nAtoms()):
        name = classifier.classify(structure.residueName(i),structure.atomName(i))
        if name not in ret:
            ret[name] = 0
        ret[name] += result.atomArea(i)
    return ret

def selectArea(commands, structure, result):
    """
    Sum SASA result over a selection of atoms

    Args:
        commands (list): A list of commands with selections using Pymol
            syntax, e.g. ``"s1, resn ala+arg"`` or ``"s2, chain A and resi 1-5"``.
            See `select-syntax`_.
        structure: A :py:class:`.Structure`.
        result: :py:class:`.Result` from sasa calculation on structure.

    Returns:
        dict: Dictionary with names of selections (``"s1"``, ``"s2"``, ...) as
        keys, and the corresponding SASA values as values.

    Raises:
        Exception: Parser failed (typically syntax error), see
           library error messages.
    """
    cdef freesasa_structure *s
    cdef freesasa_result *r
    cdef freesasa_selection *selection
    structure._get_address(<size_t> &s)
    result._get_address(<size_t> &r)
    value = dict()
    for cmd in commands:
        selection = freesasa_selection_new(cmd, s, r)
        if selection == NULL:
            raise Exception("Error parsing '%s'" % cmd)
        value[freesasa_selection_name(selection)] = freesasa_selection_area(selection)
        freesasa_selection_free(selection)
    return value

def setVerbosity(verbosity):
    """
    Set global verbosity

    Args:
        verbosity (int): Can have values :py:const:`.silent`, :py:const:`.nowarnings`
            or :py:const:`.normal`
    Raises:
        AssertionError: if verbosity has illegal value
    """
    assert(verbosity in [silent, nowarnings, normal])
    freesasa_set_verbosity(verbosity)


def getVerbosity():
    """
    Get global verbosity

    Returns:
        int: Verbosity :py:const:`.silent`, :py:const:`.nowarnings`
        or :py:const:`.normal`
    """
    return freesasa_get_verbosity()

def calcBioPDB(bioPDBStructure, parameters = Parameters(),
               classifier = None, options = Structure.defaultOptions):
    """
    Calc SASA from `BioPython` PDB structure.

    Usage::

        result, sasa_classes, residue_areas = calcBioPDB(structure, ...)

    Experimental, not thorougly tested yet

    Args:
        bioPDBStructure: A `Bio.PDB` structure
        parameters: A :py:class:`.Parameters` object (uses default if none specified)
        classifier: A :py:class:`.Classifier` object (uses default if none specified)
        options (dict): Options supported are 'hetatm', 'skip-unknown' and 'halt-at-unknown'
            (uses :py:attr:`.Structure.defaultOptions` if none specified

    Returns:
        A :py:class:`.Result` object, a dictionary with classes
        defined by the classifier and associated areas,
        and a dictionary of the type returned by :py:meth:`.Result.residueAreas`.

    Raises:
        Exception: if unknown atom is encountered and the option
            'halt-at-unknown' is active. Passes on exceptions from
            :py:func:`.calc()`, :py:func:`.classifyResults()` and
            :py:func:`.structureFromBioPDB()`.
    """
    structure = structureFromBioPDB(bioPDBStructure, classifier, options)
    result = calc(structure, parameters)

    # Hack!:
    # This calculation depends on the structure not having been deallocated,
    # calling it later will cause seg-faults. By calling it now the result
    # is stored.
    # TODO: See if there is a refactoring that solves this in a more elegant way
    # residue_areas = result.residueAreas()

    sasa_classes = classifyResults(result, structure, classifier)
    return result, sasa_classes #, residue_areas
