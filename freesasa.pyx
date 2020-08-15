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

## Used to specify the algorithm by Shrake & Rupley
ShrakeRupley = 'ShrakeRupley'

## Used to specify the algorithm by Lee & Richards
LeeRichards = 'LeeRichards'

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


cdef class Parameters:
    """
    Stores parameter values to be used by calculation.

    Default parameters are
    ::
        defaultParameters = {
            'algorithm'    : LeeRichards,
            'probe-radius' : freesasa_default_parameters.probe_radius,
            'n-points'     : freesasa_default_parameters.shrake_rupley_n_points,
            'n-slices'     : freesasa_default_parameters.lee_richards_n_slices,
            'n-threads'    : freesasa_default_parameters.n_threads
        }

    Attributes:
        defaultParamers (dict): The default parameters
    """


    cdef freesasa_parameters _c_param

    defaultParameters = {
        'algorithm'    : LeeRichards,
        'probe-radius' : freesasa_default_parameters.probe_radius,
        'n-points'     : freesasa_default_parameters.shrake_rupley_n_points,
        'n-slices'     : freesasa_default_parameters.lee_richards_n_slices,
        'n-threads'    : freesasa_default_parameters.n_threads
    }

    def __init__(self,param=None):
        """
        Initializes Parameters object.

        Args:
              param (dict): optional argument to specify parameter-values,
                  see :py:attr:`.Parameters.defaultParameters`.

        Raises:
              AssertionError: Invalid parameter values supplied
        """
        self._c_param = freesasa_default_parameters
        if param != None:
            if 'algorithm' in param:    self.setAlgorithm(param['algorithm'])
            if 'probe-radius' in param: self.setProbeRadius(param['probe-radius'])
            if 'n-points' in param:     self.setNPoints(param['n-points'])
            if 'n-slices' in param:     self.setNSlices(param['n-slices'])
            if 'n-threads' in param:    self.setNThreads(param['n-threads'])
            unknownKeys = []
            for key in param:
                if not key in self.defaultParameters:
                    unknownKeys.append(key)
                if len(unknownKeys) > 0:
                    raise AssertionError('Key(s): ',unknownKeys,', unknown')

    def setAlgorithm(self,alg):
        """
        Set algorithm.

        Args:
             alg (str): algorithm name, only allowed values are
             :py:data:`freesasa.ShrakeRupley` and :py:data:`freesasa.LeeRichards`

        Raises:
             AssertionError: unknown algorithm specified
        """
        if alg == ShrakeRupley:
            self._c_param.alg = FREESASA_SHRAKE_RUPLEY
        elif alg == LeeRichards:
            self._c_param.alg = FREESASA_LEE_RICHARDS
        else:
            raise AssertionError("Algorithm '%s' is unknown" % alg)

    def algorithm(self):
        """
        Get algorithm.

        Returns:
            str: Name of algorithm
        """
        if self._c_param.alg == FREESASA_SHRAKE_RUPLEY:
            return ShrakeRupley
        if self._c_param.alg == FREESASA_LEE_RICHARDS:
            return LeeRichards
        raise Exception("No algorithm specified, shouldn't be possible")

    def setProbeRadius(self,r):
        """
        Set probe radius.

        Args:
            r (float): probe radius in Å (>= 0)

        Raises:
            AssertionError: r < 0
        """
        assert(r >= 0)
        self._c_param.probe_radius = r

    def probeRadius(self):
        """
        Get probe radius.

        Returns:
             float: Probe radius in Å
        """
        return self._c_param.probe_radius

    def setNPoints(self,n):
        """
        Set number of test points in Shrake & Rupley algorithm.

        Args:
            n (int): Number of points (> 0).

        Raises:
            AssertionError: n <= 0.
        """
        assert(n > 0)
        self._c_param.shrake_rupley_n_points = n

    def nPoints(self):
        """
        Get number of test points in Shrake & Rupley algorithm.

        Returns:
            int: Number of points.
        """
        return self._c_param.shrake_rupley_n_points

    def setNSlices(self,n):
        """
        Set the number of slices per atom in Lee & Richards algorithm.

        Args:
            n (int): Number of slices (> 0)

        Raises:
            AssertionError: n <= 0
        """
        assert(n> 0)
        self._c_param.lee_richards_n_slices = n

    def nSlices(self):
        """
        Get the number of slices per atom in Lee & Richards algorithm.

        Returns:
            int: Number of slices.
        """
        return self._c_param.lee_richards_n_slices

    def setNThreads(self,n):
        """
        Set the number of threads to use in calculations.

        Args:
            n (int): Number of points (> 0)

        Raises:
            AssertionError: n <= 0
        """
        assert(n>0)
        self._c_param.n_threads = n

    def nThreads(self):
        """
        Get the number of threads to use in calculations.

        Returns:
            int: Number of threads.
        """
        return self._c_param.n_threads

    # not pretty, but only way I've found to pass pointers around
    def _get_address(self, size_t ptr2ptr):
        cdef freesasa_parameters **p = <freesasa_parameters**> ptr2ptr
        p[0] = &self._c_param

class ResidueArea:
    """
    Stores absolute and relative areas for a residue

    Attributes:
        residueType (str): Type of Residue
        residueNumber (str): Residue number
        hasRelativeAreas (bool): False if there was noe reference area to calculate relative areas from

        total (float): Total SASA of residue
        polar (float): Polar SASA
        apolar (float): Apolar SASA
        mainChain (float): Main chain SASA
        sideChain (float): Side chain SASA

        relativeTotal (float): Relative total SASA
        relativePolar (float): Relative polar SASA
        relativeApolar (float): Relative Apolar SASA
        relativeMainChain (float): Relative main chain SASA
        relativeSideChain (float): Relative side chain SASA
    """

    residueType = ""
    residueNumber = ""
    hasRelativeAreas = False

    total = 0
    polar = 0
    apolar = 0
    mainChain = 0
    sideChain = 0

    relativeTotal = 0
    relativePolar = 0
    relativeApolar = 0
    relativeMainChain = 0
    relativeSideChain = 0


cdef class Result:
    """
    Stores results from SASA calculation.

    The type of object returned by :py:func:`freesasa.calc()`,
    not intended to be used outside of that context.
    """

    cdef freesasa_result* _c_result
    cdef freesasa_node* _c_root_node
    cdef freesasa_structure* _c_structure

    ## The constructor
    def __init__ (self):
        self._c_result = NULL
        self._c_root_node = NULL
        self._c_structure = NULL

    ## The destructor
    def __dealloc__(self):
        if self._c_result is not NULL:
            freesasa_result_free(self._c_result)
        if self._c_root_node is not NULL:
            freesasa_node_free(self._c_root_node)

    def nAtoms(self):
        """
        Number of atoms in the results.

        Returns:
            int: Number of atoms.
        """
        if self._c_result is not NULL:
            return self._c_result.n_atoms
        return 0

    def totalArea(self):
        """
        Total SASA.

        Returns:
            The total area in Å^2.
        Raises:
            AssertionError: If no results have been associated with the object.
        """
        assert(self._c_result is not NULL)
        return self._c_result.total

    def atomArea(self,i):
        """
        SASA for a given atom.

        Args:
            i (int): index of atom.

        Returns:
            float: SASA of atom i in Å^2.

        Raise:
            AssertionError: If no results have been associated
                      with the object or if index is out of bounds
        """
        assert(self._c_result is not NULL)
        assert(i < self._c_result.n_atoms)
        return self._c_result.sasa[i]

    def residueAreas(self):
        """
        Get SASA for all residues including relative areas if available for the
        classifier used.

        Returns dictionary of results where first dimension is chain label and
        the second dimension residue number. I.e. ``result["A"]["5"]`` gives the
        :py:class:`freesasa.ResidueArea` of residue number 5 in chain A.

        Relative areas are normalized to 1, but can be larger than one for
        residues in unusual conformations or at the ends of chains.

        Returns:
            dictionary

        Raise:
            AssertionError: If no results or structure has been associated
                 with the object.
        """
        assert(self._c_result is not NULL)
        assert(self._c_structure is not NULL)

        if (self._c_root_node == NULL):
            self._c_root_node = <freesasa_node*> freesasa_tree_init(self._c_result,
                                                                    self._c_structure,
                                                                    "Structure")
        cdef freesasa_node* result_node = <freesasa_node*> freesasa_node_children(self._c_root_node)
        cdef freesasa_node* structure = <freesasa_node*> freesasa_node_children(result_node)
        cdef freesasa_node* chain
        cdef freesasa_node* residue
        cdef freesasa_nodearea* c_area
        cdef freesasa_nodearea* c_ref_area

        result = {}

        chain = <freesasa_node*> freesasa_node_children(structure)
        while (chain != NULL):
            residue  = <freesasa_node*> freesasa_node_children(chain)
            chainLabel = freesasa_node_name(chain)
            result[chainLabel] = {}

            while (residue != NULL):
                c_area = <freesasa_nodearea*> freesasa_node_area(residue)
                c_ref_area = <freesasa_nodearea*> freesasa_node_residue_reference(residue)
                residueNumber = freesasa_node_residue_number(residue).strip()
                residueType = freesasa_node_name(residue).strip()

                area = ResidueArea()

                area.residueType = residueType
                area.residueNumber = residueNumber

                area.total = c_area.total
                area.mainChain = c_area.main_chain
                area.sideChain = c_area.side_chain
                area.polar = c_area.polar
                area.apolar = c_area.apolar

                if (c_ref_area is not NULL):
                    area.hasRelativeAreas = True
                    area.relativeTotal = self._safe_div(c_area.total, c_ref_area.total)
                    area.relativeMainChain = self._safe_div(c_area.main_chain, c_ref_area.main_chain)
                    area.relativeSideChain = self._safe_div(c_area.side_chain, c_ref_area.side_chain)
                    area.relativePolar = self._safe_div(c_area.polar, c_ref_area.polar)
                    area.relativeApolar = self._safe_div(c_area.apolar, c_ref_area.apolar)

                result[chainLabel][residueNumber.strip()] = area

                residue = <freesasa_node*> freesasa_node_next(residue)

            chain = <freesasa_node*> freesasa_node_next(chain)

        return result

    def _safe_div(self,a,b):
        try:
            return a/b
        except ZeroDivisionError:
            return float('nan')

    def _get_address(self, size_t ptr2ptr):
        cdef freesasa_result **p = <freesasa_result**> ptr2ptr
        p[0] = self._c_result


cdef class Classifier:
    """
    Assigns class and radius to atom by residue and atom name.

    Subclasses derived from :py:class:`.Classifier` can be used to define custom
    atomic radii and/or classes. Can also be initialized from
    config-files_ with a custom classifier.

    If initialized without arguments the default classifier is used.

    Derived classifiers must set the member :py:attr:`.purePython` to ``True``

    Residue names should be of the format ``"ALA"``, ``"ARG"``, etc.
    Atom names should be of the format ``"CA"``, ``"N"``, etc.
    """
    # this reference is used for classification
    cdef const freesasa_classifier *_c_classifier

    # if the classifier is read from a file we store it here,
    # with a reference in _c_classifier (for the sake of const-correctness)
    cdef freesasa_classifier *_dynamic_c_classifier

    # to be used by derived classes
    purePython = False

    def __init__ (self, fileName=None):
        """Constructor.

        If no file is provided the default classifier is used.

        Args:
            fileName (str): Name of file with classifier configuration.

        Raises:
            IOError:   Problem opening/reading file
            Exception: Problem parsing provided configuration or
                       initializing defaults
        """
        cdef FILE *config

        self._c_classifier = NULL
        self._dynamic_c_classifier = NULL

        if fileName is not None:
            config = fopen(fileName, 'rb')
            if config is NULL:
                raise IOError("File '%s' could not be opened." % fileName)
            self._dynamic_c_classifier = freesasa_classifier_from_file(config)
            fclose(config)
            self._c_classifier = self._dynamic_c_classifier;
            if self._c_classifier is NULL:
                raise Exception("Error parsing configuration in '%s'." % fileName)

        else:
            self._c_classifier = &freesasa_default_classifier

    # The destructor
    def __dealloc__(self):
        if (self._isCClassifier()):
            freesasa_classifier_free(self._dynamic_c_classifier)

    @staticmethod
    def getStandardClassifier(type):
        """
        Get a standard classifier (ProtOr, OONS or NACCESS)

        Args:
            type (str): The type, can have values ``'protor'``, ``'oons'`` or ``'naccess'``

        Returns:
            :py:class:`.Classifier`: The requested classifier

        Raises:
            Exception: If type not recognized
        """
        classifier = Classifier()
        if type == 'naccess':
            classifier._c_classifier = &freesasa_naccess_classifier
        elif type == 'oons':
            classifier._c_classifier = &freesasa_oons_classifier
        elif type == 'protor':
            classifier._c_classifier = &freesasa_protor_classifier
        else:
            raise Exception("Uknown classifier '%s'" % type)
        return classifier

    # This is used internally to determine if a Classifier wraps a C
    # classifier or not (necessary when generating structures)
    # returns Boolean
    def _isCClassifier(self):
        return not self.purePython

    def classify(self, residueName, atomName):
        """Class of atom.

        Depending on the configuration these classes can be
        anything, but typically they will be ``"Polar"`` and ``"Apolar"``.
        Unrecognized atoms will get the class ``"Unknown"``.

        Args:
            residueName (str): Residue name (`"ALA"`, `"ARG"`,...).
            atomName (str): Atom name (`"CA"`, `"C"`,...).

        Returns:
            str: Class name
        """
        classIndex = freesasa_classifier_class(self._c_classifier, residueName, atomName)
        return freesasa_classifier_class2str(classIndex)

    def radius(self,residueName,atomName):
        """Radius of atom.

        This allows the classifier to be used to calculate the atomic
        radii used in calculations. Unknown atoms will get a negative
        radius.

        Args:
            residueName (str): Residue name (`"ALA"`, `"ARG"`, ...).
            atomName (str): Atom name (`"CA"`, `"C"`, ...).

        Returns:
            float: The radius in Å.
        """
        return freesasa_classifier_radius(self._c_classifier, residueName, atomName)

    # the address obtained is a pointer to const
    def _get_address(self, size_t ptr2ptr):
        cdef freesasa_classifier **p = <freesasa_classifier**> ptr2ptr
        p[0] = <freesasa_classifier*>self._c_classifier # const cast

cdef class Structure:
    """
    Represents a protein structure, including its atomic radii.

    Initialized from PDB-file. Calculates atomic radii using default
    classifier, or custom one provided as argument to initalizer

    Since it is intended to be a static structure the word 'get' is
    omitted in the getter-functions.

    The default options are:
    ::
        defaultOptions = {
          'hetatm' : False,
          'hydrogen' : False,
          'join-models' : False,
          'skip-unknown' : False,
          'halt-at-unknown' : False
          }

    Attributes:
          defaultOptions: Default options for reading structure from PDB.
              By default ignore HETATM, Hydrogens, only use first
              model. For unknown atoms try to guess the radius, if
              this fails, assign radius 0 (to allow changing the
              radius later).

    """
    cdef freesasa_structure* _c_structure

    defaultOptions = {
          'hetatm' : False,
          'hydrogen' : False,
          'join-models' : False,
          'skip-unknown' : False,
          'halt-at-unknown' : False
          }

    defaultStructureArrayOptions = {
          'hetatm' : False,
          'hydrogen' : False,
          'separate-chains' : True,
          'separate-models' : False
    }

    def __init__(self,fileName=None,classifier=None,
                 options = defaultOptions):
        """
        Constructor

        If PDB file is provided, the structure will be constructed
        based on the file. If not, this simply initializes an empty
        structure and the other arguments are ignored. In this case
        atoms will have to be added manually using addAtom().

        Args:
            fileName (str): PDB file (if `None` empty structure generated).
            classifier: An optional :py:class:`.Classifier` to calculate atomic
                radii, uses default if none provided
            options (dict): specify which atoms and models to include, default is
                :py:attr:`.Structure.defaultOptions`

        Raises:
            IOError: Problem opening/reading file.
            Exception: Problem parsing PDB file or calculating
                atomic radii.
            Exception: If option 'halt-at-unknown' selected and
                unknown atom encountered.
        """

        self._c_structure = NULL
        cdef freesasa_classifier *c = NULL
        if classifier is None:
            classifier = Classifier()
        if classifier._isCClassifier():
            classifier._get_address(<size_t>&c)

        if fileName is None:
            self._c_structure = freesasa_structure_new()
            return
        cdef FILE *input
        input = fopen(fileName,'rb')
        if input is NULL:
            raise IOError("File '%s' could not be opened." % fileName)
        structure_options = Structure._get_structure_options(options)

        if not classifier._isCClassifier(): # supress warnings
            setVerbosity(silent)

        self._c_structure = freesasa_structure_from_pdb(input, c, structure_options)

        if not classifier._isCClassifier():
            setVerbosity(normal)

        fclose(input)

        if self._c_structure is NULL:
            raise Exception("Error reading '%s'." % fileName)

        # for pure Python classifiers we use the default
        # classifier above to initialize the structure and then
        # reassign radii using the provided classifier here
        if (not classifier._isCClassifier()):
            self.setRadiiWithClassifier(classifier)


    def addAtom(self, atomName, residueName, residueNumber, chainLabel, x, y, z):
        """
        Add atom to structure.

        This function is meant to be used if the structure was not
        initialized from a PDB. Default radii will be assigned to each
        atom. This can be overriden by calling
        :py:meth:`.Structure.setRadiiWithClassifier()` afterwards.

        There are no restraints on string lengths for the arguments, but
        the atom won't be added if the default classifier doesn't
        recognize the atom and also cannot deduce its element from the
        atom name.

        Args:
            atomName (str): atom name (e.g. `"CA"`)
            residueName (str): residue name (e.g. `"ALA"`)
            residueNumber (str or int): residue number (e.g. `'12'`)
                or integer. Some PDBs have residue-numbers that aren't
                regular numbers. Therefore treated as a string primarily.
            chainLabel (str): 1-character string with chain label (e.g. 'A')
                x,y,z (float): coordinates

        Raises:
            Exception: Residue-number invalid
        """
        if (type(residueNumber) is str):
            resnum = residueNumber
        elif (type(residueNumber) is int):
            resnum = "%d" % residueNumber
        else:
            raise Exception("Residue-number invalid, must be either string or number")
        cdef const char *label = chainLabel
        ret = freesasa_structure_add_atom(self._c_structure, atomName,
                                          residueName, resnum, label[0],
                                          x, y, z)
        assert(ret != FREESASA_FAIL)

    def setRadiiWithClassifier(self,classifier):
        """
        Assign radii to atoms in structure using a classifier.

        Args:
            classifier: A :py:class:`.Classifier` to use to calculate radii.

        Raises:
            AssertionError: if structure not properly initialized
        """
        assert(self._c_structure is not NULL)
        n = self.nAtoms()
        r = []
        for i in range(0,n):
            r.append(classifier.radius(self.residueName(i), self.atomName(i)))
        self.setRadii(r)

    def setRadii(self,radiusArray):
        """
        Set atomic radii from an array

        Args:
            radiusArray (list): Array of atomic radii in Ångström, should
                have nAtoms() elements.
        Raises:
            AssertionError: if radiusArray has wrong dimension, structure
                not properly initialized, or if the array contains
                negative radii (not properly classified?)
        """
        assert(self._c_structure is not NULL)
        n = self.nAtoms()
        assert len(radiusArray) == n
        cdef double *r = <double *>malloc(sizeof(double)*n)
        assert(r is not NULL)
        for i in range(0,n):
            r[i] = radiusArray[i]
            assert(r[i] >= 0), "Error: Radius array is <= 0 for the residue: " + self.residueName(i) + " ,atom: " + self.atomName(i)
        freesasa_structure_set_radius(self._c_structure, r)

    def nAtoms(self):
        """
        Number of atoms.

        Returns:
            int: Number of atoms

        Raises:
            AssertionError: if not properly initialized
        """
        assert(self._c_structure is not NULL)
        return freesasa_structure_n(self._c_structure)

    def radius(self,i):
        """
        Radius of atom.

        Args:
            i (int): Index of atom.

        Returns:
            float: Radius in Å.

        Raises:
            AssertionError: if index out of bounds, object not properly initalized.
        """
        assert(i >= 0 and i < self.nAtoms())
        assert(self._c_structure is not NULL)
        cdef const double *r = freesasa_structure_radius(self._c_structure)
        assert(r is not NULL)
        return r[i]

    def setRadius(self, atomIndex, radius):
        """
        Set radius for a given atom

        Args:
            atomIndex (int): Index of atom
            radius (float): Value of radius

        Raises:
            AssertionError: if index out of bounds, radius
                negative, or structure not properly initialized
        """
        assert(self._c_structure is not NULL)
        assert(atomIndex >= 0 and atomIndex < self.nAtoms())
        assert(radius >= 0)
        freesasa_structure_atom_set_radius(self._c_structure, atomIndex, radius)

    def atomName(self,i):
        """
        Get atom name

        Args:
            i (int): Atom index.

        Returns:
            str: Atom name as 4-character string.

        Raises:
            AssertionError: if index out of range or Structure not properly initialized.
        """
        assert(i >= 0 and i < self.nAtoms())
        assert(self._c_structure is not NULL)
        return freesasa_structure_atom_name(self._c_structure,i)

    def residueName(self,i):
        """
        Get residue name of given atom.

        Args:
            i (int): Atom index.

        Returns:
            str: Residue name as 3-character string.

        Raises:
            AssertionError: if index out of range or Structure not properly initialized
        """
        assert(i >= 0 and i < self.nAtoms())
        assert(self._c_structure is not NULL)
        return freesasa_structure_atom_res_name(self._c_structure,i)

    def residueNumber(self,i):
        """
        Get residue number for given atom.

        Residue number will include the insertion code if there is one.

        Args:
            i (int): Atom index.

        Returns:
            str: Residue number as 5-character string (last character is either whitespace or insertion code)

        Raises:
            AssertionError: if index out of range or Structure not properly initialized
        """
        assert(i >= 0 and i < self.nAtoms())
        assert(self._c_structure is not NULL)
        return freesasa_structure_atom_res_number(self._c_structure,i)

    def chainLabel(self,i):
        """
        Get chain label for given atom.

        Args:
            i (int): Atom index.

        Returns:
            str: Chain label as 1-character string.

        Raises:
            AssertionError: if index out of range or Structure not properly initialized
        """
        assert(i >= 0 and i < self.nAtoms())
        assert(self._c_structure is not NULL)
        cdef char label[2]
        label[0] = freesasa_structure_atom_chain(self._c_structure,i)
        label[1] = '\0'
        return label

    def coord(self, i):
        """
        Get coordinates of given atom.

        Args:
            i (int): Atom index.

        Returns:
            list: array of x, y, and z coordinates

        Raises:
            AssertionError: if index out of range or Structure not properly initialized
        """
        assert(i >= 0 and i < self.nAtoms())
        assert(self._c_structure is not NULL)
        cdef const double *coord = freesasa_structure_coord_array(self._c_structure)
        return [coord[3*i], coord[3*i+1], coord[3*i+2]]

    @staticmethod
    def _get_structure_options(param):
        options = 0

        # check validity of options
        knownOptions = {'hetatm','hydrogen','join-models','separate-models',
                        'separate-chains','skip-unknown','halt-at-unknown'}
        unknownOptions = []
        for key in param:
            if not key in knownOptions:
                unknownOptions.append(key)
        if len(unknownOptions) > 0:
            raise AssertionError("Option(s): ",unknownOptions," unknown.")

        # calculate bitfield
        if 'hetatm' in param and param['hetatm']:
            options |= FREESASA_INCLUDE_HETATM
        if 'hydrogen' in param and param['hydrogen']:
            options |= FREESASA_INCLUDE_HYDROGEN
        if 'join-models' in param and param['join-models']:
            options |= FREESASA_JOIN_MODELS
        if 'separate-models' in param and param['separate-models']:
            options |= FREESASA_SEPARATE_MODELS
        if 'separate-chains' in param and param['separate-chains']:
            options |= FREESASA_SEPARATE_CHAINS
        if 'skip-unknown' in param and param['skip-unknown']:
            options |= FREESASA_SKIP_UNKNOWN
        if 'halt-at-unknown' in param and param['halt-at-unknown']:
            options |= FREESASA_HALT_AT_UNKNOWN
        return options

    def _get_address(self, size_t ptr2ptr):
        cdef freesasa_structure **p = <freesasa_structure**> ptr2ptr
        p[0] = self._c_structure

    def _set_address(self, size_t ptr2ptr):
        cdef freesasa_structure **p = <freesasa_structure**> ptr2ptr
        self._c_structure = p[0]

    ## The destructor
    def __dealloc__(self):
        if self._c_structure is not NULL:
            freesasa_structure_free(self._c_structure)


def structureArray(fileName,
                   options = Structure.defaultStructureArrayOptions,
                   classifier = None):
    """
    Create array of structures from PDB file.

    Split PDB file into several structures by either by treating
    chains separately, by treating each MODEL as a separate
    structure, or both.

    Args:
        fileName (str): The PDB file.
        options (dict): Specification for how to read the PDB-file
            (see :py:attr:`.Structure.defaultStructureArrayOptions` for
            options and default value).
        classifier: :py:class:`.Classifier` to assign atoms radii, default is used
            if none specified.

    Returns:
        list: An array of :py:class:`.Structure`

    Raises:
        AssertionError: if `fileName` is None
        AssertionError: if an option value is not recognized
        AssertionError: if neither of the options `'separate-chains'`
            and `'separate-models'` are specified.
        IOError: if can't open file
        Exception: if there are problems parsing the input
    """

    assert fileName is not None
    # we need to have at least one of these
    assert(('separate-chains' in options and options['separate-chains'] is True)
           or ('separate-models' in options and options['separate-models'] is True))
    structure_options = Structure._get_structure_options(options)
    cdef FILE *input
    input = fopen(fileName,'rb')
    if input is NULL:
        raise IOError("File '%s' could not be opened." % fileName)
    cdef int n

    verbosity = getVerbosity()

    if classifier is not None:
        setVerbosity(silent)
    cdef freesasa_structure** sArray = freesasa_structure_array(input,&n,NULL,structure_options)
    fclose(input)

    if classifier is not None:
        setVerbosity(verbosity)

    if sArray is NULL:
        raise Exception("Problems reading structures in '%s'." % fileName)
    structures = []
    for i in range(0,n):
        structures.append(Structure())
        structures[-1]._set_address(<size_t> &sArray[i])
        if classifier is not None:
            structures[-1].setRadiiWithClassifier(classifier)
    free(sArray)
    return structures


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
    result._c_structure = <freesasa_structure*> s

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

def structureFromBioPDB(bioPDBStructure, classifier=None, options = Structure.defaultOptions):
    """
    Create a freesasa structure from a Bio.PDB structure

    Experimental, not thorougly tested yet.
    Structures generated this way will not preserve whitespace in residue numbers, etc,
    as in :py:class:`.Structure`.

    Args:
        bioPDBStructure: a `Bio.PDB` structure
        classifier: an optional :py:class:`.Classifier` to specify atomic radii
        options (dict): Options supported are `'hetatm'`, `'skip-unknown'` and `'halt-at-unknown'`

    Returns:
        :py:class:`.Structure`: The structure

    Raises:
        Exception: if option 'halt-at-unknown' is selected and
            unknown atoms are encountered. Passes on exceptions from
            :py:meth:`.Structure.addAtom()` and
            :py:meth:`.Structure.setRadiiWithClassifier()`.
    """
    structure = Structure()
    if (classifier is None):
        classifier = Classifier()
    optbitfield = Structure._get_structure_options(options)

    atoms = bioPDBStructure.get_atoms()

    for a in atoms:
        r = a.get_parent()
        hetflag, resseq, icode = r.get_id()
        resname = r.get_resname()

        if (hetflag is not ' ' and not (optbitfield & FREESASA_INCLUDE_HETATM)):
            continue

        c = r.get_parent()
        v = a.get_vector()
        if (icode):
            resseq = str(resseq) + str(icode)

        if (classifier.classify(resname, a.get_fullname()) is 'Unknown'):
            if (optbitfield & FREESASA_SKIP_UNKNOWN):
                continue
            if (optbitfield & FREESASA_HALT_AT_UNKNOWN):
                raise Exception("Halting at unknown atom")

        structure.addAtom(a.get_fullname(), r.get_resname(), resseq, c.get_id(),
                          v[0], v[1], v[2])

    structure.setRadiiWithClassifier(classifier)
    return structure

def calcBioPDB(bioPDBStructure, parameters = Parameters(),
               classifier = None, options = Structure.defaultOptions):
    """
    Calc SASA from `BioPython` PDB structure.

    Usage::

        result, sasa_classes = calcBioPDB(structure, ...)

    Experimental, not thorougly tested yet

    Args:
        bioPDBStructure: A `Bio.PDB` structure
        parameters: A :py:class:`.Parameters` object (uses default if none specified)
        classifier: A :py:class:`.Classifier` object (uses default if none specified)
        options (dict): Options supported are 'hetatm', 'skip-unknown' and 'halt-at-unknown'
            (uses :py:attr:`.Structure.defaultOptions` if none specified

    Returns:
        A :py:class:`.Result` object and a dictionary with classes
        defined by the classifier and associated areas

    Raises:
        Exception: if unknown atom is encountered and the option
            'halt-at-unknown' is active. Passes on exceptions from
            :py:func:`.calc()`, :py:func:`.classifyResults()` and
            :py:func:`.structureFromBioPDB()`.
    """
    structure = structureFromBioPDB(bioPDBStructure, classifier, options)
    result = calc(structure, parameters)
    sasa_classes = classifyResults(result, structure, classifier)
    return result, sasa_classes
