import string
import warnings
from cfreesasa cimport *
from libc.stdio cimport FILE, fopen, fclose
from libc.stdlib cimport malloc, realloc

cdef class Structure:
    """
    Represents a protein structure, including its atomic radii.

    Initialized from PDB-file. Calculates atomic radii using default
    classifier, or custom one provided as argument to initalizer.

    Default options are ::

        Structure.defaultOptions = {
            'hetatm' : False,          # False: skip HETATM
                                       # True: include HETATM

            'hydrogen' : False,        # False: ignore hydrogens
                                       # True: include hydrogens

            'join-models' : False,     # False: Only use the first MODEL
                                       # True: Include all MODELs

            'skip-unknown' : False,    # False: Guess radius for unknown atoms
                                       #     based on element
                                       # True: Skip unknown atoms

            'halt-at-unknown' : False  # False: set radius for unknown atoms,
                                       #    that can not be guessed to 0.
                                       # True: Throw exception on unknown atoms.
        }

    Attributes:
        defaultOptions:  Default options for reading structure from PDB.

    """
    cdef freesasa_structure* _c_structure
    cdef const freesasa_classifier* _c_classifier
    cdef int _c_options

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

    def __init__(self, fileName=None, classifier=None,
                 options = defaultOptions):
        """
        Constructor

        If a PDB file is provided, the structure will be constructed
        based on the file. If not, this simply initializes an empty
        structure with the given classifier and options. Atoms will then
        have to be added manually using `:py:meth:`.Structure.addAtom()`.

        Args:
            fileName (str): PDB file (if `None` empty structure generated).
            classifier: An optional :py:class:`.Classifier` to calculate atomic
                radii, uses default if none provided.
                This classifier will also be used in calls to :py:meth:`.Structure.addAtom()`
                but only if it's the default classifier, one of the standard
                classifiers from :py:meth:`.Classifier.getStandardClassifier()`,
                or defined by a config-file (i.e. if it uses the underlying
                C API).
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
        self._c_classifier = NULL

        if classifier is None:
            classifier = Classifier()
        if classifier._isCClassifier():
            classifier._get_address(<size_t>&self._c_classifier)

        Structure._validate_options(options)
        self._c_options = Structure._get_bitfield_from_options(options)

        if fileName is None:
            self._c_structure = freesasa_structure_new()
        else:
            self._initFromFile(fileName, classifier)

    def _initFromFile(self, fileName, classifier):
        cdef FILE *input
        input = fopen(fileName,'rb')
        if input is NULL:
            raise IOError("File '%s' could not be opened." % fileName)

        if not classifier._isCClassifier(): # supress warnings
            setVerbosity(silent)

        self._c_structure = freesasa_structure_from_pdb(input, self._c_classifier, self._c_options)

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
        initialized from a PDB. The options and classifier passed to
        the constructor for the :py:class:`.Structure` will be used
        (see the documentation of the constructor for restrictions).
        The radii set by the classifier can be overriden by calling
        :py:meth:`.Structure.setRadiiWithClassifier()` afterwards.

        There are no restraints on string lengths for the arguments, but
        the atom won't be added if the classifier doesn't
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
            AssertionError:
        """
        if (type(residueNumber) is str):
            resnum = residueNumber
        elif (type(residueNumber) is int):
            resnum = "%d" % residueNumber
        else:
            raise Exception("Residue-number invalid, must be either string or number")

        cdef const char *label = chainLabel

        ret = freesasa_structure_add_atom_wopt(
            self._c_structure, atomName,
            residueName, resnum, label[0],
            x, y, z,
            self._c_classifier, self._c_options)

        assert(ret != FREESASA_FAIL)

    def addAtoms(self, atomNames, residueNames, residueNumbers, chainLabels, xs, ys, zs):
        """
        Add multiple atoms to structure.

        Args:
            atomNames (list): list of atom name (e.g. `["CA"]`)
            residueNames (list): list of residue name (e.g. `["ALA"]`)
            residueNumbers (list): list of residue number (e.g. `['12']`)
                or integer. Some PDBs have residue-numbers that aren't
                regular numbers. Therefore treated as a string primarily.
            chainLabels (list): list of 1-character string with chain label (e.g. ['A'])
                xs,ys,zs (list): list of coordinates

        Raises:
            AssertionError: inconsistent size of input args
        """
        assert(len(set([len(atomNames), len(residueNames), len(residueNumbers), \
            len(chainLabels), len(xs), len(ys), len(zs)])) == 1), "Inconsistent size of input args"

        for i in range(len(atomNames)):
            self.addAtom(atomNames[i], residueNames[i], residueNumbers[i], \
                chainLabels[i], xs[i], ys[i], zs[i])

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
            assert(r[i] >= 0), "Error: Radius is <= 0 (" + str(r[i]) + ") for the residue: " + self.residueName(i) + ", atom: " + self.atomName(i)
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
    def _validate_options(param):
        # check validity of options
        knownOptions = {'hetatm','hydrogen','join-models','separate-models',
                        'separate-chains', 'chain-groups', 'skip-unknown', 'halt-at-unknown'}
        unknownOptions = set(param.keys()).difference(knownOptions)
        if len(unknownOptions) > 0:
            raise AssertionError("Option(s): ",unknownOptions," unknown.")

        # 'chain-groups' additional checks
        if param.get('chain-groups', False):
            if param.get('separate-chains', False):
                raise ValueError("'chain-groups' and 'separate-chains' cannot be set simultaneously.")
            allowed_chars = set(string.ascii_lowercase + string.ascii_uppercase + '+')
            if not all([character in allowed_chars for character in param['chain-groups']]):
                raise ValueError("'{}' is not a valid chain-groups selection.".format(param['chain-groups']))

    @staticmethod
    def _get_bitfield_from_options(param):
        options = 0
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
    chains separately and/or by treating each MODEL as a separate
    structure and/or grouping chains.

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
           or ('separate-models' in options and options['separate-models'] is True)
           or (options.get('chain-groups', False)))

    Structure._validate_options(options)
    structure_options = Structure._get_bitfield_from_options(options)
    cdef FILE *input
    input = fopen(fileName,'rb')
    if input is NULL:
        raise IOError("File '%s' could not be opened." % fileName)

    verbosity = getVerbosity()

    if classifier is not None:
        setVerbosity(silent)

    cdef int n
    cdef freesasa_structure** sArray
    if options.get('separate-chains', False) or options.get('separate-models', False):
        sArray = freesasa_structure_array(input,&n,NULL,structure_options)
    else:
        sArray = <freesasa_structure **> malloc(sizeof(freesasa_structure *))
        sArray[0] = freesasa_structure_from_pdb(input, NULL, structure_options)
        n = 1
    fclose(input)

    if sArray is NULL:
        raise Exception("Problems reading structures in '%s'." % fileName)
    if sArray[0] is NULL:
        free(sArray)
        raise Exception("Problems reading structures in '%s'." % fileName)

    cdef freesasa_structure* group
    if options.get('chain-groups', False):
        # Add groups from each model
        chain_groups = options['chain-groups'].split('+')
        n_groups = len(chain_groups)
        n_total = n * (1 + n_groups)
        sArray = <freesasa_structure**> realloc(sArray, n_total * sizeof(freesasa_structure*))
        if sArray is NULL:
            raise Exception("Out of memory when allocating '%i' structures from '%s'." % (n_total, fileName))
        for i in range(0, n):
            for j, chainID_group in enumerate(chain_groups):
                idx = j + (i * n_groups) + n
                group = freesasa_structure_get_chains(sArray[i], chainID_group)
                sArray[idx] = group
        n = n_total

    if classifier is not None:
        setVerbosity(verbosity)

    structures = []
    for i in range(0, n):
        structures.append(Structure())
        structures[-1]._set_address(<size_t> &sArray[i])
        if classifier is not None:
            structures[-1].setRadiiWithClassifier(classifier)
    free(sArray)

    return structures


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
    Structure._validate_options(options)
    optbitfield = Structure._get_bitfield_from_options(options)

    atoms = bioPDBStructure.get_atoms()

    if (not optbitfield & FREESASA_JOIN_MODELS):
        models = bioPDBStructure.get_models()
        atoms = next(models).get_atoms()

    for a in atoms:
        r = a.get_parent()
        s = r.get_parent()
        hetflag, resseq, icode = r.get_id()
        resname = r.get_resname()
        atomname = a.get_fullname()
        element = a.element

        if (hetflag is not ' ' and not (optbitfield & FREESASA_INCLUDE_HETATM)):
            continue
        if (((element == "H") or (element == "D")) and not (optbitfield & FREESASA_INCLUDE_HYDROGEN)):
            continue

        c = r.get_parent()
        v = a.get_vector()
        if (icode):
            resseq = str(resseq) + str(icode)

        if (classifier.classify(resname, atomname) is 'Unknown'):
            if (optbitfield & FREESASA_SKIP_UNKNOWN):
                continue
            if (optbitfield & FREESASA_HALT_AT_UNKNOWN):
                raise Exception("Halting at unknown atom")

        structure.addAtom(atomname, r.get_resname(), resseq, c.get_id(),
                          v[0], v[1], v[2])

    structure.setRadiiWithClassifier(classifier)
    return structure
