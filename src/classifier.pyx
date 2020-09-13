from cfreesasa cimport *
from libc.stdio cimport FILE, fopen, fclose

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
