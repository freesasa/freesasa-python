.. currentmodule:: freesasa

Introduction
============

This package provides Python bindings for the `FreeSASA C Library
<https://github.com/mittinatten/freesasa>`_.

It can be installed using

.. code::

   pip install freesasa

Binaries are available for Python 2.7, 3.5, 3.6 and 3.7 for Mac OS X
and Windows, in addition to the source distribution.


Basic calculations
------------------

Using defaults everywhere a simple calculation can be carried out as
follows (assuming the file ``1ubq.pdb`` is available)

.. code:: python

    import freesasa

    structure = freesasa.Structure("1ubq.pdb")
    result = freesasa.calc(structure)
    area_classes = freesasa.classifyResults(result, structure)

    print "Total : %.2f A2" % result.totalArea()
    for key in area_classes:
        print key, ": %.2f A2" % area_classes[key]

Which would give the following output

.. code::

    Total : 4804.06 A2
    Polar : 2504.22 A2
    Apolar : 2299.84 A2

The following does a high precision L&R calculation

.. code:: python

    result = freesasa.calc(structure,
                           freesasa.Parameters({'algorithm' : freesasa.LeeRichards,
                                                'n-slices' : 100}))

Using the results from a calculation we can also integrate SASA over a selection of
atoms, using a subset of the Pymol `selection syntax`_:

.. _selection syntax: http://freesasa.github.io/doxygen/Selection.html

.. code:: python

    selections = freesasa.selectArea(('alanine, resn ala', 'r1_10, resi 1-10'),
                                     structure, result)
    for key in selections:
        print key, ": %.2f A2" % selections[key]

which gives the output

.. code::

    alanine : 120.08 A2
    r1_10 : 634.31 A2

Customizing atom classification
-------------------------------

This uses the NACCESS parameters (the file ``naccess.config`` is
available in the ``share/`` directory of the repository).

.. code:: python

    classifier = freesasa.Classifier("naccess.config")
    structure = freesasa.Structure("1ubq.pdb", classifier)
    result = freesasa.calc(structure)
    area_classes = freesasa.classifyResults(result, structure, classifier)

Classification can be customized also by extending the :py:class:`.Classifier`
interface. The code below is an illustration of a classifier that
classes nitrogens separately, and assigns radii based on element only
(and crudely).

.. code:: python

    import freesasa
    import re

    class DerivedClassifier(freesasa.Classifier):
        # this must be set explicitly in all derived classifiers
        purePython = True

        def classify(self, residueName, atomName):
            if re.match('\s*N', atomName):
                return 'Nitrogen'
            return 'Not-nitrogen'

        def radius(self, residueName, atomName):
            if re.match('\s*N',atomName): # Nitrogen
                return 1.6
            if re.match('\s*C',atomName): # Carbon
                return 1.7
            if re.match('\s*O',atomName): # Oxygen
                return 1.4
            if re.match('\s*S',atomName): # Sulfur
                return 1.8
            return 0;                     # everything else (Hydrogen, etc)

    classifier = DerivedClassifier()

    # use the DerivedClassifier to calculate atom radii
    structure = freesasa.Structure("1ubq.pdb", classifier)
    result = freesasa.calc(structure)

    # use the DerivedClassifier to classify atoms
    area_classes = freesasa.classifyResults(result,structure,classifier)

Of course, this example is somewhat contrived, if we only want the
integrated area of nitrogen atoms, the simpler choice would be

.. code:: python

    selection = freesasa.selectArea('nitrogen, symbol n', structure, result)


However, extending :py:class:`.Classifier`, as illustrated above, allows
classification to arbitrary complexity and also lets us redefine the
radii used in the calculation.

Bio.PDB
-------

FreeSASA can also calculate the SASA of a ``Bio.PDB`` structure from BioPython

.. code:: python

    from Bio.PDB import PDBParser
    parser = PDBParser()
    structure = parser.get_structure("Ubiquitin", "1ubq.pdb")
    result, sasa_classes = freesasa.calcBioPDB(structure)

If one needs more control over the analysis the structure can be
converted to a :py:class:`.Structure` using :py:func:`.structureFromBioPDB()`
and the calculation can be performed the normal way using this
structure.
