.. toctree::
   :maxdepth: 2

   intro
   functions
   classes

FreeSASA Python Module
======================

The module provides Python bindings for the `FreeSASA C Library <https://github.com/mittinatten/freesasa>`_.
It works with Python 2.7 and 3.5+, on Linux, Mac OS X and Windows. The source is available as a PyPi source
distribution and on `GitHub <https://github.com/freesasa/freesasa-python>`_.

Install the FreeSASA Python Module by

.. code::

    pip install freesasa


Developers can clone the library, and then build the module by the following

.. code::

   git clone https://github.com/freesasa/freesasa-python.git
   git submodule update --init
   python setyp.py build

Tests are run by

.. code::

    python setup.py test
