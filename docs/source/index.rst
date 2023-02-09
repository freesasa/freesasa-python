.. toctree::
   :maxdepth: 2

   intro
   functions
   classes

FreeSASA Python Module
======================

The module provides Python bindings for the `FreeSASA C Library <https://github.com/mittinatten/freesasa>`_.
Python 3.7+, on Linux, Mac OS X and Windows are officially supported (it will probably still run on older Python
versions if you build it from source, or use older PyPi packages).
The source is available as a PyPi source
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
