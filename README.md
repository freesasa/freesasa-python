FreeSASA Python module
======================
[![Appveyor build status](https://ci.appveyor.com/api/projects/status/nyo51pv2ufj2yhcj/branch/master?svg=true)](https://ci.appveyor.com/project/mittinatten/freesasa-python/branch/master)
[![Travis build status](https://travis-ci.org/freesasa/freesasa-python.svg?branch=master)](https://travis-ci.org/freesasa/freesasa-python)

The module provides Python bindings for the [FreeSASA C Library](https://github.com/mittinatten/freesasa).
It works with Python 2.7 and 3.4+, on Linux, Mac OS X and Windows. Documentation can be found at http://freesasa.github.io/python/.

Install the module by
~~~~sh
pip install freesasa
~~~~

Developers can clone the library, and then build the module by the following
~~~~sh
git submodule update --init
USE_CYTHON=1 python setup.py build
~~~~

Tests can be run using
~~~~sh
python setup.py test
~~~~
