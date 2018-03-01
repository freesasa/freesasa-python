FreeSASA Python module
======================
[![Appveyor build status](https://ci.appveyor.com/api/projects/status/nyo51pv2ufj2yhcj/branch/master?svg=true)](https://ci.appveyor.com/project/mittinatten/freesasa-python/branch/master)
[![Travis build status](https://travis-ci.org/freesasa/freesasa-python.svg?branch=master)](https://travis-ci.org/freesasa/freesasa-python)

The module provides Python bindings for the [FreeSASA C Library](https://github.com/mittinatten/freesasa).
It works with Python 2.7 and 3.4+, on Linux, Mac OS X and Windows. Still in pre-release until the documentation is ready, for now see http://freesasa.github.io/doxygen/Python.html.

It can be installed using
~~~~sh
pip install freesasa
~~~~
Adding the flag `--pre` is necessary on some platforms to install it before the stable release is available.

Developers can clone the library, and then build the module by the following
~~~~sh
git submodule update --init
python setup.py build
~~~~

Tests can be run using
~~~~sh
python setup.py test
~~~~
