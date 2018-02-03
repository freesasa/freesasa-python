FreeSASA Python module
======================

The PyPi module has FreeSASA as a submodule, and uses a separate
setup.py and config.h to avoid dependence on autotools.

After cloning, build the module by the following
~~~~sh
git submodule init
git submodule update
python setup.py build
~~~~

Tests can be run using
~~~~sh
PYTHONPATH=./build/<architecture-dependent build directory>/:$PYTHONPATH python test.py
~~~~
