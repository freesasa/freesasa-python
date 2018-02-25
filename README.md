FreeSASA Python module
======================
[![Build status](https://ci.appveyor.com/api/projects/status/nyo51pv2ufj2yhcj?svg=true)](https://ci.appveyor.com/project/mittinatten/freesasa-python)

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
