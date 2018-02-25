FreeSASA Python module
======================
[![Build status](https://ci.appveyor.com/api/projects/status/nyo51pv2ufj2yhcj?svg=true)](https://ci.appveyor.com/project/mittinatten/freesasa-python)

The module has the FreeSASA C source as a submodule. It works with Python 2.7 and 3.4+, on Linux, Mac OS X and Windows.

After cloning, build the module by the following
~~~~sh
git submodule update --init
python setup.py build
~~~~

Tests can be run using
~~~~sh
python setup.py test
~~~~
