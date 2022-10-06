from setuptools import setup, Extension
import os
import sys
from glob import glob

USE_CYTHON = False

try:
    USE_CYTHON = os.environ['USE_CYTHON']
except KeyError:
    if not os.path.isfile(os.path.join("src", "freesasa.c")):
        sys.stderr.write("No C source detected, define environment variable USE_CYTHON to build from Cython source.\n")
        sys.exit()
    else:
        print ("Define environment variable USE_CYTHON to build from Cython source")


# not using wild cards because we're leaving out xml and json
sources = list(map(lambda file: os.path.join('lib', 'src', file),
              ["classifier.c",
               "classifier_protor.c", "classifier_oons.c", "classifier_naccess.c",
               "coord.c", "freesasa.c", "lexer.c", "log.c",
               "nb.c", "node.c", "parser.c",
               "pdb.c", "rsa.c", "sasa_lr.c", "sasa_sr.c",
               "selection.c", "structure.c",
               "util.c"]))


extensions = None

ext = '.pyx' if USE_CYTHON else '.c'
sources.append(os.path.join('src', 'freesasa' + ext))

compile_args=['-DHAVE_CONFIG_H']

if os.name == 'posix':
    compile_args.append('-std=gnu99')

extension_src = [
    Extension("freesasa", sources,
              language='c',
              include_dirs=[os.path.join('lib', 'src'), '.'],
              extra_compile_args = compile_args
	      )
]

if USE_CYTHON:
    from Cython.Build import cythonize, build_ext
    extensions = cythonize(extension_src)
else:
    extensions = extension_src

long_description = \
 "This module provides bindings for the FreeSASA C library. " + \
 "See http://freesasa.github.io/python/ for documentation."

setup(
    name='freesasa',
    version= '2.2.0a4',
    description='Calculate solvent accessible surface areas of proteins',
    long_description=long_description,
    author='Simon Mitternacht',
    url='http://freesasa.github.io/',
    license='MIT',
    ext_modules=extensions,
    keywords=['structural biology', 'proteins', 'bioinformatics'],
    headers=glob(os.path.join('lib', 'src', '*')),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        # Will also build for python 2.7 but not officialy supported
    ],
    setup_requires=['cython>=0.29.13'],
    test_suite='test'
)
