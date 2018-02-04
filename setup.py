from setuptools import setup, Extension
from Cython.Build import cythonize, build_ext
from glob import glob

# not using wild cards because we're leaving out xml and json
sources = list(map(lambda file: "./lib/src/" + file,
              ["classifier.c",
               "classifier_protor.c", "classifier_oons.c", "classifier_naccess.c",
               "coord.c", "freesasa.c", "lexer.c", "log.c",
               "nb.c", "node.c", "parser.c",
               "pdb.c", "rsa.c", "sasa_lr.c", "sasa_sr.c",
               "selection.c", "structure.c",
               "util.c"]))

sources.append("./*.pyx")

defines = [
    "-DUSE_THREADS=1", "-DUSE_XML=0",
    "-DUSE_JSON=0", "-DUSE_CHECK=0",
    '-DPACKAGE="freesasa"',
    '-DPACKAGE_NAME="FreeSASA"',
    '-DPACKAGE_STRING="FreeSASA 2.0.2"',
    '-DPACKAGE_VERSION="2.0.2"'
]

extensions = [
    Extension("freesasa", sources,
              glob("./src/*.h"),
              language='c',
              extra_compile_args = defines
	      )
]

setup(
    name='freesasa',
    description='Calculate solvent accessible surface areas of proteins',
    version= '2.0.2',
    author='Simon Mitternacht',
    url='http://freesasa.github.io/',
    license='MIT',
    ext_modules=cythonize(extensions),
    keywords=['structural biology', 'proteins', 'bioinformatics'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
    ],
    setup_requires=['cython>=0.21']
)
