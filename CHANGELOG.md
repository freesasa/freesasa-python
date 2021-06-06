# Changelog

# 2.2.0 Pre-release

- Add `Result.write_pdb()`
- Add `Structre.addAtoms()`

# 2.1.0

- Added changelog
- Can access absolute and relative SASA for individual residues through `Result.residueAreas()`
- Can set options and classifier for a `Structure` initiated without an input file for later
  use in `Structure.addAtom()`
- Only build PyPi packages for Python 3.6+ (can still be built from source for older versions)
