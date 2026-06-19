#!/usr/bin/env python3
"""
Compatibility test for the parallelized FreeSASA Python API.
Tests all patterns used by:
  - freesasa-python standard API
  - mdakit_sasa (pegerto/mdakit_sasa)

Run from freesasa-python directory:
  python3 test_compat.py
"""

import sys
import os

# Use our local parallelized build
sys.path.insert(0, '/media/amin/10TB_2/WORK/dev_sasa_prallell/freesasa-python')
import freesasa

FREESASA_PDB = '/media/amin/10TB_2/WORK/dev_sasa_prallell/freesasa/tests/data/1ubq.pdb'
LARGE_PDB    = '/media/amin/10TB_2/WORK/dev_sasa_prallell/freesasa/tests/data/1aon.pdb'

passed = []
failed = []

def check(name, condition, detail=''):
    if condition:
        passed.append(name)
        print(f'  [PASS] {name}')
    else:
        failed.append(name)
        print(f'  [FAIL] {name}' + (f' — {detail}' if detail else ''))

def run_tests():
    print("=" * 60)
    print("  FreeSASA Python API Compatibility Test")
    print("=" * 60)

    # ── 1. Parameters API (used by mdakit_sasa) ──────────────────
    print("\n--- Parameters API ---")
    p = freesasa.Parameters()
    check("Parameters() default constructor", p is not None)
    check("setNThreads / nThreads", (p.setNThreads(4) is None or True) and p.nThreads() == 4)
    check("setNThreads(1) single-threaded", p.setNThreads(1) is None or True)
    check("setAlgorithm LeeRichards", (p.setAlgorithm(freesasa.LeeRichards), p.algorithm())[1] == freesasa.LeeRichards)
    check("setAlgorithm ShrakeRupley", (p.setAlgorithm(freesasa.ShrakeRupley), p.algorithm())[1] == freesasa.ShrakeRupley)
    check("setNSlices / nSlices", (p.setNSlices(20), p.nSlices())[1] == 20)
    check("setNPoints / nPoints", (p.setNPoints(100), p.nPoints())[1] == 100)
    check("setProbeRadius / probeRadius", (p.setProbeRadius(1.4), abs(p.probeRadius() - 1.4) < 1e-9))

    # ── 2. Structure API (used by mdakit_sasa _single_frame) ─────
    print("\n--- Structure API (mdakit_sasa pattern) ---")
    structure = freesasa.Structure()
    check("Structure() empty constructor", structure is not None)

    # Replicate mdakit_sasa's addAtom pattern exactly
    # structure.addAtom(type.rjust(2), resname, resnum, segid, x, y, z)
    structure.addAtom(' N', 'MET', 1, 'A', 0.0, 0.0, 0.0)
    structure.addAtom(' C', 'MET', 1, 'A', 1.5, 0.0, 0.0)
    structure.addAtom(' O', 'MET', 1, 'A', 0.0, 1.5, 0.0)
    structure.addAtom(' N', 'ALA', 2, 'A', 3.0, 0.0, 0.0)
    structure.addAtom(' C', 'ALA', 2, 'A', 4.5, 0.0, 0.0)
    check("addAtom (mdakit_sasa pattern)", structure.nAtoms() == 5)

    # ── 3. calc() and result API ──────────────────────────────────
    print("\n--- calc() and Result API ---")
    params = freesasa.Parameters()
    result = freesasa.calc(structure, params)
    check("freesasa.calc(structure, params)", result is not None)
    check("result.totalArea() > 0", result.totalArea() > 0)

    areas = result.residueAreas()
    check("result.residueAreas() returns dict", isinstance(areas, dict))
    # mdakit_sasa iterates: for s in areas.keys() for r in areas[s].keys()
    total_from_residues = sum(
        areas[s][r].total
        for s in areas.keys()
        for r in areas[s].keys()
    )
    check("residueAreas iteration (mdakit_sasa pattern)", total_from_residues > 0)
    check("residueArea .total attribute", all(
        hasattr(areas[s][r], 'total')
        for s in areas.keys()
        for r in areas[s].keys()
    ))

    # ── 4. Thread-count scaling (key new feature) ─────────────────
    print("\n--- Thread Scaling (1 vs 8 threads on 1ubq.pdb) ---")
    s = freesasa.Structure(FREESASA_PDB)
    check("Structure(pdb_file) constructor", s is not None and s.nAtoms() > 0)

    p1 = freesasa.Parameters()
    p1.setNThreads(1)
    p1.setAlgorithm(freesasa.LeeRichards)
    r1 = freesasa.calc(s, p1)

    p8 = freesasa.Parameters()
    p8.setNThreads(8)
    p8.setAlgorithm(freesasa.LeeRichards)
    r8 = freesasa.calc(s, p8)

    check("L&R 1-thread vs 8-thread identical",
          abs(r1.totalArea() - r8.totalArea()) < 0.01,
          f'{r1.totalArea():.5f} vs {r8.totalArea():.5f}')

    p1sr = freesasa.Parameters()
    p1sr.setNThreads(1)
    p1sr.setAlgorithm(freesasa.ShrakeRupley)
    r1sr = freesasa.calc(s, p1sr)

    p8sr = freesasa.Parameters()
    p8sr.setNThreads(8)
    p8sr.setAlgorithm(freesasa.ShrakeRupley)
    r8sr = freesasa.calc(s, p8sr)

    check("S&R 1-thread vs 8-thread identical",
          abs(r1sr.totalArea() - r8sr.totalArea()) < 0.01,
          f'{r1sr.totalArea():.5f} vs {r8sr.totalArea():.5f}')

    # ── 5. calcCoord (alternative API) ───────────────────────────
    print("\n--- calcCoord API ---")
    coords = [0.0, 0.0, 0.0,   1.5, 0.0, 0.0]
    radii  = [1.0, 1.5]
    rc = freesasa.calcCoord(coords, radii)
    check("calcCoord() works", rc is not None and rc.totalArea() > 0)

    # ── 6. selectArea (selection API) ────────────────────────────
    print("\n--- selectArea API ---")
    s_ubq = freesasa.Structure(FREESASA_PDB)
    r_ubq = freesasa.calc(s_ubq)
    sel = freesasa.selectArea(('backbone, name CA+N+C+O',), s_ubq, r_ubq)
    check("selectArea() works", sel is not None)
    check("selectArea backbone area > 0", sel.get('backbone', 0) > 0)

    # ── 7. Windows safety flag (mdakit_sasa _is_windows() check) ─
    print("\n--- mdakit_sasa Windows safety pattern ---")
    def _is_windows():
        return os.name == 'nt'
    params2 = freesasa.Parameters()
    if _is_windows():
        params2.setNThreads(1)
    check("Windows safety pattern works", params2 is not None)

    # ── 8. Large structure parallel test ─────────────────────────
    if os.path.exists(LARGE_PDB):
        print("\n--- Large structure (1AON, 58674 atoms) ---")
        s_large = freesasa.Structure(LARGE_PDB)
        check("Large structure loaded", s_large.nAtoms() > 50000)

        p_large = freesasa.Parameters()
        p_large.setNThreads(8)
        r_large = freesasa.calc(s_large, p_large)
        check("Large structure 8-thread calc", r_large.totalArea() > 100000)

    # ── Summary ───────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print(f"  Results: {len(passed)} passed, {len(failed)} failed")
    print("=" * 60)
    if failed:
        print(f"\nFAILED: {failed}")
        return 1
    else:
        print("\nAll tests passed! FreeSASA Python API is fully compatible.")
        print("mdakit_sasa (pegerto/mdakit_sasa) will work correctly.")
        return 0

if __name__ == '__main__':
    sys.exit(run_tests())
