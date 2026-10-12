#!/usr/bin/env python3
"""Compare a solver run directory against a reference oracle directory.

    tools/compare_reference.py <run_dir> <reference_dir> [--rtol 1e-12] [--atol 1e-15]

Every file present in <reference_dir> must exist in <run_dir> and match:
  * fields_*.txt            full-precision field dumps, compared value by value
                            (|a-b| <= atol + rtol*|b|); header lines must be identical
  * *.geo *.vec *.scl       Ensight Gold binaries, compared byte-for-byte
  * everything else         compared as exact text (.case, coordinates.txt,
                            u_velocity.txt, v_velocity.txt, omega_boundary.txt,
                            iccg2_iterations.txt)
Exit status 0 = all files match, 1 = mismatch, 2 = usage error.
"""
import argparse
import os
import sys


def parse_fields(path):
    header, values = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if len(parts) == 1:
                values.append(float(parts[0]))
            else:
                header.append((len(values), line.rstrip("\n")))
    return header, values


def compare_fields(run, ref, rtol, atol):
    hr, vr = parse_fields(run)
    hf, vf = parse_fields(ref)
    if len(vr) != len(vf):
        return "value count %d != %d" % (len(vr), len(vf))
    # header lines carry the scalar time/dt values: compare those numerically too
    if [h[0] for h in hr] != [h[0] for h in hf]:
        return "array layout differs"
    for (_, a), (_, b) in zip(hr, hf):
        ta, tb = a.split(), b.split()
        if ta[0] != tb[0] or len(ta) != len(tb):
            return "header %r != %r" % (a, b)
        for x, y in zip(ta[1:], tb[1:]):
            if x == y:
                continue
            try:
                fx, fy = float(x), float(y)
            except ValueError:
                return "header %r != %r" % (a, b)
            if abs(fx - fy) > atol + rtol * abs(fy):
                return "header %r != %r" % (a, b)
    worst, worst_i = 0.0, -1
    for i, (a, b) in enumerate(zip(vr, vf)):
        err = abs(a - b) - rtol * abs(b)
        if err > worst:
            worst, worst_i = err, i
    if worst > atol:
        name = [h[1].split()[0] for h in hf if h[0] <= worst_i][-1]
        return "max deviation beyond tolerance in %s: run=%r ref=%r" % (name, vr[worst_i], vf[worst_i])
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("reference_dir")
    ap.add_argument("--rtol", type=float, default=1e-12)
    ap.add_argument("--atol", type=float, default=1e-15)
    args = ap.parse_args()
    if not os.path.isdir(args.reference_dir):
        print("reference dir not found:", args.reference_dir)
        return 2

    failures = 0
    for name in sorted(os.listdir(args.reference_dir)):
        ref = os.path.join(args.reference_dir, name)
        run = os.path.join(args.run_dir, name)
        if not os.path.isfile(run):
            print("MISSING  %s" % name)
            failures += 1
            continue
        ext = os.path.splitext(name)[1]
        if name.startswith("fields_"):
            err = compare_fields(run, ref, args.rtol, args.atol)
        elif ext in (".geo", ".vec", ".scl"):
            with open(run, "rb") as a, open(ref, "rb") as b:
                da, db = a.read(), b.read()
            if len(da) != len(db):
                err = "size %d != %d bytes" % (len(da), len(db))
            elif da != db:
                off = next(i for i, (x, y) in enumerate(zip(da, db)) if x != y)
                err = "first differing byte at offset %d" % off
            else:
                err = None
        else:
            with open(run, errors="replace") as a, open(ref, errors="replace") as b:
                err = None if a.read() == b.read() else "text differs"
        if err:
            print("FAIL     %s: %s" % (name, err))
            failures += 1
        else:
            print("ok       %s" % name)

    if failures:
        print("\n%d file(s) differ from the reference" % failures)
        return 1
    print("\nall files match the reference")
    return 0


if __name__ == "__main__":
    sys.exit(main())
