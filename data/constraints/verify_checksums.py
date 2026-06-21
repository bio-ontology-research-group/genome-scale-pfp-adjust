#!/usr/bin/env python3
"""Verify the vendored constraint files against checksums.sha256.

Run from anywhere; paths are resolved relative to this script's directory. Exits
0 if every listed file matches its recorded SHA-256, non-zero otherwise. CI runs
this so an accidental edit to a pinned constraint file fails the build. See
PROVENANCE.md for what each file is and how to re-pin after an intended update.
"""
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, "checksums.sha256")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    if not os.path.exists(MANIFEST):
        print(f"ERROR: missing manifest {MANIFEST}", file=sys.stderr)
        return 1
    ok = True
    n = 0
    with open(MANIFEST) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            expected, name = line.split(None, 1)
            name = name.lstrip("*").strip()  # tolerate `sha256sum -b` markers
            path = os.path.join(HERE, name)
            if not os.path.exists(path):
                print(f"MISSING  {name}")
                ok = False
                continue
            actual = sha256(path)
            if actual != expected:
                print(f"MISMATCH {name}\n  expected {expected}\n  actual   {actual}")
                ok = False
            else:
                print(f"OK       {name}")
            n += 1
    if ok:
        print(f"\nAll {n} constraint files match checksums.sha256.")
        return 0
    print("\nChecksum verification FAILED.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
