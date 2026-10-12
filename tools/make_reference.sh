#!/usr/bin/env bash
# Build the original Fortran solver and (re)generate the reference oracle used
# to validate the C++ port:  reference/upwind/  and  reference/centred/
#
# Usage:  tools/make_reference.sh [output-dir]     (default: <repo>/reference)
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DEST=${1:-$ROOT/reference}
FORTRAN_DIR="$ROOT/fortran"
[ -d "$FORTRAN_DIR" ] || FORTRAN_DIR="$ROOT"   # Fortran sources still in the root during the port

make -C "$FORTRAN_DIR" -s
EXE="$FORTRAN_DIR/NS_lid_driven_cavity"

for scheme in 1:upwind 2:centred; do
  Q=${scheme%%:*}; name=${scheme##*:}
  out="$DEST/$name"
  rm -rf "$out"; mkdir -p "$out"
  ( cd "$out" && echo "$Q" | NS_DUMP_FIELDS=1 "$EXE" > stdout.txt )
  # keep the oracle small: only a few snapshots (first steps, some middle, last)
  find "$out" -regex '.*/\(velocity\|omega\|press\|fields_\)[0-9]+\.\(vec\|scl\|txt\)' \
       ! -name '*00000[123].*' ! -name '*000200.*' ! -name '*000400.*' ! -name '*001000.*' \
       ! -name '*010000.*' ! -name '*020000.*' -delete
  # stdout: keep only the ICCG2 iteration counts (the 'time =' lines are noise)
  grep -E 'Nombre|non convergence' "$out/stdout.txt" | awk '{print $NF}' > "$out/iccg2_iterations.txt"
  rm "$out/stdout.txt"
  echo "$out: $(ls "$out" | wc -l) files"
done
