#!/usr/bin/env bash
# Run an executable for both schemes (1=Upwind, 2=Centred) in scratch
# directories and compare every output against reference/{upwind,centred}.
#
# Usage: tools/check_against_reference.sh ./NS_lid_driven_cavity_cpp [compare options]
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
EXE=$(realpath "$1"); shift
status=0
for scheme in 1:upwind 2:centred; do
  Q=${scheme%%:*}; name=${scheme##*:}
  run="$ROOT/build/run_$name"
  rm -rf "$run"; mkdir -p "$run"
  ( cd "$run" && echo "$Q" | NS_DUMP_FIELDS=1 "$EXE" > stdout.txt 2> stderr.txt ) || {
    echo "== $name: executable FAILED (exit $?)"; tail -5 "$run/stderr.txt"; status=1; continue; }
  grep -E 'Nombre|non convergence' "$run/stdout.txt" | awk '{print $NF}' > "$run/iccg2_iterations.txt"
  echo "== $name ($run)"
  python3 "$ROOT/tools/compare_reference.py" "$run" "$ROOT/reference/$name" "$@" | grep -v '^ok '
  [ "${PIPESTATUS[0]}" -eq 0 ] || status=1
done
exit $status
