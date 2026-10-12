#!/usr/bin/env bash
# Check the C++ port against the Fortran outside the default configuration.
# For each variant below, copies fortran/ and src/ to build/variants/<name>/,
# changes one parameter in both copies, builds both programs, runs both schemes
# and compares every output of the C++ run with the Fortran run
# (tools/compare_reference.py, zero tolerance by default).
#
# Usage: tools/check_variants.sh [compare options, e.g. --rtol 1e-12 --atol 1e-15]
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
CXX=${CXX:-g++}
CXXFLAGS=${CXXFLAGS:--std=c++17 -O2 -Wall -Wextra}
COMPARE_OPTS=("$@")
[ ${#COMPARE_OPTS[@]} -gt 0 ] || COMPARE_OPTS=(--rtol 0 --atol 0)

# name | sed on fortran/variables_module.f90 | sed on src/variables.hpp
VARIANTS=(
  "re1000|s/nu = 1e-02/nu = 1e-03/|s/nu = 1e-02/nu = 1e-03/"
  "grid21|s/nx=11,ny=11/nx=21,ny=21/|s/nx = 11, ny = 11/nx = 21, ny = 21/"
)

# apply a sed expression and fail if it did not change the file
edit() {
  cp "$2" "$2.orig"
  sed -i "$1" "$2"
  if cmp -s "$2" "$2.orig"; then echo "error: '$1' did not change $2"; return 1; fi
  rm "$2.orig"
}

status=0
for v in "${VARIANTS[@]}"; do
  IFS='|' read -r name fsed csed <<< "$v"
  dir="$ROOT/build/variants/$name"
  rm -rf "$dir"; mkdir -p "$dir"
  cp -r "$ROOT/fortran" "$dir/fortran"; cp -r "$ROOT/src" "$dir/src"
  rm -f "$dir"/fortran/*.o "$dir"/fortran/*.mod "$dir/fortran/NS_lid_driven_cavity"
  echo "== variant $name"
  edit "$fsed" "$dir/fortran/variables_module.f90" && edit "$csed" "$dir/src/variables.hpp" &&
  make -C "$dir/fortran" -s > "$dir/build.log" 2>&1 &&
  $CXX $CXXFLAGS "$dir"/src/*.cpp -o "$dir/NS_lid_driven_cavity" >> "$dir/build.log" 2>&1 ||
    { echo "   build FAILED (see $dir/build.log)"; status=1; continue; }

  for scheme in 1:upwind 2:centred; do
    Q=${scheme%%:*}; sname=${scheme##*:}
    for impl in fortran:fortran/NS_lid_driven_cavity cpp:NS_lid_driven_cavity; do
      run="$dir/run_${sname}_${impl%%:*}"; mkdir -p "$run"
      ( cd "$run" && echo "$Q" | NS_DUMP_FIELDS=1 "$dir/${impl#*:}" > "$run.stdout" 2> "$run.stderr" ) ||
        { echo "   $sname ${impl%%:*}: executable FAILED"; status=1; continue 2; }
      grep -E 'Nombre|non convergence' "$run.stdout" | awk '{print $NF}' > "$run/iccg2_iterations.txt"
    done
    echo "-- $sname: C++ vs Fortran"
    python3 "$ROOT/tools/compare_reference.py" "$dir/run_${sname}_cpp" "$dir/run_${sname}_fortran" \
      "${COMPARE_OPTS[@]}" | grep -v '^ok '
    [ "${PIPESTATUS[0]}" -eq 0 ] || status=1
  done
done
exit $status
