#!/usr/bin/env bash
# Run each chapter's Python script and its Octave twin and compare what they print.
#
#   tools/compare-routes.sh                 every chapter that has an Octave twin
#   tools/compare-routes.sh ch05 chf2       only the chapters named
#
# The two routes must print the same lines and both must exit 0. Before the
# comparison, trailing blanks are dropped and every run of spaces becomes one
# space: Python pads a label to a width counted in characters and Octave's
# printf counts bytes, so a label holding a non-ASCII letter such as a Greek
# epsilon gets fewer padding spaces from Octave. Nothing else is forgiven.
#
# PYTHON and OCTAVE name the two interpreters, python3 and octave by default:
#
#   PYTHON=.venv/bin/python tools/compare-routes.sh
#
# Exit status is 0 when every chapter passes, 1 when any differs, and 2 when
# an interpreter is missing.

set -u
cd "$(dirname "$0")/.." || exit 2

PYTHON="${PYTHON:-python3}"
OCTAVE="${OCTAVE:-octave}"
for tool in "$PYTHON" "$OCTAVE"; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "compare-routes: $tool is not on PATH" >&2
    exit 2
  fi
done

normalize() {
  sed -E -e 's/[[:space:]]+$//' -e 's/ +/ /g'
}

stems=()
if [ "$#" -gt 0 ]; then
  for want in "$@"; do
    for f in octave/"$want"_*.m octave/"$want".m; do
      [ -f "$f" ] && stems+=("$(basename "$f" .m)")
    done
  done
else
  for f in octave/ch*.m; do
    [ -f "$f" ] && stems+=("$(basename "$f" .m)")
  done
fi
if [ "${#stems[@]}" -eq 0 ]; then
  echo "compare-routes: no Octave twin matches" >&2
  exit 1
fi

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
failed=0
for stem in "${stems[@]}"; do
  py="python/$stem.py"
  if [ ! -f "$py" ]; then
    echo "FAIL  $stem  (no python/$stem.py to compare with)"
    failed=1
    continue
  fi
  start=$SECONDS
  "$PYTHON" "$py" > "$work/py.raw" 2>&1
  py_status=$?
  py_time=$((SECONDS - start))
  start=$SECONDS
  "$OCTAVE" --no-window-system --norc --quiet "octave/$stem.m" > "$work/oct.raw" 2>&1
  oct_status=$?
  oct_time=$((SECONDS - start))
  normalize < "$work/py.raw" > "$work/py.txt"
  normalize < "$work/oct.raw" > "$work/oct.txt"
  lines=$(wc -l < "$work/py.txt" | tr -d ' ')
  if [ "$py_status" -eq 0 ] && [ "$oct_status" -eq 0 ] && cmp -s "$work/py.txt" "$work/oct.txt"; then
    echo "PASS  $stem  ($lines lines; Python ${py_time}s, Octave ${oct_time}s)"
  else
    echo "FAIL  $stem  (Python exit $py_status, Octave exit $oct_status)"
    diff "$work/py.txt" "$work/oct.txt" | head -n 20 | sed 's/^/      /'
    failed=1
  fi
done
exit "$failed"
