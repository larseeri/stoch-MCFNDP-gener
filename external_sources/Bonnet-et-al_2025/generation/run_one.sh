#!/bin/bash
# Generates one stochastic instance with stoch_gen.
# Usage: run_one.sh BASE_ID INSTANCE_ID
# Runs the configuration scenarios/BASE_ID/INSTANCE_ID.par from its folder (its paths are
# relative to it) and writes the output of stoch_gen to INSTANCE_ID.log. The executable is
# the one named by the variable STOCH_GEN, or stoch_gen on the PATH. Skips instances whose
# scenario file is already complete.
GEN="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$GEN")"
BIN="${STOCH_GEN:-$(command -v stoch_gen)}"
if [ -z "$BIN" ] || [ ! -x "$BIN" ]; then
    echo "run_one.sh: stoch_gen not found; put it on the PATH or set STOCH_GEN to its path" >&2
    exit 1
fi
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
cd "$ROOT/scenarios/$1" || exit 1
[ -f "$2.par" ] || { echo "run_one.sh: no configuration $1/$2.par" >&2; exit 1; }
if grep -q "Saved" "$2.log" 2>/dev/null && [ -s "$2.txt" ]; then exit 0; fi
rm -f "$2.txt" "$2.log" "$2_moms.dat" "$2_corrs.dat" "$2_raw.dat"
start=$(date +%s)
nice -n 10 "$BIN" +F "$2.par" > "$2.log" 2>&1
status=$?
echo "$(date '+%F %T') $2 exit=$status seconds=$(( $(date +%s) - start )) $(grep -h 'Saved' "$2.log")" >> "$GEN/progress.log"
