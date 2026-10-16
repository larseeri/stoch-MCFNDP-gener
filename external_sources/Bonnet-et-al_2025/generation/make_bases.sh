#!/bin/bash
# Regenerates the 20 base instances with the code of Bonnet et al. and checks them.
#
# Usage: generation/make_bases.sh BONNET_REPO [PYTHON]
#   BONNET_REPO  clone of https://github.com/LBonnet159/SSNDP-Realistic-Instance-Generator
#                (commit e390d35 was used; it provides the code and the benchmark networks)
#   PYTHON       interpreter with numpy, scipy, networkx and matplotlib (default: python3)
#
# For each row of design.csv, runs src/main.py of BONNET_REPO in a temporary copy, on the
# network of the row alone, with the settings of the repository's src/Config.txt except
#   doStatic=True networkNb=0 instanceNb=1 horizon=10000 discretization=1000
#   commodityNb=<commodities of the row> demandSeed=<demand_seed of the row>
# (the large horizon only avoids an abort of a travel-time check that the code also
# performs for static instances; it does not change SNDP files). The instance is then
# converted to DOW with bonnet2std.sh. If base/<base>.txt and .dow exist, the new files
# are compared with them; otherwise they are written there. The code of Bonnet et al. uses
# annotations of Python 3.10; with Python 3.9, the copy gets 'from __future__ import
# annotations', which does not change its results (the set was built this way, with
# Python 3.9.6, numpy 2.0.2, scipy 1.13.1, networkx 3.2.1 and matplotlib 3.9.4).
set -u
GEN="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$GEN")"
[ $# -ge 1 ] || { sed -n '4,8p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }
REPO="$(cd "$1" && pwd)" || exit 1
PY="${2:-python3}"
NETS="$REPO/data/Networks/Benchmark"
[ -f "$REPO/src/main.py" ] && [ -d "$NETS" ] || { echo "make_bases.sh: $1 is not a clone of SSNDP-Realistic-Instance-Generator" >&2; exit 1; }
git -C "$REPO" log -1 --format='Bonnet code: commit %h (%ad)' 2>/dev/null
"$PY" -c 'import sys, numpy, scipy, networkx, matplotlib; print("Python", sys.version.split()[0], "numpy", numpy.__version__, "scipy", scipy.__version__, "networkx", networkx.__version__)' || exit 1

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cp -R "$REPO/src" "$WORK/"
if ! "$PY" -c 'import sys; sys.exit(sys.version_info < (3, 10))'; then
    for f in "$WORK"/src/*.py; do
        { echo 'from __future__ import annotations'; cat "$f"; } > "$f.tmp" && mv "$f.tmp" "$f"
    done
fi
mkdir -p "$ROOT/base"
bad=0
while IFS=, read -r base row _ _ _ _ _ _ k seed network instance; do
    instance="${instance%$'\r'}"
    [ "$base" = base_id ] && continue
    folder=$(printf 'Row%02d' "$row")
    mkdir -p "$WORK/data/Networks/$folder"
    cp "$NETS/$network" "$WORK/data/Networks/$folder/"
    sed -e "s/^folder=.*/folder=$folder/" -e 's/^networkNb=.*/networkNb=0/' \
        -e 's/^doStatic=.*/doStatic=True/' -e 's/^instanceNb=.*/instanceNb=1/' \
        -e 's/^horizon=.*/horizon=10000/' -e 's/^discretization=.*/discretization=1000/' \
        -e "s/^commodityNb=.*/commodityNb=$k/" -e "s/^demandSeed=.*/demandSeed=$seed/" \
        "$REPO/src/Config.txt" > "$WORK/src/Config.txt"
    out="$WORK/data/Instances/$folder/$instance"
    if ! (cd "$WORK" && "$PY" src/main.py > "$WORK/$folder.log" 2>&1) || [ ! -s "$out" ]; then
        echo "$base: generation failed (expected $instance; found: $(ls "$WORK/data/Instances/$folder" 2>/dev/null))"; tail -3 "$WORK/$folder.log"; bad=1; continue
    fi
    "$ROOT/bonnet2std.sh" -q -F dow -o "$WORK/$base.dow" "$out" || { echo "$base: conversion failed"; bad=1; continue; }
    for ext in txt dow; do
        new="$out"; [ $ext = dow ] && new="$WORK/$base.dow"
        ref="$ROOT/base/$base.$ext"
        if [ ! -e "$ref" ]; then
            cp "$new" "$ref"; echo "$base.$ext: written"
        elif cmp -s "$new" "$ref"; then
            echo "$base.$ext: identical"
        else
            echo "$base.$ext: DIFFERENT from base/$base.$ext"; bad=1
        fi
    done
done < "$GEN/design.csv"
[ $bad = 0 ] && echo "all base instances reproduced" || { echo "some base instances differ or failed"; exit 1; }
