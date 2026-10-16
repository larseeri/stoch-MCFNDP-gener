#!/usr/bin/env bash
#
# bonnet2std.sh -- convert SNDP instance files produced by the SSNDP realistic
# instance generator of Bonnet, Belieres, Hewitt & Ngueveu
# (https://github.com/LBonnet159/SSNDP-Realistic-Instance-Generator)
# into complete MCFNDP instances in STD (a.k.a. Canad) format, or optionally
# DOW format, as read by stoch_gen
# (https://github.com/larseeri/stoch-MCFNDP-gener). Both formats are described
# in STD-format-desc.txt and DOW-format-desc.txt of the deterministic generator
# (https://github.com/larseeri/determ-MCFNDP-gener).
#
# Input (Bonnet format):
#
#   NODES,<n>
#   node id, cluster id[, x, y]
#   ARCS,<m>
#   arc id, origin id, destination id, unit cost, fixed cost, capacity[, time]
#   COMMODITIES,<k>
#   commodity id, origin id, destination id, quantity[, available, due]
#   [later lines (horizon=..., time windows, ...) are ignored]
#
# Mapping: node j -> node j-low+1; arc a -> arc a+1 (same order); commodity
# h -> commodity h+1 (same order). Unit costs are not commodity dependent, and
# each commodity has a single origin and a single destination. SSNDP instance
# files are accepted, but their time data are ignored (static counterpart).
#
# STD output. By default every value is rounded to the nearest integer, after
# optional scaling, because versions of stoch_gen up to commit af1d288
# (January 2026) read STD files into C++ 'int' variables (InstanceRStd::read()
# and InstanceGStd::read() in GENSTOCH/instance.cpp), as do other codes for the
# Canad instances:
#
#   volume      = round(S * quantity)
#   capacity    = round(S * capacity)
#   fixed cost  = round(C * fixed cost)
#   unit cost   = round(C * unit cost / S)
#
# where S (-s) rescales the flow unit and C (-c) rescales the cost unit; the
# objective of the converted instance is then C times the original one. Since
# commit 280ee8c (October 2026), stoch_gen reads STD values into 'double'
# variables; with -x, the values are then copied unchanged from the input file,
# without scaling or rounding. The commodity-specific capacity of every commodity on every arc is the arc
# capacity (i.e. no additional restriction) or, with -t, min(volume,
# capacity). The same readers ignore the commodity number written on each
# commodity line and fill the values by position, so every arc lists all
# commodities in increasing order, all with a positive capacity.
#
# DOW output (-F dow). stoch_gen reads DOW files into 'double' variables, so
# the values are copied unchanged from the input file (no rounding).

set -euo pipefail

prog=$(basename "$0")

usage() {
    cat <<EOF
Usage: $prog [options] [-o OUTPUT] INPUT
       $prog [options] -d OUTDIR INPUT...

Convert Bonnet et al. SNDP/SSNDP instance files to complete MCFNDP instances
in STD (default) or DOW format, readable by stoch_gen.

  -o OUTPUT   output file for a single INPUT ('-' for standard output);
              default: INPUT with its extension replaced by .std or .dow
  -d OUTDIR   write <basename>.std (or .dow) for each INPUT into OUTDIR
  -F FORMAT   std (default) or dow
  -s FACTOR   STD only: flow scale; quantities and capacities are multiplied
              and unit costs divided by FACTOR before rounding (default 1)
  -c FACTOR   STD only: cost scale; fixed and unit costs are multiplied by
              FACTOR before rounding (default 1)
  -t          STD only: commodity-specific capacities min(volume, capacity)
              instead of the arc capacity
  -x          STD only: exact values, copied unchanged from INPUT, instead of
              integers (for stoch_gen 280ee8c or later); excludes -s and -c
  -q          quiet: do not print the per-file summary
  -h          show this help

Without -x, the summary reports, for each kind of value, the largest
relative rounding error and how many non-zero values were rounded to zero.
EOF
}

out=""
outdir=""
fmt=std
fscale=1
cscale=1
tight=0
exact=0
scaled=0
quiet=0
stdOnly=""
while getopts ":o:d:F:s:c:txqh" opt; do
    case $opt in
        o) out=$OPTARG ;;
        d) outdir=$OPTARG ;;
        F) fmt=$(printf '%s' "$OPTARG" | tr '[:upper:]' '[:lower:]') ;;
        s) fscale=$OPTARG; scaled=1; stdOnly="$stdOnly -s" ;;
        c) cscale=$OPTARG; scaled=1; stdOnly="$stdOnly -c" ;;
        t) tight=1; stdOnly="$stdOnly -t" ;;
        x) exact=1; stdOnly="$stdOnly -x" ;;
        q) quiet=1 ;;
        h) usage; exit 0 ;;
        :) echo "$prog: option -$OPTARG requires an argument" >&2; usage >&2; exit 2 ;;
        *) echo "$prog: unknown option -$OPTARG" >&2; usage >&2; exit 2 ;;
    esac
done
shift $((OPTIND - 1))

if [ $# -eq 0 ]; then
    usage >&2
    exit 2
fi
case $fmt in
    std|dow) ;;
    *) echo "$prog: -F must be std or dow" >&2; exit 2 ;;
esac
if [ "$fmt" = dow ] && [ -n "$stdOnly" ]; then
    echo "$prog:$stdOnly: only available with STD output (DOW keeps the original values)" >&2
    exit 2
fi
if [ $exact -eq 1 ] && [ $scaled -eq 1 ]; then
    echo "$prog: -x keeps the values unchanged and cannot be combined with -s or -c" >&2
    exit 2
fi
for f in "$fscale" "$cscale"; do
    if ! printf '%s' "$f" | grep -Eq '^([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$' ||
       ! awk -v x="$f" 'BEGIN { exit !(x + 0 > 0) }'; then
        echo "$prog: scale factor '$f' must be a positive number" >&2
        exit 2
    fi
done
if [ -n "$out" ] && [ -n "$outdir" ]; then
    echo "$prog: -o and -d cannot be used together" >&2
    exit 2
fi
if [ -n "$out" ] && [ $# -gt 1 ]; then
    echo "$prog: -o accepts a single INPUT; use -d for several" >&2
    exit 2
fi
if [ -n "$outdir" ]; then
    mkdir -p "$outdir"
fi

# Converts one file. Arguments: input path, output path ('-' for stdout).
convert() {
    local in=$1 dest=$2 tmp summary
    tmp=$(mktemp "${TMPDIR:-/tmp}/bonnet2std.XXXXXX")

    if ! summary=$(awk -F',' -v src="$(basename "$in")" -v out="$tmp" -v fmt="$fmt" \
                       -v S="$fscale" -v C="$cscale" -v tight="$tight" -v exact="$exact" '
        function fail(msg) { printf "%s: %s\n", src, msg > "/dev/stderr"; err = 1; exit 1 }
        function isint(s)  { return s ~ /^-?[0-9]+$/ }
        function isnum(s)  { return s ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/ }
        function trim(s)   { gsub(/^[ \t]+|[ \t]+$/, "", s); return s }
        # Non-negative number in field i of the current line (as written,
        # without a leading plus sign).
        function num(i, what,    s) {
            s = trim($i)
            if (!isnum(s) || s + 0 < 0) fail("line " NR ": invalid " what " \"" s "\"")
            sub(/^\+/, "", s)
            return s
        }
        # Node id in field i of the current line, checked against NODES.
        function nodeId(i, what,    s) {
            s = trim($i)
            if (!isint(s)) fail("line " NR ": invalid " what " \"" s "\"")
            s += 0
            if (!(s in node)) fail("line " NR ": " what " " s " is not a declared node")
            return s
        }
        # Rounds x >= 0 to an integer that stoch_gen can read into an int,
        # and records the rounding error for the kind of value t.
        function rnd(x, t,    r, e) {
            r = int(x + 0.5)
            if (r > 2147483647)
                fail("scaled " t " " x " exceeds the range of a C int; use a smaller -c or -s, or -x")
            if (x > 0) {
                e = (r > x ? r - x : x - r) / x
                if (e > maxErr[t]) maxErr[t] = e
                if (r == 0) zeroed[t]++
            }
            return r
        }

        { sub(/\r$/, "") }                       # tolerate CRLF line endings
        /^[ \t]*$/ { next }

        # Section headers start with a letter; parsing stops at the first line
        # starting with a letter after the COMMODITIES section.
        /^[A-Za-z]/ {
            key = toupper(trim($1))
            if (sawC) exit
            if (key == "NODES")       { sec = "N"; declN = trim($2); sawN = 1; next }
            if (key == "ARCS")        { sec = "A"; declA = trim($2); sawA = 1; next }
            if (key == "COMMODITIES") { sec = "C"; declK = trim($2); sawC = 1; next }
            fail("unexpected line " NR " before COMMODITIES section: " $0)
        }

        sec == "N" {
            id = trim($1)
            if (!isint(id))   fail("line " NR ": invalid node id \"" id "\"")
            id += 0
            if (id in node)   fail("line " NR ": duplicate node id " id)
            node[id] = 1; nN++
            if (nN == 1 || id < low)  low = id
            if (nN == 1 || id > high) high = id
            next
        }

        sec == "A" {
            if (NF < 6) fail("line " NR ": arc line has fewer than 6 fields")
            nRead++
            o = nodeId(2, "arc origin"); d = nodeId(3, "arc destination")
            u = num(4, "unit cost"); f = num(5, "fixed cost"); c = num(6, "capacity")
            if (o == d) { selfLoops++; next }    # a loop cannot carry flow
            if ((o SUBSEP d) in seen) parallel++
            seen[o, d] = 1
            nA++
            aO[nA] = o; aD[nA] = d; aU[nA] = u; aF[nA] = f; aC[nA] = c
            next
        }

        sec == "C" {
            if (NF < 4) fail("line " NR ": commodity line has fewer than 4 fields")
            if (NF >= 6) timed = 1
            o = nodeId(2, "commodity origin"); d = nodeId(3, "commodity destination")
            q = num(4, "quantity")
            if (o == d) fail("line " NR ": commodity origin and destination are both node " o)
            if (q + 0 == 0) fail("line " NR ": commodity quantity is zero")
            nK++
            kO[nK] = o; kD[nK] = d; kQ[nK] = q
            next
        }

        { fail("line " NR ": data before NODES section") }

        END {
            if (err) exit 1
            if (!sawN)  fail("no NODES section found")
            if (!sawA)  fail("no ARCS section found")
            if (!sawC)  fail("no COMMODITIES section found (a network file has no commodities; generate an SNDP instance with the Bonnet generator)")
            if (nN == 0) fail("NODES section is empty")
            if (nA == 0) fail("ARCS section contains no usable arc")
            if (nK == 0) fail("COMMODITIES section is empty")
            if (isint(declN) && declN + 0 != nN)    fail("NODES header announces " declN " nodes but " nN " were read")
            if (isint(declA) && declA + 0 != nRead) fail("ARCS header announces " declA " arcs but " nRead " were read")
            if (isint(declK) && declK + 0 != nK)    fail("COMMODITIES header announces " declK " commodities but " nK " were read")
            if (high - low + 1 != nN)
                fail("node ids are not contiguous (" nN " ids in [" low ", " high "])")
            base = low - 1                       # output nodes are numbered 1..n

            if (fmt == "dow") {
                # Values are copied verbatim; stoch_gen reads them as doubles.
                print " MULTIGEN.DAT:" > out
                printf "%8d%12d%12d\n", nN, nA, nK > out
                for (a = 1; a <= nA; a++)
                    printf "%8d %11d %15s %15s %19s %3d %7d\n", aO[a] - base, aD[a] - base, \
                           aU[a], aC[a], aF[a], 1, a > out
                for (h = 1; h <= nK; h++)
                    printf "%8d %11d %15s\n", kO[h] - base, kD[h] - base, kQ[h] > out
            } else if (exact) {
                # Values are copied verbatim; stoch_gen 280ee8c and later read
                # them as doubles.
                printf "%8d%12d%12d\n", nN, nA, nK > out
                for (a = 1; a <= nA; a++) {
                    if (aC[a] + 0 == 0) fail("arc " a - 1 " has capacity 0; commodity capacities must be positive")
                    printf "%8d %11d %19s %15s %11d\n", aO[a] - base, aD[a] - base, aF[a], aC[a], nK > out
                    for (h = 1; h <= nK; h++)
                        printf "%8d %15s %15s\n", h, aU[a], \
                               ((tight && kQ[h] + 0 < aC[a] + 0) ? kQ[h] : aC[a]) > out
                }
                for (h = 1; h <= nK; h++) {
                    o = kO[h] - base; d = kD[h] - base
                    if (o < d) {
                        printf "%8d %11d %15s\n", h, o, kQ[h] > out
                        printf "%8d %11d %15s\n", h, d, "-" kQ[h] > out
                    } else {
                        printf "%8d %11d %15s\n", h, d, "-" kQ[h] > out
                        printf "%8d %11d %15s\n", h, o, kQ[h] > out
                    }
                }
            } else {
                for (h = 1; h <= nK; h++) {
                    vol[h] = rnd(S * kQ[h], "volumes")
                    if (vol[h] == 0) fail("commodity " h - 1 " has volume " kQ[h] ", which rounds to 0; use a larger -s")
                }
                printf "%8d%12d%12d\n", nN, nA, nK > out
                for (a = 1; a <= nA; a++) {
                    cap = rnd(S * aC[a], "capacities")
                    if (cap == 0) fail("arc " a - 1 " has capacity " aC[a] ", which rounds to 0; use a larger -s")
                    printf "%8d%12d%12d%12d%12d\n", aO[a] - base, aD[a] - base, \
                           rnd(C * aF[a], "fixed costs"), cap, nK > out
                    var = rnd(C * aU[a] / S, "unit costs")
                    for (h = 1; h <= nK; h++)
                        printf "%8d%12d%12d\n", h, var, ((tight && vol[h] < cap) ? vol[h] : cap) > out
                }
                # Supply (> 0) and demand (< 0) of each commodity, by node.
                for (h = 1; h <= nK; h++) {
                    o = kO[h] - base; d = kD[h] - base
                    if (o < d) {
                        printf "%8d%12d%12d\n", h, o,  vol[h] > out
                        printf "%8d%12d%12d\n", h, d, -vol[h] > out
                    } else {
                        printf "%8d%12d%12d\n", h, d, -vol[h] > out
                        printf "%8d%12d%12d\n", h, o,  vol[h] > out
                    }
                }
            }
            close(out)

            msg = toupper(fmt) ": " nN " nodes, " nA " arcs, " nK " commodities"
            if (selfLoops) msg = msg "; dropped " selfLoops " self-loop(s)"
            if (parallel)  msg = msg "; " parallel " parallel arc(s)"
            if (timed)     msg = msg "; SSNDP time data ignored"
            if (fmt == "std" && exact) msg = msg "; exact values (no rounding)"
            if (fmt == "std" && !exact) {
                if (S != 1 || C != 1) msg = msg "; scales -s " S " -c " C
                msg = msg "\n    max relative rounding error:"
                n = split("volumes,capacities,fixed costs,unit costs", kinds, ",")
                for (i = 1; i <= n; i++) {
                    t = kinds[i]
                    msg = msg sprintf(" %s %.3g%%", t, 100 * maxErr[t])
                    if (zeroed[t]) msg = msg " (" zeroed[t] " rounded to 0)"
                    msg = msg (i < n ? "," : "")
                }
                if (zeroed["fixed costs"] || zeroed["unit costs"])
                    msg = msg "\n    WARNING: some non-zero costs were rounded to 0; consider a larger -c"
            }
            print msg
        }
    ' "$in"); then
        rm -f "$tmp"
        return 1
    fi

    if [ "$dest" = "-" ]; then
        cat "$tmp"
        rm -f "$tmp"
    elif ! mv "$tmp" "$dest"; then
        rm -f "$tmp"
        return 1
    fi
    if [ $quiet -eq 0 ]; then
        echo "$in -> $dest: $summary" >&2
    fi
}

status=0
for in in "$@"; do
    if [ ! -f "$in" ]; then
        echo "$prog: $in: no such file" >&2
        status=1
        continue
    fi
    base=$(basename "$in")
    stem=${base%.*}
    if [ -n "$out" ]; then
        dest=$out
    elif [ -n "$outdir" ]; then
        dest=$outdir/$stem.$fmt
    else
        dest=$(dirname "$in")/$stem.$fmt
    fi
    convert "$in" "$dest" || status=1
done
exit $status
