# Stochastic MCFNDP instances on the benchmark networks of Bonnet et al.

This folder documents a set of 20 instances of the two-stage stochastic multi-commodity
capacitated fixed-charge network design problem (MCFNDP) with random demands. Each
instance has two parts:
- a deterministic **base instance** built on one of the 180 realistic hub-and-spoke
  benchmark networks of Bonnet, Belieres, Hewitt and Ngueveu (2025);
- a set of demand **scenarios** produced by the stochastic generator `stoch_gen` of
  Larsen, Bisaillon, Cordeau and Frejinger.

The set shows that the stochastic generator applies to any deterministic base instance,
not only to those of the companion deterministic generator. **Every scenario of every
instance is feasible.** Every generation run is recorded with its seed and stream, so
each instance can be reproduced exactly.

**The scenario files themselves are not included** (515 MB in total). They are rebuilt
from the configurations of this folder (see *Reproducing*), and `checksums.sha256` gives
their SHA-256 checksums so that the rebuilt files can be verified. This folder (about 4 MB)
contains everything else:
- the base instances;
- the conversion script;
- the configuration, log, target moments and correlations of every generation run;
- the design, the scripts and the log of all generation experiments.

## Contents

| Path | Content |
|---|---|
| `base/<base>.txt` | Base instance in Bonnet's SNDP format (nodes, arcs, commodities) |
| `base/<base>.dow` | The same instance in DOW format, converted with `bonnet2std.sh -F dow` (values copied unchanged) |
| `bonnet2std.sh` | Converter from Bonnet's SNDP format to the STD or DOW format (see *Conversion*) |
| `scenarios/<base>/<instance>.par` | `stoch_gen` configuration of the instance, i.e. the content of its command line, including the seed (`r`) and the stream (`R`) |
| `scenarios/<base>/<instance>.log` | `stoch_gen` output: moment-matching error, feasibility report and time stamps |
| `scenarios/<base>/<instance>_moms.dat`, `_corrs.dat` | Target moments and correlations of the random demands, as written by `stoch_gen` |
| `scenarios/<base>/rejected/` | Configuration and log of each rejected attempt (B03 only) |
| `manifest.csv` | One row per instance: base features, settings, experiment, seed, stream, times and realized statistics |
| `checksums.sha256` | SHA-256 checksums of the scenario files `<instance>.txt` and of the raw scenario files `<instance>_raw.dat` (not included) |
| `generation/design.csv` | Design of the base instances: Bonnet network, number of commodities, demand seed |
| `generation/assignment.csv` | Spread and correlation assigned to each base instance |
| `generation/experiments.csv` | **Log of all generation experiments**, accepted and rejected, with seed, stream and duration |
| `generation/slack.csv` | λ\* of each base instance (see *Base instances*) |
| `generation/make_bases.sh` | Regenerates the base instances with the code of Bonnet et al. and compares them with `base/` |
| `generation/generate_set.py`, `run_one.sh` | Generation of the scenarios: whole procedure, or one instance |
| `generation/slack.py`, `summarize.py` | Computation of λ\* and of `manifest.csv` |
| `LICENSE-SSNDP-Realistic-Instance-Generator.txt` | License of the code and networks of Bonnet et al., from which the base instances derive |

Names encode the design. In `B07_D20_A1_R90_I1_K200_V10_C50`:
- `B07` is the row of the design;
- `D20`, `A1`, `R90` and `I1` identify the Bonnet network `DR30_A1_UF5_R90_D20_N50_I1_S0`
  (density 20%, hub ratio 1%, reciprocity 90%, replicate 1);
- `K200` is the number of commodities;
- `V10` is a demand spread of ±10%;
- `C50` is a correlation of 0.5 between demands (`C00`: independent demands).

## Base instances

The benchmark networks of Bonnet et al. have 50 nodes and no commodities. They cross four
arc densities, three hub ratios and three levels of reciprocity, with five replicates of
each combination, and every arc has capacity 1000. The commodities were generated with the
code of Bonnet et al. itself:
- repository commit `e390d35`;
- static SNDP instances, with the default `quantityToCapaMean=0.1` and `quantityToCapaDev=0.5`;
- origin–destination pairs drawn among the pairs joined by a path;
- quantities drawn from a normal distribution with mean 100 and standard deviation 50,
  truncated to [10, 1000];
- a different demand seed for each instance.

The 20 networks and commodity counts follow a balanced design (`generation/design.csv`):
* density 5%, 20%, 50% and 80% (122, 490, 1225 and 1397–1960 arcs): five instances each,
  one per replicate;
* hub ratio 1%, 10% and 20%, and reciprocity 20%, 50% and 90%: six or seven instances
  each, every density covering all three levels of each;
* 50, 100 or 200 commodities: seven, six and seven instances. The 5% density has no
  200-commodity instance: in a trial, a 5%-density network with 200 commodities could not
  route the demands within its arc capacities (λ\* = 0.67, see below).

The quantities of the base instances range from about 10 to 260. The column λ\* of the
manifest is the largest factor by which all demands of the base instance can be multiplied
simultaneously and still be routed with all arcs open (`generation/slack.py`). It
measures how tight the capacities are, from 1.40 (B03) to 66.8 (B16). Since feasibility is
monotone in the demands, every scenario whose demands stay below λ\* times their base values
is feasible.

## Conversion

`bonnet2std.sh` converts SNDP or SSNDP instance files of Bonnet et al. to complete MCFNDP
instances in the STD format (default) or the DOW format, both read by `stoch_gen`:

```
bonnet2std.sh [options] [-o OUTPUT] INPUT
bonnet2std.sh [options] -d OUTDIR INPUT...
```

It renumbers the nodes from 1 and keeps the order of arcs and commodities. The two
formats differ as follows:
* **DOW** (`-F dow`, used for this set): values are copied unchanged.
* **STD**: every arc lists every commodity with the unit cost and the capacity of the
  arc, since Bonnet instances have no commodity-specific data. By default, values are
  rounded to integers, for codes that expect the integer data of the original Canad
  files. The script then reports, for each kind of value, the largest relative rounding
  error and the number of non-zero values rounded to zero. Options:
  - `-s` and `-c`: flow and cost scale factors applied before rounding;
  - `-t`: commodity capacities min(volume, capacity);
  - `-x`: exact values instead of integers, for `stoch_gen` `280ee8c` or later.

`bonnet2std.sh -h` describes all options. The script needs only `bash` and `awk`. It
checks its input and stops with a message on a malformed file, or on a network file
without commodities. The STD and DOW formats are described in `program/STD-format-desc.txt`
and `program/DOW-format-desc.txt` of https://github.com/larseeri/determ-MCFNDP-gener.

## Stochastic instances

* **Random parameters:** the commodity demands only (`stoch_gen` option `S 1`). Costs,
  capacities and origin–destination pairs are those of the base instance.
* **Distribution:** target moments of a uniform distribution on
  [(1 − s) d, (1 + s) d] around the base demand d (α = β = s), and a common correlation of
  0 or 0.5 between distinct demands. The design assigns each of the ten combinations of
  s ∈ {0.1, 0.2, 0.3, 0.4, 0.5} and correlation 0 or 0.5 to two base instances, so that
  each spread appears once in every density and the combinations are balanced across
  commodity counts (`generation/assignment.csv`). Because of the feasibility rule below,
  B03 received s = 0.4 instead of 0.5.
* **Scenarios:** twice as many as commodities (100, 200 or 400), equiprobable. Moment
  matching (Høyland–Kaut–Wallace) needs more scenarios than random parameters. Its
  tolerance is 0.001 on moments and correlations, with at most 20 trials of 50 iterations.
* **Feasibility rule:** `stoch_gen` checks every scenario with HiGHS, with all arcs open. An
  attempt is accepted only if *all* its scenarios are feasible; otherwise s is lowered by
  0.05 and a new attempt is made. Only B03 (λ\* = 1.40) needed retries: 194 of 200 feasible
  scenarios at s = 0.5, 199 at s = 0.45, all 200 at s = 0.4.
* **Seeds and streams:** each run of `stoch_gen` is a numbered *experiment*:
  - all experiments use the seed 50163844 (option `r`);
  - experiment e uses the stream 47 + 100 (e − 1) (option `R`);
  - experiments 1–20 are the first attempts, in the order of the design (streams 47 to 1947);
  - retries are numbered from 21 on, after all experiments of the previous round.

  Of the 22 experiments, 20 were accepted (`generation/experiments.csv`).
* **Accuracy:** moment matching converged in its first trial for every instance.
  - The mean ratio of scenario demand to base demand is 1.000.
  - Its standard deviation is within 0.6% of the target s/√3.
  - The average correlation between distinct demands is 0.4994–0.4998 for a target of
    0.5, and 0.000 for a target of 0.
  - Moment matching does not impose the bounds of the uniform distribution, so some
    demands fall slightly outside [(1 − s) d, (1 + s) d], for instance between 0.40 d and
    1.58 d for s = 0.5. All demands are positive.

### Generation times

The set was generated on the following setup:
- a MacBook Pro with an Apple M1 Pro processor (8 cores) and 16 GB of memory, running macOS 26.6;
- `stoch_gen` commit `280ee8c`, compiled with Apple clang 21.0.0 and linked statically with HiGHS 1.11.0;
- three runs in parallel at low priority.

Run times range from 2 seconds to 79 minutes (median about one minute). Most of the time
goes to the feasibility checks: moment matching takes at most 1 second per instance,
while the checks account for 89–99% of every run longer than one minute. Only four
instances, all with 200 commodities and at least 1225 arcs, took more than ten minutes.
The large runs use about 650 MB of memory each.

The time records are in two files:
- `generation/experiments.csv` gives the duration of every experiment;
- the manifest splits it into moment matching (`hkw_seconds`) and feasibility checks
  (`check_seconds`), from the time stamps of the logs.

### Format of a stochastic instance

`<instance>.txt` is the output of `stoch_gen` for a DOW base instance: the scenarios
one after another, each starting with a line `SCENARIO <i>` (from 0) followed by a
complete DOW instance:

```
 SCENARIO        0
 MULTIGEN.DAT:
          50            490            200       <- nodes, arcs, commodities
           1             31         39.586       1000.000     791727.588              1             76
          ...                                    <- origin, destination, unit cost, capacity,
                                                    fixed cost, 1, arc index (490 lines)
          25              6         35.727       <- origin, destination, demand (200 lines)
          ...
 SCENARIO        1
          ...
```
(excerpt of `scenarios/B07_D20_A1_R90_I1_K200/B07_D20_A1_R90_I1_K200_V10_C50.txt`)

Only the demands differ from one scenario to the next. Other conventions:
- nodes are numbered from 1 (node j of the Bonnet file is node j + 1);
- arcs are listed by origin node, and the arc index in the last column refers to the
  order of `base/<base>.dow`;
- commodities keep the order of the base instance;
- `stoch_gen` prints values with three decimals.

## Instances

The last three columns give the duration of the accepted experiment, in seconds, in
total and in feasibility checks, and the size of the scenario file once regenerated.

| Instance | Arcs | Commodities | Total demand | λ\* | Spread | Correlation | Experiment | Stream | Scenarios | Time (s) | Checks (s) | Size (MB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `B01_D5_A1_R20_I0_K100_V40_C00` | 122 | 100 | 10,131 | 1.53 | ±40% | 0.0 | 1 | 47 | 200 | 33 | 32 | 3.4 |
| `B02_D5_A10_R20_I1_K50_V10_C00` | 122 | 50 | 4,832 | 1.62 | ±10% | 0.0 | 2 | 147 | 100 | 5 | 5 | 1.5 |
| `B03_D5_A20_R20_I2_K100_V40_C50` | 122 | 100 | 10,045 | 1.40 | ±40% (design ±50%) | 0.5 | 22 | 2147 | 200 | 7 | 5 | 3.4 |
| `B04_D5_A20_R50_I3_K50_V20_C00` | 122 | 50 | 5,373 | 1.85 | ±20% | 0.0 | 4 | 347 | 100 | 2 | 1 | 1.5 |
| `B05_D5_A20_R90_I4_K50_V30_C50` | 122 | 50 | 5,112 | 1.42 | ±30% | 0.5 | 5 | 447 | 100 | 2 | 2 | 1.5 |
| `B06_D20_A1_R50_I0_K100_V50_C50` | 490 | 100 | 10,845 | 11.75 | ±50% | 0.5 | 6 | 547 | 200 | 54 | 51 | 11.0 |
| `B07_D20_A1_R90_I1_K200_V10_C50` | 490 | 200 | 20,344 | 5.61 | ±10% | 0.5 | 7 | 647 | 400 | 232 | 220 | 23.7 |
| `B08_D20_A10_R50_I2_K200_V20_C00` | 490 | 200 | 20,754 | 8.19 | ±20% | 0.0 | 8 | 747 | 400 | 406 | 394 | 23.7 |
| `B09_D20_A10_R90_I3_K100_V30_C50` | 490 | 100 | 10,061 | 7.91 | ±30% | 0.5 | 9 | 847 | 200 | 72 | 69 | 11.0 |
| `B10_D20_A20_R20_I4_K50_V40_C00` | 490 | 50 | 5,121 | 11.71 | ±40% | 0.0 | 10 | 947 | 100 | 14 | 13 | 5.3 |
| `B11_D50_A1_R20_I0_K50_V50_C00` | 1225 | 50 | 4,933 | 45.48 | ±50% | 0.0 | 11 | 1047 | 100 | 25 | 23 | 12.8 |
| `B12_D50_A10_R20_I1_K200_V30_C00` | 1225 | 200 | 21,160 | 19.98 | ±30% | 0.0 | 12 | 1147 | 400 | 461 | 420 | 53.9 |
| `B13_D50_A10_R50_I2_K200_V20_C50` | 1225 | 200 | 19,680 | 9.67 | ±20% | 0.5 | 13 | 1247 | 400 | 1753 | 1720 | 53.9 |
| `B14_D50_A20_R50_I3_K100_V10_C00` | 1225 | 100 | 10,214 | 24.09 | ±10% | 0.0 | 14 | 1347 | 200 | 89 | 81 | 26.1 |
| `B15_D50_A20_R90_I4_K200_V40_C50` | 1225 | 200 | 20,072 | 9.24 | ±40% | 0.5 | 15 | 1447 | 400 | 4738 | 4706 | 53.9 |
| `B16_D80_A1_R50_I0_K50_V20_C50` | 1636 | 50 | 4,639 | 66.81 | ±20% | 0.5 | 16 | 1547 | 100 | 34 | 31 | 17.1 |
| `B17_D80_A1_R90_I1_K200_V10_C50` | 1960 | 200 | 21,199 | 29.42 | ±10% | 0.5 | 17 | 1647 | 400 | 2753 | 2692 | 84.2 |
| `B18_D80_A10_R20_I2_K50_V40_C50` | 1397 | 50 | 5,127 | 49.29 | ±40% | 0.5 | 18 | 1747 | 100 | 23 | 20 | 14.6 |
| `B19_D80_A10_R90_I3_K100_V30_C00` | 1960 | 100 | 10,524 | 47.28 | ±30% | 0.0 | 19 | 1847 | 200 | 147 | 131 | 41.2 |
| `B20_D80_A20_R50_I4_K200_V50_C00` | 1643 | 200 | 22,111 | 30.89 | ±50% | 0.0 | 20 | 1947 | 400 | 661 | 609 | 71.2 |

`manifest.csv` gives, for each instance:
- the base features: density, hub ratio, reciprocity, network replicate, nodes, arcs,
  commodities, demand seed, total demand, λ\* and Bonnet network;
- the settings: design and final spread, correlation, experiment, number of attempts
  for the base, seed, stream and scenarios requested;
- the results: scenarios saved, moment-matching trials and error, times of moment
  matching and of the feasibility checks, and file size;
- the realized statistics of the ratio demand / base demand: mean, standard deviation
  (average over commodities) and its target, minimum and maximum, and the average
  correlation between distinct demands.

`generation/experiments.csv` lists every experiment, including the two rejected attempts
for B03, with its spread, correlation, seed, stream, numbers of requested and feasible
scenarios, status, duration and completion time.

## Reproducing

Requirements:
- `bash` and `awk` for the shell scripts;
- Python 3 for `generate_set.py`, with `numpy` for `summarize.py` and `highspy` for `slack.py`;
- the `stoch_gen` executable;
- to regenerate the base instances only, a clone of the repository of Bonnet et al.
  and its Python dependencies.

1. **`stoch_gen`.** Build HiGHS 1.11.0 as a static library. Then build `stoch_gen` in
   the folder `program/` of https://github.com/larseeri/stoch-MCFNDP-gener, commit
   `280ee8c` or later, with:
   ```
   make HIGHS_SOLVER=1 HIGHS_HOME=<HiGHS installation> CXXFLAGS="-DHIGHS -std=gnu++17" LDFLAGS="-L<HiGHS installation>/lib -lhighs -lz"
   ```
   Put the executable on the `PATH`, or name it in the variable `STOCH_GEN`, e.g.
   `export STOCH_GEN=/path/to/stoch_gen`.
2. **Base instances** (optional, since they are included). Run
   `generation/make_bases.sh <clone of SSNDP-Realistic-Instance-Generator> [python]`.
   - The script runs the code of Bonnet et al. on each row of `generation/design.csv`, in
     a temporary copy, and converts the result with `bonnet2std.sh -F dow`.
   - It compares both files with `base/`, or writes them there if they are missing.
   - The code of Bonnet et al. requires Python 3.10. With Python 3.9, the script makes it
     run by adding `from __future__ import annotations`, which does not change the results.
   - The set was built with Python 3.9.6, numpy 2.0.2, scipy 1.13.1, networkx 3.2.1 and
     matplotlib 3.9.4.
3. **Scenarios.** Two ways:
   - **Replay the whole procedure:** run `python3 generation/generate_set.py --log replay.csv`.
     - The script applies the numbering of experiments, the seed and streams, the
       feasibility rule and the retries, and records every experiment in `replay.csv`.
     - Option `-j` sets the number of parallel runs (default 3). Run at most two or three
       at a time on a laptop. The original generation took about 90 minutes with three
       parallel runs on the setup above.
     - It resumes from its log if it is interrupted.
     - At the end, `diff <(cut -d, -f1-10 generation/experiments.csv) <(cut -d, -f1-10 replay.csv)`
       must report no difference. (Run without `--log`, the script finds every experiment
       of `generation/experiments.csv` completed and does nothing.)
   - **Regenerate one instance:** run `generation/run_one.sh <base> <instance>`. It runs
     `stoch_gen` on `scenarios/<base>/<instance>.par`, from the folder of the
     configuration.
   - **Rerun a rejected attempt:** from `scenarios/<base>/`, run
     `$STOCH_GEN +F rejected/<instance>.par`. Its outputs are written in
     `scenarios/<base>/`.

   Then check the regenerated files with `shasum -a 256 -c checksums.sha256`
   (`sha256sum -c` on Linux). Add `--ignore-missing` if only some instances were
   regenerated. Regeneration rewrites the logs, where only the time stamps change, and
   the files of target moments and correlations, which are identical.
4. **Manifest and λ\*.** Run `python3 generation/slack.py base/B*.txt > generation/slack.csv`,
   then, once the scenario files are regenerated, `python3 generation/summarize.py`.

## Notes

* **Exact reproduction.** All steps are deterministic. Rerunning an experiment
  reproduces its scenario file byte for byte. This was checked on experiments 3, 4, 5,
  10, 21 and 22, with the configurations and scripts of this folder for experiments 3,
  4, 5, 21 and 22. All the base instances were also regenerated byte for byte with
  `generation/make_bases.sh`. These checks used the build described above. Other
  compilers or platforms may round floating-point operations differently, so their
  output may differ in the last printed digits.
* **Version of `stoch_gen`.** Use commit `280ee8c` (3 October 2026) or later. Earlier
  versions compute the target variance of the uniform distribution as D²(α + β)/12
  instead of D²(α + β)²/12. The configurations set the number of trials explicitly
  (`t 20`), so the default of the version used does not matter.
* **File paths.** `stoch_gen` cannot read file names containing spaces, so the
  configurations use paths relative to their own folder and must be run from it
  (`generation/run_one.sh` does this).

## Sources

* L. Bonnet, S. Belieres, M. Hewitt and S. U. Ngueveu. A complex network analysis
  approach for generating realistic instances of the scheduled service network design
  problem. HAL open science hal-05345829, 2025. Code and networks:
  https://github.com/LBonnet159/SSNDP-Realistic-Instance-Generator (MIT License,
  reproduced in `LICENSE-SSNDP-Realistic-Instance-Generator.txt`).
* K. Høyland, M. Kaut and S. W. Wallace. A heuristic for moment-matching scenario
  generation. *Computational Optimization and Applications* 24(2):169–185, 2003.
* Stochastic generator: https://github.com/larseeri/stoch-MCFNDP-gener. Deterministic
  generator and format descriptions: https://github.com/larseeri/determ-MCFNDP-gener.
