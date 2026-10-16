"""Generates the stochastic instances of the set, one per base instance.

Usage: python3 generate_set.py [-j JOBS] [--log FILE]     (from any folder)
  -j JOBS     number of parallel runs (default 3)
  --log FILE  log of the experiments (default: experiments.csv in this folder)

Policy (agreed for the set):
  * every run of stoch_gen is a generation experiment, numbered 1, 2, ...;
    all use the seed 50163844 (option r) and experiment e uses the stream
    47 + 100 (e - 1) (option R);
  * experiments 1 to 20 are the first attempts, in the order of design.csv,
    with the spread and correlation of assignment.csv;
  * an attempt is accepted only if all its scenarios are feasible; otherwise
    the spread alpha = beta is lowered by 0.05 and a new experiment is run;
    the retries of a round are numbered after all experiments of the previous
    round, in the order of design.csv.

Each finished experiment is appended to the log, from which the plan is
rebuilt when the script is restarted, so an interrupted generation resumes
where it stopped. The log experiments.csv of this folder records the
generation of the set and is complete; to replay the generation, give
another log file, e.g. --log replay.csv, then compare it with
experiments.csv (all columns but seconds and finished must agree). Accepted instances are in ../scenarios/<base>/; the
configuration and log of a rejected attempt are moved to
../scenarios/<base>/rejected/ and its other outputs are deleted (the
experiment can be reproduced from its configuration). The script runs
stoch_gen through run_one.sh, at low priority.
"""
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

GEN = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(GEN)
LOG = os.path.join(GEN, "experiments.csv")
SEED = 50163844
FIRST_STREAM, STREAM_STEP = 47, 100
SPREAD_STEP = 5                       # percent
FIELDS = ["experiment", "base_id", "instance_id", "spread", "correlation", "seed", "stream",
          "scenarios_requested", "scenarios_feasible", "status", "seconds", "finished"]


def stream(e):
    return FIRST_STREAM + STREAM_STEP * (e - 1)


def instance_id(base, v, c):
    return f"{base}_V{v}_C{c:02d}"


def write_par(base, k, v, c, e):
    sid = instance_id(base, v, c)
    lines = ["F D", f"I ../../base/{base}.dow", f"O ./{sid}.txt", "S 1",
             "K", f"XDD {c / 100:.1f}", f"C ./{sid}_corrs.dat",
             "G", "T U", f"A {v / 100:.2f}", f"B {v / 100:.2f}", f"M ./{sid}_moms.dat",
             f"n {2 * k}", "t 20", "i 50", "m 0.001", "c 0.001", "l 1",
             f"o ./{sid}_raw.dat", f"r {SEED}", f"R {stream(e)}", ""]
    with open(os.path.join(ROOT, "scenarios", base, sid + ".par"), "w") as f:
        f.write("\n".join(lines))
    return sid


def run(exp):
    """Runs one experiment and returns its log record."""
    e, base, k, v, c = exp
    sid = write_par(base, k, v, c, e)
    start = time.time()
    subprocess.run([os.path.join(GEN, "run_one.sh"), base, sid], check=False)
    folder = os.path.join(ROOT, "scenarios", base)
    text = open(os.path.join(folder, sid + ".log")).read() if os.path.exists(os.path.join(folder, sid + ".log")) else ""
    m = re.search(r"Saved (\d+) feasible scenarios out of (\d+)", text)
    feasible = int(m.group(1)) if m else ""
    requested = 2 * k
    if not m:
        status = "failed"
    elif feasible == requested:
        status = "accepted"
    else:
        status = "rejected"
        rej = os.path.join(folder, "rejected")
        os.makedirs(rej, exist_ok=True)
        for ext in (".par", ".log"):
            shutil.move(os.path.join(folder, sid + ext), os.path.join(rej, sid + ext))
        for ext in (".txt", "_moms.dat", "_corrs.dat", "_raw.dat"):
            p = os.path.join(folder, sid + ext)
            if os.path.exists(p):
                os.remove(p)
    return {"experiment": e, "base_id": base, "instance_id": sid, "spread": f"{v / 100:.2f}",
            "correlation": f"{c / 100:.1f}", "seed": SEED, "stream": stream(e),
            "scenarios_requested": requested, "scenarios_feasible": feasible, "status": status,
            "seconds": round(time.time() - start), "finished": time.strftime("%Y-%m-%d %H:%M:%S")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=3, help="number of parallel runs")
    ap.add_argument("--log", default=LOG, help="log of the experiments (default: %(default)s)")
    args = ap.parse_args()
    jobs, log = args.j, os.path.abspath(args.log)
    design = list(csv.DictReader(open(os.path.join(GEN, "design.csv"))))
    assign = {r["base_id"]: r for r in csv.DictReader(open(os.path.join(GEN, "assignment.csv")))}
    done = {}
    if os.path.exists(log):
        for r in csv.DictReader(open(log)):
            done[int(r["experiment"])] = r
    else:
        with open(log, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=FIELDS).writeheader()

    # Round 1: one experiment per base, in design order, with the design settings.
    current = [(i + 1, r["base_id"], int(r["commodities"]),
                int(assign[r["base_id"]]["spread_pct"]), int(assign[r["base_id"]]["correlation_pct"]))
               for i, r in enumerate(design)]
    rnd = 1
    while current:
        todo = [x for x in current if x[0] not in done]
        if todo:
            # largest instances first, to shorten the round
            size = {r["base_id"]: int(r["arcs"]) * int(r["commodities"]) ** 2 for r in design}
            todo.sort(key=lambda x: -size[x[1]])
            print(f"round {rnd}: {len(todo)} experiment(s) to run", flush=True)
            with ThreadPoolExecutor(max_workers=jobs) as pool:
                for fut in as_completed([pool.submit(run, x) for x in todo]):
                    rec = fut.result()
                    if rec["status"] != "failed":       # a failed run is not recorded, so it is retried
                        with open(log, "a", newline="") as f:
                            csv.DictWriter(f, fieldnames=FIELDS).writerow(rec)
                    done[rec["experiment"]] = {k: str(v) for k, v in rec.items()}
                    print(f"experiment {rec['experiment']:>3} {rec['instance_id']}: {rec['status']} "
                          f"({rec['scenarios_feasible']}/{rec['scenarios_requested']}, {rec['seconds']} s)",
                          flush=True)
        failed = [x for x in current if done[x[0]]["status"] == "failed"]
        if failed:
            sys.exit("stoch_gen failed in experiment(s) " + ", ".join(str(x[0]) for x in failed) +
                     "; see the logs in scenarios/<base>/")
        nxt = max(x[0] for x in current) + 1
        retry = []
        for x in sorted((x for x in current if done[x[0]]["status"] == "rejected"),
                        key=lambda x: [r["base_id"] for r in design].index(x[1])):
            e, base, k, v, c = x
            if v - SPREAD_STEP <= 0:
                sys.exit(f"{base}: infeasible scenarios even with a spread of {v}%")
            retry.append((nxt, base, k, v - SPREAD_STEP, c))
            nxt += 1
        current = retry
        rnd += 1
    rows = sorted(csv.DictReader(open(log)), key=lambda r: int(r["experiment"]))
    with open(log, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print("all instances accepted", flush=True)


if __name__ == "__main__":
    main()
