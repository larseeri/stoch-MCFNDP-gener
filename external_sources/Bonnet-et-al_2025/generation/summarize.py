"""Builds manifest.csv for the stochastic instances.

Usage: python summarize.py [ROOT]      (ROOT defaults to the parent folder)

For every configuration scenarios/<base>/<id>.par, reads the stoch_gen log
and the stacked DOW scenario file <id>.txt, and reports the features of the
base instance (from generation/design.csv and generation/slack.csv), the
design and final spreads, the experiment that produced the instance and the
number of attempts for its base (from generation/experiments.csv), the
requested and
saved numbers of scenarios, the moment-matching error, and statistics of
the realized demands relative to the base demands: mean and standard
deviation of the ratio demand/base demand (target 1 and spread/sqrt(3) for
a uniform distribution), its extreme values, and the average correlation
between the demands of distinct commodities. It also reports, from the time
stamps of the log, the seconds spent in moment matching (HKW algorithm) and
in the feasibility checks of the scenarios. Requires numpy and the scenario
files (see the README for their regeneration).
"""
import csv
import datetime
import math
import os
import re
import sys

import numpy as np


def base_volumes(path):
    tok = open(path).read().split()
    n, m, k = int(tok[1]), int(tok[2]), int(tok[3])
    p = 4 + 7 * m
    return np.array([float(tok[p + 3 * h + 2]) for h in range(k)])


def scenario_demands(path, k):
    """Demands of every saved scenario, one row per scenario."""
    rows, cur, state = [], [], None
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "SCENARIO":
                if cur:
                    rows.append(cur)
                cur, state = [], "head"
            elif parts[0] == "MULTIGEN.DAT:":
                continue
            elif state == "head":
                state = "arcs"
            elif len(parts) == 3:
                cur.append(float(parts[2]))
    if cur:
        rows.append(cur)
    d = np.array(rows)
    assert d.shape[1] == k, (path, d.shape)
    return d


def log_seconds(log):
    """Seconds of moment matching and of feasibility checks, from the time stamps of a log."""
    stamps = [(m.start(), datetime.datetime.strptime(m.group(), "%Y-%m-%d %H:%M:%S"))
              for m in re.finditer(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d", log)]

    def after(label):
        pos = log.find(label)
        return next(t for p, t in stamps if p > pos)

    def before(label):
        pos = log.find(label)
        return [t for p, t in stamps if p < pos][-1]

    hkw = before("End generation of raw scenarios") - after("Begin generation of raw scenarios")
    check = after("Saved ") - after("Verifying feasibility")
    return int(hkw.total_seconds()), int(check.total_seconds())


def base_features(root):
    """Features of the base instances: design.csv plus, if present, slack.csv."""
    feats = {}
    for r in csv.DictReader(open(f"{root}/generation/design.csv")):
        feats[r["base_id"]] = {"density_pct": r["density"], "hub_ratio_pct": r["hub_ratio"],
                               "reciprocity_pct": r["reciprocity"], "network_replicate": r["network_replicate"],
                               "nodes": r["nodes"], "arcs": r["arcs"], "commodities": r["commodities"],
                               "demand_seed": r["demand_seed"], "bonnet_network": r["bonnet_network"]}
    for r in csv.DictReader(open(f"{root}/generation/assignment.csv")):
        feats[r["base_id"]]["design_spread"] = f"{int(r['spread_pct']) / 100:.2f}"
    if os.path.exists(f"{root}/generation/slack.csv"):
        for line in open(f"{root}/generation/slack.csv"):
            p = line.strip().split(",")
            feats[p[0][:-4]].update({"total_demand": p[4], "lambda_star": p[5]})
    return feats


def experiments(root):
    """Experiment number of each instance, and number of attempts of each base."""
    exp, attempts = {}, {}
    path = f"{root}/generation/experiments.csv"
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            exp[r["instance_id"]] = r["experiment"]
            attempts[r["base_id"]] = attempts.get(r["base_id"], 0) + 1
    return exp, attempts


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    feats = base_features(root)
    exp, attempts = experiments(root)
    out = []
    for base in sorted(d for d in os.listdir(f"{root}/scenarios")
                       if not d.startswith(".") and os.path.isdir(f"{root}/scenarios/{d}")):
        vol = base_volumes(f"{root}/base/{base}.dow")
        k = len(vol)
        folder = f"{root}/scenarios/{base}"
        for par in sorted(f for f in os.listdir(folder) if f.endswith(".par")):
            sid = par[:-4]
            opts = {(p := line.split(maxsplit=1))[0]: p[1].strip() if len(p) > 1 else ""
                    for line in open(f"{folder}/{par}") if line.strip()}
            spread, corr = float(opts["A"]), float(opts["XDD"])
            row = {"instance_id": sid, "base_id": base, **feats[base], "spread": f"{spread:.2f}",
                   "correlation": f"{corr:.1f}",
                   "experiment": exp.get(sid, ""), "attempts": attempts.get(base, ""),
                   "seed": opts["r"], "stream": opts["R"], "scenarios_requested": int(opts["n"])}
            log = open(f"{folder}/{sid}.log").read() if os.path.exists(f"{folder}/{sid}.log") else ""
            saved = re.search(r"Saved (\d+) feasible scenarios out of (\d+)", log)
            if not saved:
                row["status"] = "missing"
                out.append(row)
                continue
            if not os.path.exists(f"{folder}/{sid}.txt"):
                sys.exit(f"{folder}/{sid}.txt not found: regenerate the scenario files first (see the README)")
            errs = re.findall(r"after trial\s+(\d+):\s+([0-9.eE+-]+)", log)
            d = scenario_demands(f"{folder}/{sid}.txt", k)
            ratio = d / vol
            c = np.corrcoef(d, rowvar=False)
            off = c[~np.eye(k, dtype=bool)]
            hkw_s, check_s = log_seconds(log)
            row.update({
                "status": "ok", "scenarios_saved": int(saved.group(1)),
                "hkw_seconds": hkw_s, "check_seconds": check_s,
                "hkw_trials": int(errs[-1][0]) if errs else "", "hkw_error": float(errs[-1][1]) if errs else "",
                "file_MB": round(os.path.getsize(f"{folder}/{sid}.txt") / 1e6, 2),
                "ratio_mean": round(float(ratio.mean()), 4),
                "ratio_sd": round(float(ratio.std(axis=0).mean()), 4),
                "ratio_sd_target": round(spread / math.sqrt(3), 4),
                "ratio_min": round(float(ratio.min()), 4), "ratio_max": round(float(ratio.max()), 4),
                "corr_mean": round(float(off.mean()), 4)})
            assert row["scenarios_saved"] == d.shape[0], (sid, row["scenarios_saved"], d.shape)
            out.append(row)
    fields = ["instance_id", "base_id", "density_pct", "hub_ratio_pct", "reciprocity_pct", "network_replicate",
              "nodes", "arcs", "commodities", "demand_seed", "total_demand", "lambda_star", "bonnet_network",
              "design_spread", "spread", "correlation", "experiment", "attempts", "scenarios_requested",
              "scenarios_saved",
              "status", "hkw_trials", "hkw_error", "hkw_seconds", "check_seconds", "seed", "stream", "file_MB", "ratio_mean", "ratio_sd",
              "ratio_sd_target", "ratio_min", "ratio_max", "corr_mean"]
    with open(f"{root}/manifest.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out)
    done = sum(r["status"] == "ok" for r in out)
    print(f"{len(out)} configurations, {done} completed; manifest written to {root}/manifest.csv")


main()
