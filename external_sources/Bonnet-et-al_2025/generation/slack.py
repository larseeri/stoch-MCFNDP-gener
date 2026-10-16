"""Largest uniform demand multiplier that can be routed with all arcs open.

Usage: python slack.py BONNET_INSTANCE.txt...
For each SNDP instance (Bonnet format), solves
    max lambda  s.t.  a multicommodity flow routes lambda * d_k for every
                      commodity k within the arc capacities (all arcs open).
lambda* >= 1 means the base instance is feasible; if lambda* >= 1 + spread,
every scenario whose demands stay within base * (1 + spread) is feasible.
Prints: file, nodes, arcs, commodities, total demand, lambda*. Requires
highspy and numpy.
"""
import sys

import highspy
import numpy as np


def read(path):
    nodes, arcs, comms, sec = 0, [], [], None
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line[0].isalpha():
            if sec == "C":
                break
            sec = {"NODES": "N", "ARCS": "A", "COMMODITIES": "C"}[line.split(",")[0]]
            continue
        p = line.split(",")
        if sec == "N":
            nodes += 1
        elif sec == "A":
            arcs.append((int(p[1]), int(p[2]), float(p[5])))
        else:
            comms.append((int(p[1]), int(p[2]), float(p[3])))
    return nodes, arcs, comms


def lam_star(n, arcs, comms):
    m, k = len(arcs), len(comms)
    inf = highspy.kHighsInf
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("solver", "ipm")
    nv = m * k + 1                                    # x[a,h] then lambda
    h.addVars(nv, np.zeros(nv), np.full(nv, inf))
    h.changeColsCost(1, np.array([nv - 1], dtype=np.int32), np.array([-1.0]))
    out = [[] for _ in range(n)]
    inn = [[] for _ in range(n)]
    for a, (o, d, _) in enumerate(arcs):
        out[o].append(a)
        inn[d].append(a)
    lo, up, starts, index, value = [], [], [], [], []
    for hh, (o, d, q) in enumerate(comms):
        for j in range(n):
            idx = [a * k + hh for a in out[j]] + [a * k + hh for a in inn[j]]
            val = [1.0] * len(out[j]) + [-1.0] * len(inn[j])
            if j == o:
                idx.append(nv - 1)
                val.append(-q)
            elif j == d:
                idx.append(nv - 1)
                val.append(q)
            starts.append(len(index))
            lo.append(0.0)
            up.append(0.0)
            index += idx
            value += val
    for a, (_, _, cap) in enumerate(arcs):
        starts.append(len(index))
        lo.append(-inf)
        up.append(cap)
        index += [a * k + hh for hh in range(k)]
        value += [1.0] * k
    h.addRows(len(lo), np.array(lo), np.array(up), len(index), np.array(starts, dtype=np.int32),
              np.array(index, dtype=np.int32), np.array(value))
    h.run()
    status = h.modelStatusToString(h.getModelStatus())
    assert status == "Optimal", status
    return -h.getInfo().objective_function_value


for path in sys.argv[1:]:
    n, arcs, comms = read(path)
    print(f"{path.split('/')[-1]},{n},{len(arcs)},{len(comms)},{sum(q for *_, q in comms):.3f},"
          f"{lam_star(n, arcs, comms):.4f}", flush=True)
