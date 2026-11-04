"""Compares the Python model with results of the original C++ code (issue #1).

Usage, from the repository root:
    python3 scripts/compare_cpp.py CPP_RESULTS.jsonl OUT.jsonl [n_seeds]

CPP_RESULTS.jsonl has one line per C++ run, with the keys "alg_strength",
"dead" (number of captures), "polarization_start" (when the predator appears)
and "polarisation_finale" (surviving prey at the end of the attack). The Python
runs use the alignment strengths of that file between 0.5 and 5 and write the
same quantities.
"""
import json
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, ".")
from predprey.params import Params  # noqa: E402
from predprey.simulation import polarization, simulate  # noqa: E402


def run(args):
    alg, seed = args
    r = simulate(Params(alg_strength=alg), seed=seed)
    u = r["prey_u"][-1]
    return {"alg_strength": alg, "seed": seed, "dead": len(r["kills"]),
            "polarisation_finale": float(polarization(u[~np.isnan(u[:, 0])])),
            "polarization_start": float(r["polarization_start"])}


def summary(rows, key):
    by_alg = defaultdict(list)
    for r in rows:
        by_alg[r["alg_strength"]].append(r[key])
    return {a: (np.mean(v), np.std(v)) for a, v in sorted(by_alg.items())}


if __name__ == "__main__":
    cpp_file, out_file = sys.argv[1], sys.argv[2]
    n_seeds = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    cpp = [json.loads(line) for line in open(cpp_file)]
    algs = sorted({r["alg_strength"] for r in cpp if 0.5 <= r["alg_strength"] <= 5})

    with Pool() as pool, open(out_file, "w") as f:
        py = []
        for r in pool.imap_unordered(run, [(a, s) for a in algs for s in range(n_seeds)]):
            py.append(r)
            f.write(json.dumps(r) + "\n")
            f.flush()

    cpp = [r for r in cpp if r["alg_strength"] in algs]
    print("mean ± std over runs")
    print("alg   | captures C++ | captures Py | start pol. C++ | start pol. Py | final pol. C++ | final pol. Py")
    cd, pd = summary(cpp, "dead"), summary(py, "dead")
    cs, ps = summary(cpp, "polarization_start"), summary(py, "polarization_start")
    cp, pp = summary(cpp, "polarisation_finale"), summary(py, "polarisation_finale")
    for a in algs:
        print(f"{a:<5} | {cd[a][0]:5.1f} ± {cd[a][1]:4.1f} | {pd[a][0]:5.1f} ± {pd[a][1]:4.1f}"
              f" | {cs[a][0]:6.3f} ± {cs[a][1]:5.3f} | {ps[a][0]:6.3f} ± {ps[a][1]:5.3f}"
              f" | {cp[a][0]:6.3f} ± {cp[a][1]:5.3f} | {pp[a][0]:6.3f} ± {pp[a][1]:5.3f}")
