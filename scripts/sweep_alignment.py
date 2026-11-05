"""Captures and polarization against alignment strength, to reproduce Fig 2 of the paper (issue #5).

Usage, from the repository root:
    python3 -m scripts.sweep_alignment OUT.jsonl ALG_STRENGTHS [n_runs]

ALG_STRENGTHS is a comma-separated list such as 0.5,1,2 or a range start:stop:step
such as 0:1.6:0.1 (stop included). Runs use the parameters of table S1 of the paper
and are spread over all cores. One JSON line is appended to OUT.jsonl per run, so an
interrupted sweep can be restarted with the same command: finished runs are skipped.
Run i of every alignment strength uses seed i.
"""

import json
import os
import sys
from multiprocessing import Pool

import numpy as np

from predprey.params import Params
from predprey.simulation import simulate

# Table S1 of Klamser & Romanczuk (2021), where it differs from the C++ defaults of Params
TABLE_S1 = {"attack_time": 120.0, "kill_rate": 1 / 3}


def parse_algs(text):
    if ":" in text:
        start, stop, step = (float(x) for x in text.split(":"))
        return [round(a, 6) for a in np.arange(start, stop + step / 2, step)]
    return [float(a) for a in text.split(",")]


def run(args):
    alg, seed = args
    p = Params(alg_strength=alg, **TABLE_S1)
    r = simulate(p, seed=seed)
    return {
        "alg_strength": alg,
        "seed": seed,
        "captures": len(r["kills"]),
        "capture_rate": len(r["kills"]) / p.attack_time,  # gamma_c of the paper
        "polarization_start": float(r["polarization_start"]),
        "cluster_size": r["cluster_size"],
    }


if __name__ == "__main__":
    out_file, algs = sys.argv[1], parse_algs(sys.argv[2])
    n_runs = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    done = set()
    if os.path.exists(out_file):
        for line in open(out_file):
            r = json.loads(line)
            done.add((r["alg_strength"], r["seed"]))
    todo = [(a, s) for s in range(n_runs) for a in algs if (a, s) not in done]
    print(f"{len(todo)} runs to do, {len(done)} already in {out_file}")

    os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)
    with Pool() as pool, open(out_file, "a") as f:
        for k, r in enumerate(pool.imap_unordered(run, todo), 1):
            f.write(json.dumps(r) + "\n")
            f.flush()
            print(
                f"{k}/{len(todo)}  alg {r['alg_strength']:<5} captures {r['captures']}", flush=True
            )
