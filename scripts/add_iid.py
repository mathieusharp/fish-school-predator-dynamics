"""Adds the inter-individual distance (IID) to a sweep written before it was recorded.

Usage, from the repository root:
    python3 -m scripts.add_iid SWEEP.jsonl

SWEEP.jsonl is written by scripts/sweep_alignment.py at commit c3addae or later (it has
"nn_distance_start"). For each run, the equilibration is run again with the same seed,
which gives the same school when the predator appears; the attack is skipped. The
polarization and nearest-neighbour distance found again must equal those of the file,
which checks that the school is the same. The file is then rewritten with "iid_start"
added to each line, other values unchanged. About 25 s per run, spread over all cores.
"""

import json
import os
import sys
from multiprocessing import Pool

from predprey.params import Params
from predprey.simulation import simulate
from scripts.sweep_alignment import TABLE_S1


def equilibrate(run):
    # attack of one time step: only the equilibration matters
    p = Params(alg_strength=run["alg_strength"], **(TABLE_S1 | {"attack_time": 0.02}))
    r = simulate(p, seed=run["seed"])
    for key in ("polarization_start", "nn_distance_start"):
        if float(r[key]) != run[key]:
            raise RuntimeError(f"different school for {run['alg_strength']}, {run['seed']}: {key}")
    return run | {"iid_start": float(r["iid_start"])}


if __name__ == "__main__":
    path = sys.argv[1]
    runs = [json.loads(line) for line in open(path)]
    if any("iid_start" in r for r in runs):
        sys.exit(f"{path} already has iid_start")
    print(f"{len(runs)} runs to equilibrate again", flush=True)
    done = []
    with Pool() as pool:
        for k, r in enumerate(pool.imap(equilibrate, runs), 1):
            done.append(r)
            if k % 20 == 0 or k == len(runs):
                print(f"{k}/{len(runs)}", flush=True)
    # same order as before; written to a new file first so that an interruption loses nothing
    with open(path + ".tmp", "w") as f:
        for r in done:
            f.write(json.dumps(r) + "\n")
    os.replace(path + ".tmp", path)
    print(f"iid_start added to {path}")
