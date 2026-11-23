"""Captures against the attack direction relative to the school heading (fountain effect,
Bartashevich et al. 2024).

Usage, from the repository root:
    python3 -m scripts.sweep_angle OUT.jsonl ALG_STRENGTHS ANGLES [n_runs [dim [noise]]]

ALG_STRENGTHS as in scripts/sweep_alignment.py. ANGLES in degrees, same syntax: 0 means
the predator starts in front of the school, 180 behind it (pred_angle of Params). Each
angle gets a uniform noise of +-noise degrees (default 10) instead of the almost random
direction of the paper. Same output and restart as scripts/sweep_alignment.py, with the
keys pred_angle and pred_angle_noise (radians) added.
"""

import json
import os
import sys
from multiprocessing import Pool

import numpy as np

from scripts.sweep_alignment import code_version, parse_algs, run

if __name__ == "__main__":
    out_file, algs, angles = sys.argv[1], parse_algs(sys.argv[2]), parse_algs(sys.argv[3])
    n_runs = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    dim = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    noise = np.radians(float(sys.argv[6]) if len(sys.argv) > 6 else 10.0)
    done = set()
    if os.path.exists(out_file):
        for line in open(out_file):
            r = json.loads(line)
            done.add((r["alg_strength"], round(np.degrees(r["pred_angle"]), 6), r["seed"]))
    todo = [
        (a, s, dim, {"pred_angle": float(np.radians(g)), "pred_angle_noise": float(noise)})
        for s in range(n_runs)
        for g in angles
        for a in algs
        if (a, g, s) not in done
    ]
    version = code_version()
    print(f"{len(todo)} runs to do, {len(done)} already in {out_file}, code {version}")
    if version.endswith("-dirty"):
        print("warning: uncommitted changes, results may not be reproducible from a commit")

    os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)
    with Pool() as pool, open(out_file, "a") as f:
        for k, r in enumerate(pool.imap_unordered(run, todo), 1):
            f.write(json.dumps({**r, "commit": version}) + "\n")
            f.flush()
            angle = np.degrees(r["pred_angle"])
            print(
                f"{k}/{len(todo)}  alg {r['alg_strength']:<5} angle {angle:<5.0f} "
                f"captures {r['captures']}",
                flush=True,
            )
