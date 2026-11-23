"""Captures and polarization against alignment strength, to reproduce Fig 2 of the paper (issue #5).

Usage, from the repository root:
    python3 -m scripts.sweep_alignment OUT.jsonl ALG_STRENGTHS [n_runs [dim [NAME=VALUE ...]]]

ALG_STRENGTHS is a comma-separated list such as 0.5,1,2 or a range start:stop:step
such as 0:1.6:0.1 (stop included). Runs use the parameters of table S1 of the paper
and are spread over all cores. One JSON line is appended to OUT.jsonl per run, so an
interrupted sweep can be restarted with the same command: finished runs are skipped.
Run i of every alignment strength uses seed i. dim is 2 (default) or 3 (issue #10).
Each line also records the commit of the code ("-dirty" if predprey/ or scripts/ were modified).
NAME=VALUE sets other fields of Params (VALUE in JSON: 2.5, true, false), for example
surface=true pred_below=true; they are saved in each line. Use one OUT.jsonl per setting.

Each run also attacks a non-fleeing copy of the school (control school, issue #4), which
gives the escape ratio of Fig 2D, and records the structure of the school when the predator
appears: inter-individual distance (mean over all pairs, inset of Fig 2C) and mean
nearest-neighbour distance. Files written by older versions lack some of these keys: write
to a new file rather than completing them (scripts/add_iid.py adds the inter-individual
distance to a file that has all the others).
"""

import json
import os
import subprocess
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


def code_version():
    """Commit of the code that runs the sweep, with -dirty if predprey/ or scripts/ are modified."""
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "predprey", "scripts"])
    return commit if clean.returncode == 0 else commit + "-dirty"


def run(args):
    """args: (alg_strength, seed, dim) or (alg_strength, seed, dim, other Params values)."""
    alg, seed, dim, *extra = args
    extra = extra[0] if extra else {}
    p = Params(dim=dim, alg_strength=alg, control=True, **TABLE_S1, **extra)
    r = simulate(p, seed=seed)
    return {
        **extra,
        "alg_strength": alg,
        "seed": seed,
        "dim": dim,
        "captures": len(r["kills"]),
        "capture_rate": len(r["kills"]) / p.attack_time,  # gamma_c of the paper
        "captures_nf": len(r["kills_nf"]),  # same school, not fleeing (issue #4)
        "capture_rate_nf": len(r["kills_nf"]) / p.attack_time,  # gamma_c,NF
        "polarization_start": float(r["polarization_start"]),
        "iid_start": float(r["iid_start"]),  # inset of Fig 2C
        "nn_distance_start": float(r["nn_distance_start"]),
        "cluster_size": r["cluster_size"],
    }


if __name__ == "__main__":
    out_file, algs = sys.argv[1], parse_algs(sys.argv[2])
    n_runs = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    dim = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    extra = {}
    for arg in sys.argv[5:]:
        name, value = arg.split("=", 1)
        if name not in Params.__dataclass_fields__:
            sys.exit(f"unknown parameter {name}")
        extra[name] = json.loads(value)
    done = set()
    if os.path.exists(out_file):
        for line in open(out_file):
            r = json.loads(line)
            if any(k not in r for k in ("captures_nf", "nn_distance_start", "iid_start")):
                sys.exit(f"{out_file} was written by an older version: use a new file")
            if any(r.get(name) != value for name, value in extra.items()):
                sys.exit(f"{out_file} was written with other parameters: use a new file")
            done.add((r["alg_strength"], r["seed"]))
    todo = [(a, s, dim, extra) for s in range(n_runs) for a in algs if (a, s) not in done]
    version = code_version()
    print(f"{len(todo)} runs to do, {len(done)} already in {out_file}, code {version}")
    if version.endswith("-dirty"):
        print("warning: uncommitted changes, results may not be reproducible from a commit")

    os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)
    with Pool() as pool, open(out_file, "a") as f:
        for k, r in enumerate(pool.imap_unordered(run, todo), 1):
            f.write(json.dumps({**r, "commit": version}) + "\n")
            f.flush()
            print(
                f"{k}/{len(todo)}  alg {r['alg_strength']:<5} captures {r['captures']}", flush=True
            )
