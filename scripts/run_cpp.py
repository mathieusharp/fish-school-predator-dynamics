"""Runs the original C++ code and writes one JSON line per run (issue #1).

Usage, from the repository root (see README, "Running the C++ code"):
    python3 scripts/run_cpp.py CPP_DIR OUT.jsonl ALG_STRENGTHS [n_runs] [--table-s1]

CPP_DIR is the compiled clone of https://github.com/PaPeK/PredatorPrey,
ALG_STRENGTHS a comma-separated list such as 0.5,1,2,5. Lines are appended to
OUT.jsonl with the keys read by compare_cpp.py and plot_comparison.py.
Without --table-s1 the attack lasts 60 time units with kill rate 1 (defaults of
the C++ code); with it, 120 time units and kill rate 1/3 (Table S1 of the paper,
as in scripts/sweep_alignment.py).
Runs are spread over all cores, each worker in its own copy of CPP_DIR, since
the C++ code always writes to the same out_xx.h5.
The C++ seed comes from the clock, so runs cannot be replayed.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import h5py
import numpy as np

TABLE_S1 = {"attack_time": 120, "kill_rate": 1 / 3}
DEFAULTS = {"attack_time": 60, "kill_rate": 1}


def polarization(vel):
    u = vel / np.linalg.norm(vel, axis=1, keepdims=True)
    return float(np.linalg.norm(u.mean(axis=0)))


def init_worker(cpp_dir):
    work = tempfile.mkdtemp(prefix="predatorprey_")
    shutil.copytree(cpp_dir, work, dirs_exist_ok=True)
    os.chdir(work)
    sys.path.insert(0, work)
    # forked workers share NumPy's state: without this they would draw the same
    # individual alignment strengths (heteroPopulation)
    np.random.seed()


def run(args):
    alg, params = args
    import SwarmDynByPy as swarmPy  # from the copy of CPP_DIR

    dic = swarmPy.get_base_params(200, params["attack_time"])  # equilibration 200
    dic["alg_strength"] = alg
    dic["kill_rate"] = params["kill_rate"]
    swarmPy.heteroPopulation(dic)
    swarmPy.solve_parameter_dependencies(dic)
    subprocess.run("make cl > /dev/null", shell=True, check=True)
    subprocess.run(swarmPy.dic2swarmdyn_command(dic) + " > /dev/null", shell=True, check=True)
    with h5py.File("out_xx.h5", "r") as f:
        part = np.array(f["/part"])  # (frames, N, [x, y, vx, vy, ...])
    pos, vel = part[..., :2], part[..., 2:4]
    # dead prey and prey outside the attacked cluster are at (0, 0)
    alive = np.any(pos != 0, axis=-1)
    return {
        "alg_strength": alg,
        "dead": int(alive[0].sum() - alive[-1].sum()),
        "polarization_start": polarization(vel[0][alive[0]]),
        "polarisation_finale": polarization(vel[-1][alive[-1]]),
    }


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--table-s1"]
    params = TABLE_S1 if "--table-s1" in sys.argv else DEFAULTS
    cpp_dir, out_file = os.path.abspath(args[0]), os.path.abspath(args[1])
    algs = [float(a) for a in args[2].split(",")]
    n_runs = int(args[3]) if len(args) > 3 else 4
    tasks = [(alg, params) for alg in algs for _ in range(n_runs)]
    with Pool(initializer=init_worker, initargs=(cpp_dir,)) as pool, open(out_file, "a") as f:
        for r in pool.imap_unordered(run, tasks):
            print(json.dumps(r), flush=True)
            f.write(json.dumps(r) + "\n")
            f.flush()
