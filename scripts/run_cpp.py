"""Runs the original C++ code and writes one JSON line per run (issue #1).

Usage, from the repository root (see README, "Running the C++ code"):
    python3 scripts/run_cpp.py CPP_DIR OUT.jsonl ALG_STRENGTHS [n_runs]

CPP_DIR is the compiled clone of https://github.com/PaPeK/PredatorPrey,
ALG_STRENGTHS a comma-separated list such as 0.5,1,2,5. Lines are appended to
OUT.jsonl with the keys read by compare_cpp.py and plot_comparison.py.
Runs are sequential: the C++ code always writes to the same out_xx.h5.
The C++ seed comes from the clock, so runs cannot be replayed.
"""
import json
import os
import subprocess
import sys

import h5py
import numpy as np


def polarization(vel):
    u = vel / np.linalg.norm(vel, axis=1, keepdims=True)
    return float(np.linalg.norm(u.mean(axis=0)))


def run(alg):
    import SwarmDynByPy as swarmPy  # from CPP_DIR
    dic = swarmPy.get_base_params(200, 60)  # equilibration 200, attack 60 time units
    dic["alg_strength"] = alg
    swarmPy.heteroPopulation(dic)
    swarmPy.solve_parameter_dependencies(dic)
    subprocess.run("make cl > /dev/null", shell=True, check=True)
    subprocess.run(swarmPy.dic2swarmdyn_command(dic) + " > /dev/null", shell=True, check=True)
    with h5py.File("out_xx.h5", "r") as f:
        part = np.array(f["/part"])  # (frames, N, [x, y, vx, vy, ...])
    pos, vel = part[..., :2], part[..., 2:4]
    alive = np.any(pos != 0, axis=-1)  # dead prey and prey outside the attacked cluster are at (0, 0)
    return {"alg_strength": alg,
            "dead": int(alive[0].sum() - alive[-1].sum()),
            "polarization_start": polarization(vel[0][alive[0]]),
            "polarisation_finale": polarization(vel[-1][alive[-1]])}


if __name__ == "__main__":
    cpp_dir, out_file = sys.argv[1], os.path.abspath(sys.argv[2])
    algs = [float(a) for a in sys.argv[3].split(",")]
    n_runs = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    os.chdir(cpp_dir)
    sys.path.insert(0, ".")
    with open(out_file, "a") as f:
        for alg in algs:
            for _ in range(n_runs):
                r = run(alg)
                print(json.dumps(r), flush=True)
                f.write(json.dumps(r) + "\n")
