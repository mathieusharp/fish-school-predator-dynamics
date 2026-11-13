"""Python model against the original C++ code, Table S1 parameters (issue #1).

Usage, from the repository root:
    python3 -m scripts.validate_cpp

Reads report/data/cpp_sweep.jsonl (scripts/run_cpp.py with --table-s1, 40 runs per
alignment strength) and report/data/sweep_control.jsonl (scripts/sweep_alignment.py) and
prints a Welch test per alignment strength. The C++ curves are drawn in Fig 1 of the report
by scripts/report_figures.py, which uses load and curves.
"""

import json
from collections import defaultdict

import numpy as np
from scipy.stats import ttest_ind

from scripts.report_figures import DATA
from scripts.sweep_alignment import TABLE_S1

N, BETA = 400, 2  # chi = beta N (<Phi^2> - <Phi>^2), beta = 1/(D v0), as in report_figures


def load(path, captures):
    """{alg: {"rate": [...], "pol": [...]}}, captures read from the key `captures`."""
    data = defaultdict(lambda: {"rate": [], "pol": []})
    for line in open(path):
        r = json.loads(line)
        a = round(r["alg_strength"], 1)
        data[a]["rate"].append(r[captures] / TABLE_S1["attack_time"])
        data[a]["pol"].append(r["polarization_start"])
    return {a: {k: np.array(v) for k, v in d.items()} for a, d in sorted(data.items())}


def curves(data):
    rate = [d["rate"] for d in data.values()]
    pol = [d["pol"] for d in data.values()]
    return {
        "rate": np.array([r.mean() for r in rate]),
        "rate_se": np.array([r.std(ddof=1) / np.sqrt(len(r)) for r in rate]),
        "pol": np.array([p.mean() for p in pol]),
        "pol_sd": np.array([p.std(ddof=1) for p in pol]),
        "chi": np.array([BETA * N * p.var() for p in pol]),
    }


def print_tests(cpp, py):
    print("alg  | gamma C++  Python  p     | Phi C++  Python  p")
    for a in cpp:
        pr = ttest_ind(cpp[a]["rate"], py[a]["rate"], equal_var=False).pvalue
        pp = ttest_ind(cpp[a]["pol"], py[a]["pol"], equal_var=False).pvalue
        print(
            f"{a:.1f}  | {cpp[a]['rate'].mean():.3f}      {py[a]['rate'].mean():.3f}   {pr:.3f} "
            f"| {cpp[a]['pol'].mean():.3f}    {py[a]['pol'].mean():.3f}   {pp:.3f}"
        )


if __name__ == "__main__":
    cpp = load(f"{DATA}cpp_sweep.jsonl", "dead")
    py = load(f"{DATA}sweep_control.jsonl", "captures")
    print_tests(cpp, py)
