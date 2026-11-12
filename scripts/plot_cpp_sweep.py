"""Python model against the original C++ code, Table S1 parameters (issue #1).

Usage, from the repository root:
    python3 -m scripts.plot_cpp_sweep

Reads report/data/cpp_sweep.jsonl (scripts/run_cpp.py with --table-s1, 40 runs per
alignment strength) and report/data/sweep_control.jsonl (scripts/sweep_alignment.py),
prints a Welch test per alignment strength and writes
    report/figures/cpp_comparison.pdf, .png   polarization, susceptibility, capture rate
"""

import json
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

from scripts.report_figures import BLUE, DATA, FIGURES, MUTED, ORANGE, panel
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


def figure(algs, cpp, py):
    line = {"marker": "o", "ms": 3.5, "lw": 1.5, "capsize": 2, "elinewidth": 0.8}
    style = {"C++": {"color": BLUE}, "Python": {"color": ORANGE, "marker": "s", "ls": "--"}}
    fig, ax = plt.subplots(1, 3, figsize=(7, 2.4), layout="constrained")
    for name, c in (("C++", cpp), ("Python", py)):
        s = line | style[name]
        ax[0].errorbar(algs, c["pol"], c["pol_sd"], label=name, **s)
        ax[1].plot(
            algs, c["chi"], **{k: v for k, v in s.items() if k not in ("capsize", "elinewidth")}
        )
        ax[2].errorbar(algs, c["rate"], c["rate_se"], **s)
    panel(ax[0], "A", "Polarization", r"$\Phi$ (± s.d.)")
    panel(ax[1], "B", "Susceptibility", r"$\chi$")
    panel(ax[2], "C", "Capture rate", r"capture rate $\gamma_c$ (± s.e.)")
    ax[0].legend(loc="upper left", frameon=False, fontsize=7)
    for a in ax:
        a.set_xlabel(r"alignment strength $\mu_\mathrm{alg}$")
        a.axvline(0.9, color=MUTED, lw=0.8, ls=":", zorder=0)  # transition in the paper
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIGURES}cpp_comparison.{ext}", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    cpp = load(f"{DATA}cpp_sweep.jsonl", "dead")
    py = load(f"{DATA}sweep_control.jsonl", "captures")
    print_tests(cpp, py)
    figure(np.array(list(cpp)), curves(cpp), curves(py))
