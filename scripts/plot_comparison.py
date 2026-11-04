"""Figure and Welch tests comparing the Python model with the C++ code (issue #1).

Usage, from the repository root:
    python3 scripts/plot_comparison.py CPP.jsonl PY.jsonl OUT.png

Both files are written by scripts/compare_cpp.py (Python) and by the C++ runner
(one JSON line per run, keys "alg_strength", "dead", "polarization_start").
"""
import json
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import ttest_ind  # noqa: E402

SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
COLORS = {"C++": "#2a78d6", "Python": "#eb6834"}
MARKERS = {"C++": "o", "Python": "s"}
KEYS = [("dead", "Captures in 60 time units", "captures"),
        ("polarization_start", "Polarization when the predator appears", "polarization")]


def load(path):
    data = defaultdict(lambda: defaultdict(list))
    for line in open(path):
        r = json.loads(line)
        for key, _, _ in KEYS:
            data[key][r["alg_strength"]].append(r[key])
    return data


if __name__ == "__main__":
    runs = {"C++": load(sys.argv[1]), "Python": load(sys.argv[2])}
    algs = sorted(set(runs["C++"]["dead"]) & set(runs["Python"]["dead"]))

    # Welch t-test: p = probability of a gap at least this large if both codes were the same
    print("alg   | " + " | ".join(f"{label:>12} p" for _, _, label in KEYS))
    for a in algs:
        ps = [ttest_ind(runs["C++"][key][a], runs["Python"][key][a], equal_var=False).pvalue
              for key, _, _ in KEYS]
        print(f"{a:<5} | " + " | ".join(f"{p:14.3f}" for p in ps))

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), facecolor=SURFACE)
    for ax, (key, title, ylabel) in zip(axes, KEYS):
        ax.set_facecolor(SURFACE)
        for k, (code, data) in enumerate(runs.items()):
            x = np.array(algs) * (1 + 0.025 * (2 * k - 1))  # small shift so error bars do not overlap
            mean = [np.mean(data[key][a]) for a in algs]
            se = [np.std(data[key][a], ddof=1) / np.sqrt(len(data[key][a])) for a in algs]
            n = len(data[key][algs[0]])
            ax.errorbar(x, mean, yerr=se, color=COLORS[code], marker=MARKERS[code], ms=8, lw=2,
                        capsize=3, markeredgecolor=SURFACE, markeredgewidth=1.5,
                        label=f"{code} ({n} runs)")
        ax.axvline(0.83, color=MUTED, lw=1, ls=":")
        ax.text(0.85, 0.03, "transition ≈ 0.83", transform=ax.get_xaxis_transform(),
                color=MUTED, fontsize=8)
        ax.set_xscale("log")
        ax.set_xticks(algs)
        ax.set_xticklabels([f"{a:g}" for a in algs])
        ax.minorticks_off()
        ax.set_xlabel("alignment strength (log scale)", color=MUTED)
        ax.set_ylabel(ylabel, color=MUTED)
        ax.set_title(title, fontsize=10, loc="left", color=INK)
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.tick_params(colors=MUTED, labelsize=8)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
        for s in ["left", "bottom"]:
            ax.spines[s].set_color("#c8c7c0")
    axes[0].set_ylim(bottom=0)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK, loc="lower right")
    fig.text(0.01, 0.01, "Mean ± standard error over runs. N = 400 prey, default parameters of "
             "Klamser & Romanczuk (2021).", fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(sys.argv[3], dpi=150)
