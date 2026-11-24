"""Figure comparing several alignment sweeps, e.g. without surface, with surface, attack from below.

Usage, from the repository root:
    python3 -m scripts.plot_sweeps OUT.png LABEL=SWEEP.jsonl [LABEL=SWEEP.jsonl ...]

Each SWEEP.jsonl is written by scripts/sweep_alignment.py. The figure shows, against alignment
strength, the capture rate (mean +- standard error), the escape ratio R_esc (bootstrap over
runs, only for sweeps with the control school) and the polarization when the predator appears.
A table of the same numbers is printed.
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from scripts.plot_sweep import INK, MUTED, SURFACE, escape_ratio, load, style  # noqa: E402

# Okabe-Ito palette (colour-blind safe)
COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7"]

if __name__ == "__main__":
    out = sys.argv[1]
    sweeps = [arg.split("=", 1) for arg in sys.argv[2:]]
    fig, axes = plt.subplots(3, 1, figsize=(6.5, 8.5), sharex=True, facecolor=SURFACE)
    for (label, path), color in zip(sweeps, COLORS, strict=False):
        data = load(path)
        algs = np.array(sorted(data["capture_rate"]))
        rate = [np.array(data["capture_rate"][a]) for a in algs]
        pol = [np.array(data["polarization_start"][a]) for a in algs]
        mean = np.array([r.mean() for r in rate])
        se = np.array([r.std(ddof=1) / np.sqrt(len(r)) for r in rate])
        axes[0].errorbar(algs, mean, se, color=color, marker="o", ms=4, capsize=2, label=label)
        axes[2].plot(algs, [p.mean() for p in pol], color=color, marker="o", ms=4)
        print(
            f"\n{label} ({path})\nalg   | runs | capture rate      | R_esc           | polarization"
        )
        paired = all(len(data["capture_rate_nf"][a]) == len(data["capture_rate"][a]) for a in algs)
        esc = []
        for a, r, m, s, p in zip(algs, rate, mean, se, pol, strict=True):
            e = escape_ratio(r, np.array(data["capture_rate_nf"][a])) if paired else (np.nan,) * 2
            esc.append(e)
            print(
                f"{a:<5} | {len(r):4} | {m:.4f} ± {s:.4f} | {e[0]:6.3f} ± {e[1]:.3f} | {p.mean():.3f}"
            )
        if paired:
            esc = np.array(esc)
            axes[1].errorbar(algs, esc[:, 0], esc[:, 1], color=color, marker="o", ms=4, capsize=2)
        print(f"lowest capture rate at {algs[np.argmin(mean)]}")
    style(axes[0], "capture rate γc")
    style(axes[1], "escape ratio R_esc")
    style(axes[2], "polarization Φ")
    axes[1].axhline(0, color=INK, lw=0.8)
    axes[2].set_xlabel("alignment strength μ_alg", color=MUTED)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"\nsaved {out}")
