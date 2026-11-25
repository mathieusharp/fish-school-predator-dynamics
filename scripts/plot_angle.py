"""Figure of the first attack against the attack direction (issue #14).

Usage, from the repository root:
    python3 -m scripts.plot_angle SWEEP.jsonl OUT.png [window]

SWEEP.jsonl is written by scripts/sweep_angle.py, with the capture times (kill_times).
Over 120 time units the predator attacks the school many times from any direction, so the
initial direction only shows in the first attack. For each alignment strength and attack
angle, the figure shows the time to the first capture (runs without capture count as the
whole attack time, so the mean is a lower bound) and the captures in the first `window` time
units (default 20), for the fleeing school and for the control school.
"""

import json
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from predprey.params import Params  # noqa: E402
from scripts.plot_sweep import MUTED, SURFACE, style  # noqa: E402
from scripts.sweep_alignment import TABLE_S1  # noqa: E402

COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7"]


def mean_se(x):
    x = np.asarray(x, dtype=float)
    return x.mean(), x.std(ddof=1) / np.sqrt(len(x))


if __name__ == "__main__":
    path, out = sys.argv[1], sys.argv[2]
    window = float(sys.argv[3]) if len(sys.argv) > 3 else 20.0
    attack_time = TABLE_S1.get("attack_time", Params().attack_time)
    runs = defaultdict(list)
    for line in open(path):
        r = json.loads(line)
        if "kill_times" not in r:
            sys.exit(f"{path} has no capture times: run the sweep again with the current code")
        runs[(r["alg_strength"], round(np.degrees(r["pred_angle"])))].append(r)
    algs = sorted({a for a, _ in runs})
    angles = sorted({g for _, g in runs})

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 6.5), sharex=True, facecolor=SURFACE)
    print(
        "alg   | angle | runs | first capture (fleeing, control) | captures in window (fl., ctrl)"
    )
    for a, color in zip(algs, COLORS, strict=False):
        rows = []
        for g in angles:
            rs = runs[(a, g)]
            first = [min(r["kill_times"], default=attack_time) for r in rs]
            first_nf = [min(r["kill_times_nf"], default=attack_time) for r in rs]
            early = [sum(t < window for t in r["kill_times"]) for r in rs]
            early_nf = [sum(t < window for t in r["kill_times_nf"]) for r in rs]
            rows.append([*mean_se(first), *mean_se(first_nf), *mean_se(early), *mean_se(early_nf)])
            f, fs, c, cs, e, es, n, ns = rows[-1]
            print(
                f"{a:<5} | {g:5} | {len(rs):4} | {f:5.1f} ± {fs:4.1f}, {c:5.1f} ± {cs:4.1f}"
                f"       | {e:4.2f} ± {es:4.2f}, {n:4.2f} ± {ns:4.2f}"
            )
        rows = np.array(rows)
        kw = dict(color=color, marker="o", ms=4, capsize=2)
        axes[0].errorbar(angles, rows[:, 0], rows[:, 1], label=f"μ_alg = {a}", **kw)
        axes[0].errorbar(angles, rows[:, 2], rows[:, 3], ls="--", alpha=0.5, **kw)
        axes[1].errorbar(angles, rows[:, 4], rows[:, 5], **kw)
        axes[1].errorbar(angles, rows[:, 6], rows[:, 7], ls="--", alpha=0.5, **kw)
    style(axes[0], "time to first capture")
    style(axes[1], f"captures in the first {window:g} time units")
    axes[1].set_xlabel("attack angle (0 = from the front, 180 = from behind)", color=MUTED)
    axes[1].set_xticks(angles)
    axes[0].legend(frameon=False, fontsize=8, title="solid: fleeing, dashed: control")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"saved {out}")
