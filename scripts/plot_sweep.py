"""Figure of the alignment sweep: is the capture rate minimal at the order transition? (issue #5)

Usage, from the repository root:
    python3 -m scripts.plot_sweep SWEEP.jsonl OUT.png

SWEEP.jsonl is written by scripts/sweep_alignment.py. The transition is placed at the
maximum of the susceptibility chi = N * variance of the polarization over runs.
If the sweep has the control school (issue #4), the figure also shows its capture rate
gamma_c,NF and the escape ratio R_esc = 1 - gamma_c / gamma_c,NF of Fig 2D, with a
bootstrap over runs (each run and its control share the same starting state).
If it has the nearest-neighbour distance when the predator appears, the figure shows it
and the correlation between this distance and the capture rate (inset of Fig 2C).
"""

import json
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
BLUE, ORANGE = "#2a78d6", "#eb6834"
PAPER_TRANSITION = 0.9  # mu_alg,c at D = 0.5 in the paper


def load(path):
    data = defaultdict(lambda: defaultdict(list))
    for line in open(path):
        r = json.loads(line)
        for key in (
            "capture_rate",
            "capture_rate_nf",
            "polarization_start",
            "cluster_size",
            "nn_distance_start",
        ):
            if key in r:
                data[key][r["alg_strength"]].append(r[key])
    return data


def escape_ratio(rate, rate_nf, n_boot=2000):
    """R_esc = 1 - mean(gamma_c) / mean(gamma_c,NF), and its bootstrap standard error."""
    boot = np.random.default_rng(0).integers(len(rate), size=(n_boot, len(rate)))
    samples = 1 - rate[boot].mean(axis=1) / rate_nf[boot].mean(axis=1)
    return 1 - rate.mean() / rate_nf.mean(), samples.std(ddof=1)


def style(ax, ylabel):
    ax.set_facecolor(SURFACE)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.tick_params(colors=MUTED, labelsize=8)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color("#c8c7c0")


if __name__ == "__main__":
    data = load(sys.argv[1])
    algs = np.array(sorted(data["capture_rate"]))
    rate = [np.array(data["capture_rate"][a]) for a in algs]
    pol = [np.array(data["polarization_start"][a]) for a in algs]
    size = [np.mean(data["cluster_size"][a]) for a in algs]
    chi = np.array([n * np.var(p, ddof=1) for n, p in zip(size, pol, strict=True)])
    transition = algs[np.argmax(chi)]
    lowest = algs[np.argmin([r.mean() for r in rate])]
    # control school only if every run has it, so that runs and controls stay paired
    control = all(
        len(data["capture_rate_nf"][a]) == len(r) for a, r in zip(algs, rate, strict=True)
    )
    if control:
        rate_nf = [np.array(data["capture_rate_nf"][a]) for a in algs]
        esc = np.array([escape_ratio(r, n) for r, n in zip(rate, rate_nf, strict=True)])
    distance = all(
        len(data["nn_distance_start"][a]) == len(r) for a, r in zip(algs, rate, strict=True)
    )
    if distance:
        dist = [np.array(data["nn_distance_start"][a]) for a in algs]
        # Pearson correlation with the capture rate, over the points and over all runs
        r_points = np.corrcoef([r.mean() for r in rate], [d.mean() for d in dist])[0, 1]
        r_runs = np.corrcoef(np.concatenate(rate), np.concatenate(dist))[0, 1]

    print(
        "alg   | runs | capture rate      | polarization  | chi    " + control * "| no flee | R_esc"
    )
    for k, (a, r, p, c) in enumerate(zip(algs, rate, pol, chi, strict=True)):
        se = r.std(ddof=1) / np.sqrt(len(r))
        line = f"{a:<5} | {len(r):4} | {r.mean():.4f} ± {se:.4f} | {p.mean():.3f} ± {p.std():.3f} | {c:6.2f}"
        if control:
            line += f" | {rate_nf[k].mean():.4f}  | {esc[k, 0]:.3f} ± {esc[k, 1]:.3f}"
        print(line)
    print(f"susceptibility maximum at {transition:g}, lowest capture rate at {lowest:g}")
    if control:
        print(f"highest escape ratio at {algs[np.argmax(esc[:, 0])]:g}")
    if distance:
        print(
            "correlation between capture rate and nearest-neighbour distance: "
            f"R = {r_points:.2f} over the points, {r_runs:.2f} over all runs"
        )

    n_panels = 3 + control + distance
    fig, axes = plt.subplots(
        n_panels, 1, figsize=(6.5, 2.7 * n_panels), sharex=True, facecolor=SURFACE
    )
    se = [r.std(ddof=1) / np.sqrt(len(r)) for r in rate]
    axes[0].errorbar(
        algs, [r.mean() for r in rate], yerr=se, color=ORANGE, marker="o", ms=5, lw=2, capsize=3
    )
    if control:
        se_nf = [r.std(ddof=1) / np.sqrt(len(r)) for r in rate_nf]
        axes[0].errorbar(
            algs,
            [r.mean() for r in rate_nf],
            yerr=se_nf,
            color=MUTED,
            marker="s",
            ms=4,
            lw=1.5,
            ls="--",
            capsize=3,
            label="control school (no fleeing)",
        )
        axes[0].plot([], [], color=ORANGE, marker="o", lw=2, label="fleeing school")
        axes[0].legend(fontsize=7, frameon=False)
    axes[0].set_title("Capture rate (captures per time unit)", fontsize=10, loc="left", color=INK)
    style(axes[0], "capture rate γc")
    axes[1].errorbar(
        algs,
        [p.mean() for p in pol],
        yerr=[p.std() for p in pol],
        color=BLUE,
        marker="o",
        ms=5,
        lw=2,
        capsize=3,
    )
    axes[1].set_title("Polarization when the predator appears", fontsize=10, loc="left", color=INK)
    style(axes[1], "polarization Φ")
    axes[2].plot(algs, chi, color=BLUE, marker="o", ms=5, lw=2)
    axes[2].set_title(
        "Susceptibility (N × variance of Φ over runs)", fontsize=10, loc="left", color=INK
    )
    style(axes[2], "χ")
    if control:
        axes[3].errorbar(
            algs, esc[:, 0], yerr=esc[:, 1], color=ORANGE, marker="o", ms=5, lw=2, capsize=3
        )
        axes[3].set_title(
            "Escape ratio R_esc = 1 − γc / γc,NF (Fig 2D)", fontsize=10, loc="left", color=INK
        )
        style(axes[3], "R_esc")
        axes[3].axhline(0, color=MUTED, lw=1)
    if distance:
        axes[-1].errorbar(
            algs,
            [d.mean() for d in dist],
            yerr=[d.std() for d in dist],
            color=BLUE,
            marker="o",
            ms=5,
            lw=2,
            capsize=3,
        )
        axes[-1].set_title(
            f"Nearest-neighbour distance when the predator appears (R with γc = {r_points:.2f})",
            fontsize=10,
            loc="left",
            color=INK,
        )
        style(axes[-1], "distance")
    axes[-1].set_xlabel("alignment strength μ_alg", color=MUTED)
    for ax in axes:
        ax.axvline(transition, color="#c03fc0", lw=1.5, ls="-.")
        ax.axvline(PAPER_TRANSITION, color=MUTED, lw=1, ls=":")
    axes[0].text(
        transition,
        1.02,
        f" χ max = {transition:g}",
        transform=axes[0].get_xaxis_transform(),
        color="#c03fc0",
        fontsize=8,
    )
    axes[0].text(
        PAPER_TRANSITION,
        0.03,
        f" paper ≈ {PAPER_TRANSITION}",
        transform=axes[0].get_xaxis_transform(),
        color=MUTED,
        fontsize=8,
    )
    n_runs = min(len(r) for r in rate)
    fig.text(
        0.01,
        0.01,
        f"Mean ± standard error ({n_runs}+ runs per point; polarization, distance: ± standard deviation; "
        "R_esc: bootstrap over runs).\n"
        "N = 400, D = 0.5, parameters of table S1.",
        fontsize=7,
        color=MUTED,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(sys.argv[2], dpi=150)
