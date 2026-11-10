"""Figures, summary table and key numbers for the reports (issue #7).

Usage, from the repository root:
    python3 -m scripts.report_figures [--snapshots]

Reads report/data/sweep_control.jsonl (written by scripts/sweep_alignment.py) and writes
    report/figures/sweep.pdf, .png      Fig 2 of the paper at D = 0.5, six panels
    report/data/sweep_summary.csv       one line per alignment strength
and prints the numbers quoted in report/RESULTS.md. Takes a few seconds.
With --snapshots it also runs three simulations (about 2 min) and writes
    report/figures/snapshots.pdf, .png  the school when the predator appears, in the
                                        three phases
"""

import csv
import sys
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import ttest_ind  # noqa: E402

from predprey.params import Params  # noqa: E402
from predprey.simulation import simulate  # noqa: E402
from scripts.plot_sweep import PAPER_TRANSITION, escape_ratio, load  # noqa: E402
from scripts.sweep_alignment import TABLE_S1  # noqa: E402

DATA, FIGURES = "report/data/", "report/figures/"
INK, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#e6e5e0", "#c8c7c0"
BLUE, ORANGE, MAGENTA, RED = "#2a78d6", "#eb6834", "#c03fc0", "#d62728"
# alignment strengths grouped by phase, as in issue #5
PHASES = {"disordered": (0.0, 0.2), "transition": (0.6, 0.8), "ordered": (1.4, 1.6)}

plt.rcParams.update(
    {
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,  # editable text in the PDF
    }
)


def summarize(data, p):
    """One dict per alignment strength, with the quantities of Fig 2."""
    beta = 1 / (p.Dphi * p.speed0)  # chi = beta N (<Phi^2> - <Phi>^2), Eq S44d of S1 Text
    rows = []
    for a in sorted(data["capture_rate"]):
        rate = np.array(data["capture_rate"][a])
        rate_nf = np.array(data["capture_rate_nf"][a])
        pol = np.array(data["polarization_start"][a])
        dist = np.array(data["nn_distance_start"][a])
        n = np.mean(data["cluster_size"][a])
        esc, esc_se = escape_ratio(rate, rate_nf)
        rows.append(
            {
                "alg_strength": a,
                "runs": len(rate),
                "capture_rate": rate.mean(),
                "capture_rate_se": rate.std(ddof=1) / np.sqrt(len(rate)),
                "capture_rate_nf": rate_nf.mean(),
                "capture_rate_nf_se": rate_nf.std(ddof=1) / np.sqrt(len(rate)),
                "escape_ratio": esc,
                "escape_ratio_se": esc_se,
                "polarization": pol.mean(),
                "polarization_sd": pol.std(ddof=1),
                "susceptibility": beta * n * pol.var(),
                "nn_distance": dist.mean(),
                "nn_distance_sd": dist.std(ddof=1),
            }
        )
    return rows


def phase_values(data, key, phase):
    lo, hi = PHASES[phase]
    return np.concatenate([v for a, v in data[key].items() if lo - 1e-9 <= a <= hi + 1e-9])


def print_numbers(data, rows):
    col = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    algs = col["alg_strength"]
    print(f"runs: {int(col['runs'].sum())} ({int(col['runs'].min())} per alignment strength)")
    print(f"susceptibility maximum at mu_alg = {algs[np.argmax(col['susceptibility'])]:g}")
    print(
        f"polarization: {col['polarization'][0]:.2f} at 0 to {col['polarization'][-1]:.2f} at 1.6"
    )
    print(f"lowest capture rate at mu_alg = {algs[np.argmin(col['capture_rate'])]:g}")
    print(f"highest escape ratio at mu_alg = {algs[np.argmax(col['escape_ratio'])]:g}")
    print(f"largest nn distance at mu_alg = {algs[np.argmax(col['nn_distance'])]:g}")
    print("\nby phase (all runs of the phase pooled):")
    for phase in PHASES:
        rate = phase_values(data, "capture_rate", phase)
        nf = phase_values(data, "capture_rate_nf", phase)
        dist = phase_values(data, "nn_distance_start", phase)
        print(
            f"  {phase:<10} mu_alg {PHASES[phase][0]}-{PHASES[phase][1]}: "
            f"gamma_c {rate.mean():.4f}, gamma_c,NF {nf.mean():.4f}, "
            f"R_esc {1 - rate.mean() / nf.mean():.3f}, nn distance {dist.mean():.3f}"
        )
    trans = phase_values(data, "capture_rate", "transition")
    for phase in ("disordered", "ordered"):
        other = phase_values(data, "capture_rate", phase)
        t = ttest_ind(trans, other, equal_var=False)
        print(
            f"  transition vs {phase}: {100 * (1 - trans.mean() / other.mean()):.1f} % fewer "
            f"captures, Welch t = {t.statistic:.2f}, p = {t.pvalue:.1g}"
        )
    rate_all = np.concatenate([data["capture_rate"][a] for a in algs])
    nf_all = np.concatenate([data["capture_rate_nf"][a] for a in algs])
    t = ttest_ind(rate_all, nf_all, equal_var=False)
    print(
        f"  all runs: fleeing {rate_all.mean():.4f} vs not fleeing {nf_all.mean():.4f}, "
        f"R_esc {1 - rate_all.mean() / nf_all.mean():.3f}, Welch p = {t.pvalue:.1g}"
    )
    dist_all = np.concatenate([data["nn_distance_start"][a] for a in algs])
    print(
        "\nPearson correlation between capture rate and nn distance: "
        f"R = {np.corrcoef(col['capture_rate'], col['nn_distance'])[0, 1]:.2f} over the "
        f"{len(algs)} means, {np.corrcoef(rate_all, dist_all)[0, 1]:.2f} over the "
        f"{len(rate_all)} runs"
    )


def write_csv(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})


def mark_transition(ax, transition):
    ax.axvline(transition, color=MAGENTA, lw=1, ls="-.", zorder=0)
    ax.axvline(PAPER_TRANSITION, color=MUTED, lw=0.8, ls=":", zorder=0)


def panel(ax, letter, title, ylabel):
    ax.set_title(f"{letter}  {title}", loc="left", color=INK, fontweight="bold")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def sweep_figure(rows, data):
    col = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    algs = col["alg_strength"]
    transition = algs[np.argmax(col["susceptibility"])]
    line = {"marker": "o", "ms": 3.5, "lw": 1.5, "capsize": 2, "elinewidth": 0.8}
    fig, ax = plt.subplots(3, 2, figsize=(7, 7), layout="constrained")
    ax = ax.ravel()

    ax[0].errorbar(algs, col["polarization"], col["polarization_sd"], color=BLUE, **line)
    panel(ax[0], "A", "Polarization when the predator appears", "polarization Φ")
    ax[1].plot(algs, col["susceptibility"], color=BLUE, marker="o", ms=3.5, lw=1.5)
    panel(ax[1], "B", "Susceptibility", "χ")
    # the two vertical lines are named once, where panel A is empty
    xaxis = ax[0].get_xaxis_transform()
    ax[0].text(transition, 0.95, "χ max ", color=MAGENTA, transform=xaxis, ha="right")
    ax[0].text(PAPER_TRANSITION, 0.05, " paper ≈ 0.9", color=MUTED, transform=xaxis)

    ax[2].errorbar(
        algs,
        col["capture_rate_nf"],
        col["capture_rate_nf_se"],
        color=MUTED,
        ls="--",
        label="not fleeing (control), γc,NF",
        **(line | {"marker": "s"}),
    )
    ax[2].errorbar(
        algs, col["capture_rate"], col["capture_rate_se"], color=ORANGE, label="fleeing, γc", **line
    )
    ax[2].set_ylim(0.11, None)  # room for the legend below the curves
    ax[2].legend(loc="lower left", frameon=False, fontsize=7)
    panel(ax[2], "C", "Capture rate", "captures per time unit γc")
    ax[3].errorbar(algs, col["escape_ratio"], col["escape_ratio_se"], color=ORANGE, **line)
    ax[3].axhline(0, color=MUTED, lw=0.8)
    panel(ax[3], "D", "Escape ratio", "R_esc = 1 − γc / γc,NF")

    ax[4].errorbar(algs, col["nn_distance"], col["nn_distance_sd"], color=BLUE, **line)
    panel(ax[4], "E", "Nearest-neighbour distance", "distance (units of r_d)")
    r = np.corrcoef(col["capture_rate"], col["nn_distance"])[0, 1]
    ax[5].scatter(
        col["nn_distance"],
        col["capture_rate"],
        s=14,
        color=ORANGE,
        edgecolor="white",
        lw=0.5,
        zorder=3,
    )
    for a, x, y in zip(algs, col["nn_distance"], col["capture_rate"], strict=True):
        if a in (0.0, 0.8, 1.6):
            ax[5].annotate(
                f"μ_alg = {a:g}",
                (x, y),
                xytext=(0, 7) if a == 1.6 else (-5, 0),
                textcoords="offset points",
                color=MUTED,
                fontsize=7,
                ha="left" if a == 1.6 else "right",
                va="center",
            )
    panel(ax[5], "F", f"Capture rate vs distance (R = {r:.2f})", "γc")
    ax[5].set_xlabel("nearest-neighbour distance")

    for k in range(5):
        mark_transition(ax[k], transition)
        ax[k].set_xlabel("alignment strength μ_alg")
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIGURES}sweep.{ext}", dpi=200)
    plt.close(fig)


SNAPSHOT_ALGS = (0.2, 0.8, 1.6)


def snapshot(alg):
    # attack of one time step: the first frame is the school when the predator appears
    r = simulate(Params(alg_strength=alg, **(TABLE_S1 | {"attack_time": 0.02})), seed=0)
    return (
        r["prey_pos"][0],
        r["prey_u"][0],
        r["pred_pos"][0],
        r["pred_u"][0],
        r["polarization_start"],
    )


def snapshot_figure():
    with Pool(len(SNAPSHOT_ALGS)) as pool:
        shots = pool.map(snapshot, SNAPSHOT_ALGS)
    fig, ax = plt.subplots(1, 3, figsize=(7, 2.8), layout="constrained")
    names = ("disordered", "transition", "ordered")
    # same scale in the three panels, centred on each school
    half = max(np.nanmax(np.abs(np.vstack([s[0], s[2]]) - np.nanmean(s[0], axis=0))) for s in shots)
    for k, (alg, (pos, u, x, v, pol)) in enumerate(zip(SNAPSHOT_ALGS, shots, strict=True)):
        keep = ~np.isnan(pos[:, 0])
        centre = pos[keep].mean(axis=0)
        ax[k].quiver(*pos[keep].T, *u[keep].T, color=INK, scale=40, width=0.005)
        ax[k].quiver(*x, *v, color=RED, scale=10, width=0.014)
        ax[k].set_xlim(centre[0] - half - 1, centre[0] + half + 1)
        ax[k].set_ylim(centre[1] - half - 1, centre[1] + half + 1)
        ax[k].set_aspect("equal")
        ax[k].set_xticks([])
        ax[k].set_yticks([])
        for s in ax[k].spines.values():
            s.set_visible(False)
        ax[k].set_title(
            f"{'ABC'[k]}  μ_alg = {alg:g}, {names[k]}\nΦ = {pol:.2f}",
            loc="left",
            fontweight="bold",
        )
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIGURES}snapshots.{ext}", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    p = Params(**TABLE_S1)
    data = load(f"{DATA}sweep_control.jsonl")
    rows = summarize(data, p)
    write_csv(rows, f"{DATA}sweep_summary.csv")
    print_numbers(data, rows)
    sweep_figure(rows, data)
    if "--snapshots" in sys.argv:
        snapshot_figure()
