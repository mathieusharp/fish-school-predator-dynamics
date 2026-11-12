"""Figures, summary table and key numbers for the reports (issue #7).

Usage, from the repository root:
    python3 -m scripts.report_figures [--snapshots]

Reads report/data/sweep_control.jsonl (written by scripts/sweep_alignment.py) and writes
    report/figures/sweep.pdf, .png      Fig 2 of the paper at D = 0.5, six panels
    report/data/sweep_summary.csv       one line per alignment strength
and prints the numbers quoted in report/RESULTS.md. Takes a few seconds.
With --snapshots it also runs three simulations (about 30 s) and writes
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
from scipy.stats import ttest_ind, ttest_rel  # noqa: E402

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
        # Helvetica (its free clone Nimbus Sans), as in the captions of the PNAS template
        "font.family": "sans-serif",
        "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Nimbus Sans",
        "mathtext.it": "Nimbus Sans:italic",
        "mathtext.bf": "Nimbus Sans:bold",
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 3,  # Nimbus Sans is OpenType (CFF): Type 42 embeds it with a wrong font type
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
        iid = np.array(data["iid_start"][a])
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
                "iid": iid.mean(),
                "iid_sd": iid.std(ddof=1),
                "nn_distance": dist.mean(),
                "nn_distance_sd": dist.std(ddof=1),
            }
        )
    return rows


def phase_values(data, key, phase):
    return window(data, key, *PHASES[phase])


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
    print(f"largest IID at mu_alg = {algs[np.argmax(col['iid'])]:g}")
    print(f"largest nn distance at mu_alg = {algs[np.argmax(col['nn_distance'])]:g}")
    print("\nby phase (all runs of the phase pooled):")
    for phase in PHASES:
        rate = phase_values(data, "capture_rate", phase)
        nf = phase_values(data, "capture_rate_nf", phase)
        dist = phase_values(data, "nn_distance_start", phase)
        iid = phase_values(data, "iid_start", phase)
        print(
            f"  {phase:<10} mu_alg {PHASES[phase][0]}-{PHASES[phase][1]}: "
            f"gamma_c {rate.mean():.4f}, gamma_c,NF {nf.mean():.4f}, "
            f"R_esc {1 - rate.mean() / nf.mean():.3f}, IID {iid.mean():.2f}, "
            f"nn distance {dist.mean():.3f}"
        )
    trans = phase_values(data, "capture_rate", "transition")
    for phase in ("disordered", "ordered"):
        other = phase_values(data, "capture_rate", phase)
        t = ttest_ind(trans, other, equal_var=False)
        print(
            f"  transition vs {phase}: {100 * (1 - trans.mean() / other.mean()):.1f} % fewer "
            f"captures, Welch t = {t.statistic:.2f}, p = {t.pvalue:.1g}"
        )
    # the window 0.6-0.8 surrounds our minimum: same test centred on the paper's transition
    late = window(data, "capture_rate", 0.8, 1.0)
    for phase in ("disordered", "ordered"):
        other = phase_values(data, "capture_rate", phase)
        t = ttest_ind(late, other, equal_var=False)
        print(
            f"  window 0.8-1.0 vs {phase}: {100 * (1 - late.mean() / other.mean()):.1f} % fewer, "
            f"p = {t.pvalue:.1g}"
        )
    # control school: is there a dip without fleeing?
    nf_trans = window(data, "capture_rate_nf", 0.6, 0.9)
    print(f"  control school, lowest at mu_alg = {algs[np.argmin(col['capture_rate_nf'])]:g}")
    for phase in ("disordered", "ordered"):
        other = phase_values(data, "capture_rate_nf", phase)
        t = ttest_ind(nf_trans, other, equal_var=False)
        print(
            f"  control school 0.6-0.9 vs {phase}: "
            f"{100 * (1 - nf_trans.mean() / other.mean()):.1f} % fewer, p = {t.pvalue:.1g}"
        )
    rate_all = np.concatenate([data["capture_rate"][a] for a in algs])
    nf_all = np.concatenate([data["capture_rate_nf"][a] for a in algs])
    print(
        f"  all runs: fleeing {rate_all.mean():.4f} vs not fleeing {nf_all.mean():.4f}, "
        f"R_esc {1 - rate_all.mean() / nf_all.mean():.3f}, "
        f"paired t-test p = {ttest_rel(rate_all, nf_all).pvalue:.1g}"
    )
    slope = weighted_slope(algs, col["escape_ratio"], col["escape_ratio_se"])
    print(f"  R_esc against mu_alg: slope {slope[0]:.3f} ± {slope[1]:.3f} (weighted fit)")
    high = algs >= 0.6
    slope = weighted_slope(algs[high], col["escape_ratio"][high], col["escape_ratio_se"][high])
    print(f"  same for mu_alg >= 0.6: slope {slope[0]:.3f} ± {slope[1]:.3f}")

    print("\nlocation uncertainty (bootstrap over runs, 2000 resamples):")
    for name, where in bootstrap_locations(data, algs).items():
        share = ", ".join(f"{a:g}: {100 * f:.0f} %" for a, f in where if f >= 0.01)
        print(f"  {name}: {share}")
    c = np.polyfit(algs, col["capture_rate"], 2)
    print(f"  vertex of a parabola fitted to the mean capture rates: {-c[1] / (2 * c[0]):.2f}")

    print("\nPearson correlation with the capture rate:")
    for key, name in (("iid", "IID"), ("nn_distance", "nn distance")):
        runs = np.concatenate([data[f"{key}_start"][a] for a in algs])
        within = np.mean(
            [np.corrcoef(data["capture_rate"][a], data[f"{key}_start"][a])[0, 1] for a in algs]
        )
        print(
            f"  {name}: R = {np.corrcoef(col['capture_rate'], col[key])[0, 1]:.2f} over the "
            f"{len(algs)} means, {np.corrcoef(rate_all, runs)[0, 1]:.2f} over the "
            f"{len(rate_all)} runs, {within:.2f} on average within one alignment strength"
        )


def window(data, key, lo, hi):
    return np.concatenate([v for a, v in data[key].items() if lo - 1e-9 <= a <= hi + 1e-9])


def weighted_slope(x, y, se):
    """Slope of a straight line fitted with weights 1/se, and its standard error."""
    coef, cov = np.polyfit(x, y, 1, w=1 / se, cov="unscaled")
    return coef[0], np.sqrt(cov[0, 0])


def bootstrap_locations(data, algs, n_boot=2000):
    """Share of resamples where the maximum of chi and the minimum of the capture rate fall
    at each alignment strength."""
    rng = np.random.default_rng(0)
    pol = [np.array(data["polarization_start"][a]) for a in algs]
    rate = [np.array(data["capture_rate"][a]) for a in algs]
    chi_max, rate_min = [], []
    for _ in range(n_boot):
        pick = [rng.integers(len(r), size=len(r)) for r in rate]
        chi_max.append(algs[np.argmax([p[i].var() for p, i in zip(pol, pick, strict=True)])])
        rate_min.append(algs[np.argmin([r[i].mean() for r, i in zip(rate, pick, strict=True)])])
    return {
        name: [(a, np.mean(np.array(x) == a)) for a in algs]
        for name, x in (("maximum of chi", chi_max), ("minimum of capture rate", rate_min))
    }


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
    fig, ax = plt.subplots(2, 3, figsize=(7.2, 4.4), layout="constrained")
    ax = ax.ravel()

    ax[0].errorbar(algs, col["polarization"], col["polarization_sd"], color=BLUE, **line)
    panel(ax[0], "A", "Polarization", r"$\Phi$ (± s.d.)")
    ax[1].plot(algs, col["susceptibility"], color=BLUE, marker="o", ms=3.5, lw=1.5)
    panel(ax[1], "B", "Susceptibility", r"$\chi$")
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
        label=r"not fleeing (control), $\gamma_{c,\mathrm{NF}}$",
        **(line | {"marker": "s"}),
    )
    ax[2].errorbar(
        algs,
        col["capture_rate"],
        col["capture_rate_se"],
        color=ORANGE,
        label=r"fleeing, $\gamma_c$",
        **line,
    )
    ax[2].set_ylim(0.10, None)  # room for the legend below the curves
    ax[2].legend(loc="lower left", frameon=False, fontsize=7)
    panel(ax[2], "C", "Capture rate", r"capture rate (± s.e.)")
    ax[3].errorbar(algs, col["escape_ratio"], col["escape_ratio_se"], color=ORANGE, **line)
    ax[3].axhline(0, color=MUTED, lw=0.8)
    panel(ax[3], "D", "Escape ratio", r"$R_\mathrm{esc} = 1 - \gamma_c / \gamma_{c,\mathrm{NF}}$")

    ax[4].errorbar(algs, col["iid"], col["iid_sd"], color=BLUE, **line)
    panel(ax[4], "E", "Inter-individual distance", r"IID (units of $r_d$)")
    ax[5].scatter(
        col["iid"],
        col["capture_rate"],
        s=14,
        color=ORANGE,
        edgecolor="white",
        lw=0.5,
        zorder=3,
    )
    for a, x, y in zip(algs, col["iid"], col["capture_rate"], strict=True):
        if a in (0.0, 0.8, 1.6):
            ax[5].annotate(
                rf"$\mu_\mathrm{{alg}}$ = {a:g}",
                (x, y),
                xytext={1.6: (0, 7), 0.0: (-12, 12)}.get(a, (-5, 0)),
                textcoords="offset points",
                color=MUTED,
                fontsize=7,
                ha="left" if a == 1.6 else "right",
                va="center",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5) if a == 0.0 else None,
            )
    panel(ax[5], "F", "Capture rate vs IID", "capture rate")
    ax[5].set_xlabel(r"IID (units of $r_d$)")

    for k in range(5):
        mark_transition(ax[k], transition)
        ax[k].set_xlabel(r"alignment strength $\mu_\mathrm{alg}$")
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
