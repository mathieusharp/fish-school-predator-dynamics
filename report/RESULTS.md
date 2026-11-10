# Results for the reports

Shared material for writing the reports (issue #7): figures, numbers to quote, setup for the Methods, comparison with the paper and points for the Discussion. Everything below comes from `report/data/` and can be regenerated in a few seconds (see [Reproduce](#reproduce)).

- `figures/`: figures for Overleaf. Upload the PDF version (vector); the PNG is for quick viewing.
- `data/sweep_control.jsonl`: the 680 runs of the alignment sweep, one JSON line per run.
- `data/sweep_summary.csv`: one line per alignment strength (means, errors, χ, R_esc, distance).
- `data/versions.txt`: Python packages and commit used for the sweep.

Paper: Klamser & Romanczuk (2021), PLoS Comput Biol 17(3): e1008832, "the paper" below. "S1 Text" and "Table S1" are its supporting information.

## Status

| Result | Issue | Where | Status |
|---|---|---|---|
| Python model checked against the original C++ code | #1 | Mathieu's computer (`results/cpp.jsonl`, `results/py.jsonl`, `results/comparison.png`) | **to add here** (see [To do](#to-do)) |
| Capture rate minimal at the order transition (Fig 2C) | #5 | `figures/sweep.pdf` C | ready |
| Control school and escape ratio (Fig 2D) | #4 | `figures/sweep.pdf` C, D | ready |
| School structure vs capture rate (inset of Fig 2C) | #4 | `figures/sweep.pdf` E, F | ready, but **not the paper's measure** (see [Comparison](#comparison-with-the-paper)) |
| Snapshots of the three phases | | `figures/snapshots.pdf` | ready |
| Capture rate level of the C++ code | #5 | numbers below (Mathieu, comment in #5) | ready |
| Animation | #6 | `docs/run.gif`, `python3 -m scripts.animate` | ready |
| Phase diagram in the (μ_alg, D) plane (Fig 2A) | | | not done, not needed for report 1 |

## Setup (for the Methods)

Parameters of Table S1 of the paper, with their name in `predprey/params.py`:

| Parameter | Paper | Value | Code |
|---|---|---|---|
| angular noise | D | 0.5 | `Dphi` |
| alignment strength | μ_alg | 0 to 1.6, step 0.1 (swept) | `alg_strength` |
| distance strength | μ_d | 2 | `rep_strength` |
| distance steepness | m_d | 2 | `rep_steepness = -2` (sign convention of the C++ code) |
| preferred distance | r_d | 1 | `rep_range` |
| prey speed | v_0 | 1 | `speed0` |
| flee strength | μ_flee | 4 (0 for the control school) | `flee_strength` |
| predator speed | v_p | 2 | `pred_speed` |
| pursuit strength | μ_p | 2 | `pred_strength` |
| attack rate | γ_a | 1/3 | `kill_rate` |
| catch radius | r_catch | 3 | `kill_range` |
| number of prey | N | 400 | `N` |
| time step | dt | 0.02 | `dt` |
| equilibration time | T_eq | 200 | `pred_time` |
| attack time | T_simu | 120 | `attack_time` |

Values of the C++ code that are not in Table S1: spread of the individual alignment strengths 0.075 (`alg_sigma`, C++ `pava_sig`), predator placed 4.5 away from the school (`pred_radius`) at a random angle centred on the back of the school, ± 3 rad, so from almost any direction (`pred_angle = π`, `pred_angle_noise = 3`, see the comment in `params.py`), two prey closer than 6 are in the same cluster (`cludist`), equilibration restarted up to 3 times if the largest cluster has less than 90 % of the prey.

Protocol:
- 17 alignment strengths × 40 runs = 680 runs. Run i of every alignment strength uses seed i.
- Each run: equilibration for T_eq = 200 without predator, then the predator attacks the largest cluster for T_simu = 120. Here the largest cluster always held all 400 prey.
- Measured when the predator appears, on the attacked cluster: polarization Φ = |mean heading| and mean nearest-neighbour distance.
- Capture rate γc = captures / T_simu, as in the paper ("number of prey captured per time unit").
- Control school (as in the paper): a copy of the school when the predator appears, with the same predator, attacked with μ_flee = 0. Its prey still have the predator as a Voronoi neighbour; they just do not flee. It has its own random stream, so the attacked school is identical with or without it (the C++ code shares one stream between both).
- Susceptibility χ = βN (⟨Φ²⟩ − ⟨Φ⟩²) over the 40 runs, β = 1/(D v_0) = 2 (Eq S44d of S1 Text). The transition is placed at its maximum.
- Escape ratio R_esc = 1 − γc / γc,NF (paper, Methods), from the mean capture rates at each alignment strength. Error: bootstrap over runs (2000 resamples), each run kept with its control.
- Errors: γc and γc,NF ± standard error over runs; Φ and distance ± standard deviation over runs.
- Run time: about 3 h on 8 cores.

## Results (numbers to quote)

From `python3 -m scripts.report_figures`:

- **Transition**: χ peaks at μ_alg = 0.8. Φ goes from 0.05 (μ_alg = 0) to 0.83 (μ_alg = 1.6), and from 0.27 to 0.50 between 0.8 and 0.9.
- **Capture rate**: lowest at μ_alg = 0.8 (γc = 0.143 ± 0.006), the maximum of χ.

| Phase | μ_alg | γc | γc,NF | R_esc | nearest-neighbour distance |
|---|---|---|---|---|---|
| disordered | 0 – 0.2 | 0.172 | 0.181 | 0.05 | 0.73 |
| transition | 0.6 – 0.8 | **0.146** | 0.173 | 0.16 | 0.78 |
| ordered | 1.4 – 1.6 | 0.171 | 0.194 | 0.12 | 0.67 |

(120 runs per phase, pooled.)

- **15 % fewer captures at the transition** than in the disordered phase (Welch t-test, t = −5.3, p = 3·10⁻⁷) and than in the ordered phase (−14.7 %, t = −4.8, p = 3·10⁻⁶).
- **Fleeing helps**: over all runs, γc = 0.163 against γc,NF = 0.183 for the control school, so R_esc = 0.11 (Welch p = 9·10⁻¹⁹).
- **The control school has no minimum at the transition**: γc,NF stays between 0.17 and 0.19 below μ_alg = 1, then rises to 0.19 – 0.20 from μ_alg = 1.2. The dip of γc comes from fleeing, not from the structure of the school alone.
- **Escape ratio**: close to 0 at μ_alg = 0 – 0.1 (and at 0.5), 0.10 – 0.18 from μ_alg = 0.6 to 1.5, 0.07 at 1.6. Its highest mean is at μ_alg = 1.2 (0.18 ± 0.04), but the points from 0.7 to 1.5 are within error of each other: no clear peak.
- **Structure**: the nearest-neighbour distance is largest at μ_alg = 0.7 – 0.8 (0.79, against 0.72 in disorder and 0.66 in order): the school is loosest near the transition. Pearson correlation with γc: R = −0.70 over the 17 means, but only −0.11 over the 680 runs.

### Figure captions (drafts)

**`figures/snapshots.pdf`**: School of 400 prey when the predator (red) appears, after an equilibration of 200 time units, in the disordered phase (A, μ_alg = 0.2), at the transition (B, μ_alg = 0.8) and in the ordered phase (C, μ_alg = 1.6). Arrows give the heading of each prey; Φ is the polarization. Same seed and same scale in the three panels.

**`figures/sweep.pdf`**: Order transition and predation against the alignment strength μ_alg (N = 400, D = 0.5, parameters of Table S1, 40 runs per point). (A) Polarization Φ when the predator appears (± s.d.). (B) Susceptibility χ. (C) Capture rate of the fleeing school γc and of its non-fleeing control γc,NF (± s.e.). (D) Escape ratio R_esc = 1 − γc/γc,NF (± bootstrap s.e.). (E) Mean nearest-neighbour distance when the predator appears (± s.d.). (F) Capture rate against nearest-neighbour distance, one point per alignment strength. Dash-dotted line: maximum of χ (μ_alg = 0.8); dotted line: transition of the paper (μ_alg ≈ 0.9).

## Comparison with the paper

| | Paper (Fig 2, D = 0.5) | Ours | |
|---|---|---|---|
| Order transition | μ_alg ≈ 0.9 (given in the text for the evolution runs, at D = 0.5) | χ max at 0.8 (grid step 0.1) | agrees within the grid step |
| Minimum of γc at the transition | "distinct minimum at the critical point" | at 0.8, 15 % below both phases | **reproduced** |
| Level of γc | about 0.08 – 0.105 (read off Fig 2C, issue #5) | 0.14 – 0.18 | **1.6 × higher**, see below |
| R_esc (Fig 2D) | "no peak at the transition but a continuous increase with alignment strength", maximum "in the strongly ordered region" | rises from about 0 to 0.10 – 0.18, no clear peak, no further increase above μ_alg ≈ 0.7 | no peak at the transition, as in the paper; the continued increase is not seen up to 1.6 |
| Structure vs γc (inset Fig 2C) | inter-individual distance (IID, mean distance over **all pairs**), R = −0.69 | nearest-neighbour distance, R = −0.70 over the 17 means | **not the same measure**: do not present as reproduced |
| Information transfer (Fig 2B) | correlation of velocity fluctuations | not computed | |

## Points for the Discussion

- **Level of the capture rate.** Our γc is about 1.6 times that of Fig 2C. Mathieu ran the original C++ code with the same parameters (comment in #5): C++ 0.169 ± 0.016 (10 runs), 0.157 ± 0.007 (20 runs), 0.161 ± 0.015 (10 runs) at μ_alg = 0.1, 0.8 and 1.5, against 0.172, 0.146 and 0.171 in Python. So the gap is between the published code and the paper, not in our reimplementation. The paper does not say how γc is counted (the attack rate "also accounts for handling time"; caught prey removed or not).
- **Inset of Fig 2C.** The paper's IID is the distance between prey pairs averaged over all pairs; we measured the nearest-neighbour distance. The paper does not say over what its R = −0.69 is computed (means per parameter point, possibly over the whole (μ_alg, D) plane since the insets are maps, or single runs). Our correlation is strong over the means but weak over single runs.
- **One slice of the phase diagram.** The paper's Fig 2 covers the (μ_alg, D) plane; we only test D = 0.5.
- **Neighbour counter (#2).** We follow the C++ code: the social forces are divided by a counter increased for only one end of each Voronoi edge, whereas Eq 2 of the paper averages over all neighbours. The option `cpp_neighbour_count=False` gives the paper's version; its effect on the results has not been measured.
- **When measures are taken.** Φ and the distance are measured once, when the predator appears; the paper does not say when it measures them.
- **Control school.** Same starting state and predator as in the C++ code, but its own random stream (see Setup).

## Reproduce

From the repository root, with the package installed (README):

```bash
python3 -m scripts.sweep_alignment report/data/sweep_control.jsonl 0:1.6:0.1 40   # about 3 h on 8 cores
python3 -m scripts.report_figures --snapshots   # figures, sweep_summary.csv and the numbers above (about 30 s)
```

The sweep in `data/` was made at commit c3addae (`data/versions.txt`). Without `--snapshots`, `report_figures` only redraws the sweep figure (a few seconds).

## To do

- **Mathieu (#1):** add the comparison with the C++ code: `results/cpp.jsonl`, `results/py.jsonl` in `report/data/`, `results/comparison.png` in `report/figures/` (PDF if possible), and its main numbers in this file (transition of the C++ code, Welch tests). This is our validation of the reimplementation.
- **Optional:** the paper's IID (mean distance over all pairs) when the predator appears, to compare properly with the inset of Fig 2C. Needs a new sweep (about 1.5 h if only the equilibration is run).
- The PDF and LaTeX source of each report go in `report/` too (Moodle instructions).
