# Results for the reports

Shared material for writing the reports (issue #7): what we reproduce of the paper, setup for the Methods, numbers to quote, figure captions, comparison with the paper and points for the Discussion. Every number below comes from `report/data/` and is printed again by `python3 -m scripts.report_figures` (see [Reproduce](#reproduce)). An independent audit of the code, the data and this file against the paper was made on 2026-11-10; its findings are included below.

- `figures/`: figures for Overleaf. Upload the PDF version (vector); the PNG is for quick viewing.
- `data/sweep_control.jsonl`: the 680 runs of the alignment sweep, one JSON line per run.
- `data/sweep_summary.csv`: one line per alignment strength (means, errors, χ, R_esc, distances).
- `data/versions.txt`: Python packages and commits used.

Paper: Klamser & Romanczuk (2021), PLoS Comput Biol 17(3): e1008832, "the paper" below. "S1 Text" and "Table S1" are its supporting information.

## What we reproduce

**Wording for the report.** The paper has two parts. (1) Fig 2: for the group, the capture rate is lowest at the order–disorder transition (the "group optimum"). (2) Figs 3–4: under evolution, individuals do not stay there; they evolve far into the ordered phase, so criticality is optimal for the group but not evolutionarily stable. That second part is the paper's main conclusion. We reproduce the first part only, so write "we reproduce the group optimum of Fig 2", not "the main result". The second part must still be presented in the Introduction / literature review.

**Fig 2 of the paper is four colour maps** over the angular noise D (0 to 1) and the alignment strength μ_alg (0 to 1.6), 40 runs per point. We reproduce the **vertical slice at D = 0.5** (the dashed line of Fig 2A), over the same μ_alg range.

| Result of the paper | Status | Needed for report 1 |
|---|---|---|
| Fig 2A, polarization Φ | D = 0.5 slice reproduced (`figures/sweep.pdf` A) | yes, done |
| Fig 2B inset, susceptibility χ | D = 0.5 slice reproduced (B) | yes, done |
| Fig 2B main panel, information transfer (correlation of velocity fluctuations between neighbours) | not done | nice to have (equilibration only, about 1 h of computation) |
| Fig 2C, capture rate γc | D = 0.5 slice reproduced, same shape, level about 1.8 × higher (C) | yes, done |
| Fig 2C inset, inter-individual distance IID and its correlation with γc | D = 0.5 slice reproduced (E, F) | yes, done |
| Fig 2D, escape ratio R_esc, with a non-fleeing control school | D = 0.5 slice, partly: noisy, lower than the paper in the ordered phase (C, D) | yes, done with caveats |
| Fig 2D inset, γc,NF − γc | not plotted, computable from the data (`sweep_summary.csv`) | nice to have |
| Fig 1, model sketch | not made; `figures/snapshots.pdf` shows the school in the three phases | nice to have |
| Figs 3–4, evolution of the alignment strength (ESS far above the transition, self-sorting, effect of the flee strength) | not done | not to reproduce, but to describe in the literature review |
| S1–S6 Figs (forces, fitness gradients, robustness checks) | not done | no |
| Python model checked against the original C++ code (issue #1) | done: full sweep with Table S1 parameters, `data/cpp_sweep.jsonl`, green curves of `figures/sweep.pdf` (see [Validation against the C++ code](#validation-against-the-c-code)) | yes, done |
| Literature review and planned approach (required by the assignment) | not in the repository | yes (Introduction) |

## Setup (for the Methods)

### Model as implemented

We reimplement in Python the published C++ code (https://github.com/PaPeK/PredatorPrey), with its default options, rather than the equations as printed in the paper. Where the two differ, the code is what produced the paper's figures (its transition is where the paper's Fig 2A puts it). **Describe the model as implemented, not by copying the paper's equations**:

- **Prey**: 2D, constant speed v_0, heading turned by social forces, flee force and angular noise (Euler–Maruyama integration, time step dt). Social interactions with the Voronoi neighbours, from a Delaunay triangulation that includes the predator: distance regulation (repulsion below r_d, attraction above) and alignment.
- **Neighbour counter (issue #2)**: the social forces are divided by a counter that is increased for only one of the two prey of each Voronoi edge, so by about half the number of neighbours, whereas Eqs S2–S3 of S1 Text divide by the number of neighbours |N_i|. The social forces are thus about twice those of the printed equations. The option `cpp_neighbour_count=False` gives the paper's version; its effect on the results has not been measured.
- **Individual alignment strengths**: drawn around μ_alg with spread 0.075 (C++ `pava_sig`; same value as the mutation strength σ_m of Table S1, but the paper does not say whether Fig 2 used a spread), then rounded to 0.1 (C++ `discretize = 1`). At μ_alg = 0 a few prey get 0.1 or 0.2.
- **Initial condition**: prey on concentric circles, all with the same heading (C++ IC = 3), then 200 time units of equilibration without predator. Equilibration is restarted (up to 3 times) if the largest cluster (prey closer than 6 are in the same cluster) holds less than 90 % of the prey; here it always held all 400.
- **Predator**: appears 4.5 away from the attacked cluster, at a random angle centred on the back of the school, ± 3 rad, so from almost any direction. It moves at constant speed v_p in two states: it first heads for the centre of the cluster; once prey are within 3 r_catch = 9 in front of it, it hunts, heading for the mean position of its frontal Voronoi neighbours weighted by their catch probability. The pursuit force has a fixed size μ_p, whereas Eq S10 makes it proportional to the catch probabilities.
- **Captures**: a frontal Voronoi neighbour at distance r < r_catch is caught during dt with probability (r_catch − r)/r_catch · γ_a · dt divided by the number of neighbours **within r_catch**; Eq S7 divides by the number of all frontal neighbours. The audit tested the paper's version on 4 runs: 5–8 % fewer expected captures, so this does not explain the difference in level with the paper. Caught prey are removed (the paper does not say).
- **Control school** (as in the paper, Fig 2D): a copy of the school and of the predator when the predator appears, attacked with μ_flee = 0. Its prey still have the predator as a Voronoi neighbour; they just do not flee. It has its own random stream, so the attacked school is identical with or without it (the C++ code shares one stream between both).

### Parameters

Table S1 of the paper, with their name in `predprey/params.py`:

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

### Protocol and measures

- 17 alignment strengths × 40 runs = 680 runs (the paper also uses N_s = 40 runs per point). Run i of every alignment strength uses seed i.
- Measured **once, when the predator appears**, on the attacked cluster: polarization Φ = |mean heading|, inter-individual distance IID (distance between prey averaged over all pairs, as in the paper's Methods) and mean nearest-neighbour distance.
- Capture rate γc = captures / T_simu, as in the paper ("number of prey captured per time unit"); γc,NF for the control school.
- Susceptibility χ = βN (⟨Φ²⟩ − ⟨Φ⟩²) over the 40 runs, with β = 1/(D v_0) = 2 (Eq S44d of S1 Text, derived there for the polarization vector; we use its norm Φ). The transition is placed at its maximum.
- Escape ratio R_esc = 1 − γc / γc,NF (paper, Methods), from the mean capture rates at each alignment strength. Error: bootstrap over runs (2000 resamples), each run kept with its control.
- Errors: γc and γc,NF ± standard error over runs; Φ and distances ± standard deviation over runs.
- Run time: about 3 h on 8 cores for the sweep, plus about 1 h to add the IID (`scripts/add_iid.py`, which runs the equilibration again with the same seeds and checks that it finds the same schools).

## Results (numbers to quote)

- **Transition**: χ peaks at μ_alg = 0.8. Φ goes from 0.05 (μ_alg = 0) to 0.83 (μ_alg = 1.6), and from 0.27 to 0.50 between 0.8 and 0.9.
- **Capture rate**: lowest at μ_alg = 0.8 (γc = 0.143 ± 0.006), the maximum of χ.
- **Uncertainty on these locations** (bootstrap over runs): the maximum of χ falls at 0.8 in 47 % of resamples, 0.9 in 33 %, 1.0 in 20 %; the minimum of γc at 0.8 in 52 %, 0.7 in 24 %, 0.6 in 15 %, 0.9–1.0 in 8 %. A parabola fitted to the mean γc has its vertex at 0.80. So write "both lie within about 0.7–1.0", not "exactly at the same value".

| Phase | μ_alg | γc | γc,NF | R_esc | IID | nearest-neighbour distance |
|---|---|---|---|---|---|---|
| disordered | 0 – 0.2 | 0.172 | 0.181 | 0.05 | 12.4 | 0.73 |
| transition | 0.6 – 0.8 | **0.146** | 0.173 | 0.16 | 14.5 | 0.78 |
| ordered | 1.4 – 1.6 | 0.171 | 0.194 | 0.12 | 10.6 | 0.67 |

(120 runs per phase, pooled.)

- **15 % fewer captures at the transition** than in the disordered phase (Welch t-test, t = −5.3, p = 3·10⁻⁷) and than in the ordered phase (15 % fewer, t = −4.8, p = 3·10⁻⁶). The window 0.6–0.8 was chosen around our minimum; centred on the paper's transition instead (0.8–1.0), the result holds: 13 % fewer than in disorder (p = 4·10⁻⁶) and in order (p = 3·10⁻⁵).
- **Fleeing helps**: over all runs, γc = 0.163 against γc,NF = 0.183 for the control school, R_esc = 0.11 (paired t-test, p = 3·10⁻¹⁸).
- **Part of the dip is there without fleeing.** The control school is lowest at μ_alg = 0.9 (γc,NF = 0.170). Over 0.6–0.9 it is 11 % below the ordered phase (p = 5·10⁻⁵) but only 4 % below the disordered phase (p = 0.1, not significant). So on the ordered side the dip comes mostly from the structure of the school (it is there without fleeing), on the disordered side mostly from fleeing. This agrees with the paper, which attributes the optimum mainly to the spatial structure of the school.
- **Escape ratio**: about 0.01–0.03 at μ_alg = 0 – 0.1, 0.07–0.11 at 0.2–0.4, 0.10–0.18 from 0.6 to 1.5, 0.07 at 1.6. Point errors are 0.04–0.05 for values around 0.1, so the shape is uncertain. A weighted straight-line fit gives a mild increase over 0–1.6 (slope 0.054 ± 0.022) and no trend above 0.6 (slope −0.054 ± 0.045). No peak at the transition.
- **Structure**: the school is most spread out near the transition. The IID rises from 12.3 at μ_alg = 0 to a maximum of 15.3 at 0.9, then falls to 10.2 at 1.6 (ordered schools are the most compact). It also varies most between runs there (s.d. 4.1 at 0.9, against 0.3 at 0). The nearest-neighbour distance is largest at μ_alg = 0.7 – 0.8 (0.79, against 0.72 in disorder and 0.66 in order).
- **Correlation with γc** (Pearson): for the IID, R = −0.80 over the 17 means, −0.12 over the 680 runs, 0.01 on average within one alignment strength. For the nearest-neighbour distance: R = −0.70 over the 17 means, −0.11 over the 680 runs. Both are associations **across alignment strengths only**: at a fixed alignment strength, the correlation is about 0. Do not present them as a run-by-run effect.

### Validation against the C++ code

Same sweep as ours, run with the original code (`scripts/run_cpp.py --table-s1`, 17 alignment strengths × 40 runs, 2026-11-11); compared with `python3 -m scripts.validate_cpp` (Welch test per alignment strength).

- **Polarization**: differs at p < 0.05 at 1 of the 17 alignment strengths (0.9, p = 0.03), as expected by chance. χ peaks at 0.9 for C++ (bootstrap: 0.8–0.9 in 94 % of resamples) and 0.8 for Python.
- **Capture rate**: same mean over all runs, 0.161 (C++) against 0.163 (Python), p = 0.42. 4 of 17 points at p < 0.05 (0.5, 0.6, 0.7, 1.6), none after Bonferroni correction (p < 0.003). Python is lower than C++ at 0.6–0.7 (about 0.02, 2.5 s.e.) with the same polarization there.
- **Dip in the C++ code alone**: lowest at 0.8–1.0 (0.147–0.149). Over 0.8–1.0 it is 13 % below the disordered phase (p = 1·10⁻⁶) and 7 % below the ordered phase (p = 0.02); Python gives 13 % and 13 % over the same window. A low C++ point at 1.6 (0.150) makes the ordered side shallower. The window 0.6–0.8 used above shows no significant dip in C++: the report uses 0.8–1.0.
- So the published code reproduces the group optimum and the level of γc of our model (about 1.8 × the paper's Fig 2C).

### Figure captions (drafts)

**`figures/snapshots.pdf`**: School of 400 prey when the predator (red) appears, after an equilibration of 200 time units, in one run of each phase: disordered (A, μ_alg = 0.2), transition (B, μ_alg = 0.8) and ordered (C, μ_alg = 1.6). Arrows give the heading of each prey; Φ is the polarization of that run (the mean over 40 runs at μ_alg = 0.8 is 0.27). Same seed and same scale in the three panels.

**`figures/sweep.pdf`**: Order transition and predation against the alignment strength μ_alg (N = 400, D = 0.5, parameters of Table S1, 40 runs per point). (A) Polarization Φ when the predator appears (± s.d.). (B) Susceptibility χ. (C) Capture rate of the fleeing school γc and of its non-fleeing control γc,NF (± s.e.; the y-axis does not start at 0). (D) Escape ratio R_esc = 1 − γc/γc,NF (± bootstrap s.e.). (E) Inter-individual distance IID (mean distance over all pairs of prey) when the predator appears (± s.d.). (F) Mean capture rate against mean IID, one point per alignment strength. Dash-dotted line: maximum of χ (μ_alg = 0.8); dotted line: transition given in the paper's text (μ_alg ≈ 0.9).

## Comparison with the paper

Values "read off Fig 2" were read by the audit from the colour maps along D ≈ 0.45–0.53; they are approximate. Check them by eye on the paper's Fig 2 before quoting them.

| | Paper (D = 0.5 slice of Fig 2) | Ours | |
|---|---|---|---|
| Order transition | ≈ 0.9 in the text (evolution section, at D = 0.5); ≈ 0.84 read off Fig 2A; 0.83 in a comment of the C++ code (`RunSingle.py`) | χ max at 0.8 (grid step 0.1, bootstrap 0.8–1.0) | agrees |
| Minimum of γc at the transition | "distinct minimum at the critical point"; read off Fig 2C: 0.098 at μ_alg = 0, 0.086 at 0.8, 0.093–0.095 at 1.4–1.6 | at 0.8 (0.7–1.0 with the bootstrap) | **reproduced** |
| Depth of the dip | about 9–11 % below disorder, 6–8 % below order (read off Fig 2C) | 15 % below both | deeper in ours |
| Level of γc | 0.086–0.098 along the slice | 0.143–0.177 | **about 1.8 × higher** (ratio 1.65–1.88 point by point); our 17 means correlate with the paper's slice at r = 0.82 |
| R_esc (Fig 2D) | "no peak at the transition but a continuous increase with alignment strength", maximum "in the strongly ordered region"; read off Fig 2D: about 0.06–0.07 in disorder, 0.13–0.15 at 0.8–0.9, 0.17–0.20 for μ_alg ≥ 1 | about 0 to 0.10 in disorder, 0.10–0.18 above 0.6, 0.12 in the ordered phase (1.4–1.6) | no peak, as in the paper; **lower than the paper in the ordered phase** (about 3 s.e.) and no clear increase above 0.6 |
| Structure vs γc (inset Fig 2C) | IID, R = −0.69 (written C(γc, IID) in the text), probably over the whole (μ_alg, D) map | IID, R = −0.80 over the 17 means of the slice | **reproduced** in sign and size (more spread out school, fewer captures), but over one slice only |
| Information transfer (Fig 2B) | peaks at the transition | not computed | |

## Points for the Discussion

- **Level of the capture rate.** Our γc is about 1.8 times that of Fig 2C, with the same shape. The full C++ sweep with the Table S1 parameters gives the same mean level as ours (0.161 against 0.163, see [Validation against the C++ code](#validation-against-the-c-code)), so the gap is between the published code and the paper's figure, not in our reimplementation. Possible causes: how γc is counted in the paper (the attack rate "also accounts for handling time"; removal of caught prey not stated); the different normalisation of Eq S7 explains only 5–8 %. A normalisation by 200 instead of 120 time units would give a ratio of 1.67, close to ours, but this is a guess we cannot test.
- **Escape ratio.** Lower than the paper in the ordered phase, and noisy (errors of about 0.05 on values of about 0.1). 40 runs per point are enough for γc but not for the shape of R_esc.
- **Inset of Fig 2C.** The IID anticorrelates with γc across alignment strengths (R = −0.80; paper −0.69), but not within one alignment strength (R ≈ 0): schools are more spread out near the transition and are caught less there, but at a given alignment strength a more spread out school is not caught less. The association may come from a common cause (the alignment strength) rather than from the structure itself; the paper does not say whether its correlation holds within one parameter point. The paper's R = −0.69 is probably computed over its whole (μ_alg, D) map; ours is over one slice, so the values are not strictly comparable.
- **One slice of the phase diagram.** We only test D = 0.5.
- **Model as implemented vs paper equations.** Neighbour counter, target selection, pursuit force, two-state predator (see Setup): we follow the C++ code, which produced the paper's figures. The neighbour counter roughly doubles the social forces compared with the printed equations.
- **One measurement per run.** Φ and the distances are measured once, when the predator appears, after 200 time units from an ordered start. Near the transition the school may not have fully relaxed (critical slowing down, discussed in S1 Text, Sec. V.2). χ from 40 single snapshots has a relative error of about 25 %.
- **Choice of windows.** The transition window 0.6–0.8 was chosen after seeing the data; the result holds for a window centred on the paper's transition (see Results), which is the one used in the report (0.8–1.0) because it also holds for the C++ code.

## Reproduce

From the repository root, with the package installed (README):

```bash
python3 -m scripts.sweep_alignment report/data/sweep_control.jsonl 0:1.6:0.1 40   # about 3 h on 8 cores
python3 -m scripts.report_figures --snapshots   # figures, sweep_summary.csv and the numbers above (about 1 min)
python3 scripts/run_cpp.py local/PredatorPrey report/data/cpp_sweep.jsonl 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6 40 --table-s1   # C++ sweep, about 45 min on 8 cores
python3 -m scripts.validate_cpp     # Welch tests against the C++ code
```

The sweep in `data/` was made at commit c3addae and the IID added by `scripts/add_iid.py` (`data/versions.txt`); a new sweep records the IID directly. Without `--snapshots`, `report_figures` only redraws the sweep figure.

The animation of the README (`docs/run.gif`) uses the default parameters of the C++ code (attack 60, γ_a = 1, μ_alg = 2), not those of Table S1: say so if it is shown.

## To do

- **Optional:** information transfer of Fig 2B; effect of the neighbour counter (#2).
- The PDF and LaTeX source of each report go in `report/` too (Moodle instructions).
