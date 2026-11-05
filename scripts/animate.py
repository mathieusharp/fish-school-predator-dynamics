"""Animation of one run : gif for the README, mp4 for the video (1280x720, 30 fps).

Prey as arrows, predator in red, caught prey disappear.
Run from the repository root: python3 -m scripts.animate
"""

import matplotlib

matplotlib.use("Agg")  # no window, only files
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter

from predprey.params import Params
from predprey.simulation import simulate

# r = simulate(Params(N=100, pred_time=20, attack_time=20), seed=1)  # small run while testing
r = simulate(Params(), seed=1)  # full run

P, U = r["prey_pos"], r["prey_u"]  # shape (frames, N, 2) and NaN if prey is caught
X, V = r["pred_pos"], r["pred_u"]  # shape (frames, 2)

fig, ax = plt.subplots(figsize=(12.8, 7.2), dpi=100)  # 1280*720
ax.set_aspect("equal")

prey = ax.quiver(
    P[0, :, 0], P[0, :, 1], U[0, :, 0], U[0, :, 1], color="black", scale=60
)  # x,y,u,v for every prey
pred = ax.quiver(X[0, 0], X[0, 1], V[0, 0], V[0, 1], color="red", scale=40)  # x,y,u,v for pred


def update(f):
    """Moves the arrows to frame f and zooms so that all prey and the predator stay in view"""
    prey.set_offsets(P[f])
    prey.set_UVC(U[f, :, 0], U[f, :, 1])
    pred.set_offsets(X[f])
    pred.set_UVC(V[f, 0], V[f, 1])
    x_min = min(np.nanmin(P[f, :, 0]), X[f, 0])
    x_max = max(np.nanmax(P[f, :, 0]), X[f, 0])
    y_min = min(np.nanmin(P[f, :, 1]), X[f, 1])
    y_max = max(np.nanmax(P[f, :, 1]), X[f, 1])
    center = np.array([(x_min + x_max) / 2, (y_min + y_max) / 2])
    half_width = max(25, (x_max - x_min) / 2, (y_max - y_min) / 2 * 16 / 9) + 5
    ax.set_xlim(center[0] - half_width, center[0] + half_width)
    ax.set_ylim(
        center[1] - half_width * 9 / 16, center[1] + half_width * 9 / 16
    )  # for a 16:9 screen
    return prey, pred


anim = FuncAnimation(fig, update, frames=len(P[:, 0, 0]))
anim.save("results/run.gif", writer=PillowWriter(fps=15), dpi=40)
anim.save("results/run.mp4", writer=FFMpegWriter(fps=30))
