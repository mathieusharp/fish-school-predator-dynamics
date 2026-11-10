"""Tests of the predprey package. Run from the repository root: python3 -m pytest"""

from dataclasses import replace

import numpy as np

from predprey.geometry import (
    cluster_labels,
    delaunay_edges,
    largest_cluster,
    norm,
    ray_to_alpha_shape,
    turn,
    unique_edges,
)
from predprey.params import Params
from predprey.simulation import (
    alignment_strengths,
    attack,
    iid,
    init_school,
    nn_distance,
    prey_step,
    simulate,
)


def test_delaunay_edges_square_with_centre():
    points = np.array([[0, 0], [2, 0], [2, 2], [0, 2], [1, 1]], dtype=float)
    edges = delaunay_edges(points)
    assert (edges[:, 0] < edges[:, 1]).all()
    assert len(np.unique(edges, axis=0)) == len(edges)
    # the centre is a Voronoi neighbour of the four corners
    assert {tuple(e) for e in edges if 4 in e} == {(0, 4), (1, 4), (2, 4), (3, 4)}


def test_unique_edges_matches_numpy():
    rng = np.random.default_rng(0)
    edges = np.sort(rng.integers(0, 50, (500, 2)), axis=1)
    expected, first = np.unique(edges, axis=0, return_index=True)
    assert np.array_equal(unique_edges(edges, 50), first)
    assert np.array_equal(edges[unique_edges(edges, 50)], expected)


def test_turn_rotates_and_keeps_unit_length():
    u = np.array([[1.0, 0.0], [0.0, 1.0]])
    rotated = turn(u, [np.pi / 2, np.pi])
    assert np.allclose(rotated, [[0, 1], [0, -1]])
    assert np.allclose(norm(rotated), 1)


def test_clusters_two_groups():
    rng = np.random.default_rng(0)
    big = rng.random((30, 2))
    small = rng.random((10, 2)) + 100
    points = np.vstack([big, small])
    labels = cluster_labels(points, cludist=6)
    assert len(np.unique(labels)) == 2
    assert set(largest_cluster(points, cludist=6)) == set(range(30))


def test_ray_to_alpha_shape():
    # a square block of prey from x = 10 to x = 14
    xs, ys = np.meshgrid(np.arange(10, 15), np.arange(-2, 3))
    points = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    origin = np.array([0.0, 0.0])
    assert np.isclose(ray_to_alpha_shape(points, origin, np.array([1.0, 0.0]), 3), 10)
    assert ray_to_alpha_shape(points, origin, np.array([-1.0, 0.0]), 3) == np.inf


def test_init_school():
    p = Params(N=100)
    pos, u = init_school(p, np.random.default_rng(0))
    assert pos.shape == (100, 2)
    assert np.allclose(u, u[0])  # everyone heads the same way
    assert np.isclose(norm(u[0]), 1)


def test_alignment_strengths():
    p = Params(N=1000, alg_strength=2, alg_sigma=0.1)
    alg = alignment_strengths(p, np.random.default_rng(0))
    assert (alg >= 0).all()
    assert np.allclose(alg, np.round(alg, 1))
    assert abs(alg.mean() - 2) < 0.02


def test_close_prey_repel():
    p = Params(Dphi=0)  # no noise
    # Two prey closer than rep_range. As in the C++ code, only one end of each
    # Voronoi edge counts it, so only one of them turns away: distance still grows.
    pos = np.array([[0.0, 0.0], [0.0, 0.5]])
    u = np.tile([1.0, 0.0], (2, 1))
    new_pos, new_u = prey_step(pos, u, np.full(2, 2.0), p, np.random.default_rng(0))
    assert norm(new_pos[1] - new_pos[0]) > 0.5
    assert np.allclose(norm(new_u), 1)


def test_prey_flee_from_predator():
    p = Params(Dphi=0)
    pos = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.5]])
    u = np.tile([1.0, 0.0], (3, 1))
    pred = np.array([-2.0, -0.5])  # behind and slightly below the school
    _, new_u = prey_step(pos, u, np.zeros(3), p, np.random.default_rng(0), pred)
    assert new_u[0, 1] > 0  # the closest prey turns away from the predator


def test_simulate_is_reproducible():
    p = Params(N=50, pred_time=2, attack_time=2)
    a = simulate(p, seed=3)
    b = simulate(p, seed=3)
    assert np.allclose(a["pred_pos"], b["pred_pos"])
    assert np.array_equal(a["kills"], b["kills"])
    assert a["prey_pos"].shape == (10, 50, 2)  # 2 time units, one frame every 0.2


def test_control_school_leaves_the_attacked_school_unchanged():
    p = Params(N=60, pred_time=5, attack_time=10, kill_rate=3.0)
    a = simulate(p, seed=4)
    b = simulate(replace(p, control=True), seed=4)
    assert "kills_nf" not in a and "kills_nf" in b
    assert np.array_equal(a["kills"], b["kills"])
    assert np.array_equal(a["prey_pos"], b["prey_pos"], equal_nan=True)
    assert np.array_equal(a["pred_pos"], b["pred_pos"])


def test_control_school_starts_when_the_predator_appears_and_does_not_flee():
    p = Params(N=60, pred_time=5, attack_time=10, kill_rate=3.0, control=True)
    r = simulate(p, seed=4)
    # first frame = prey (NaN outside the attacked cluster) and predator when it appears
    pos, u = r["prey_pos"][0].copy(), r["prey_u"][0].copy()
    cluster = np.flatnonzero(~np.isnan(pos[:, 0]))
    twin_rng = np.random.default_rng(np.random.SeedSequence(4).spawn(1)[0])
    kills, frames = attack(
        pos,
        u,
        r["alg"],
        cluster,
        r["pred_pos"][0],
        r["pred_u"][0],
        replace(p, flee_strength=0.0),
        twin_rng,
        record=False,
    )
    assert np.array_equal(kills, r["kills_nf"])
    assert frames["t"] == []


def test_nn_distance():
    grid = 2.0 * np.array([[i, j] for i in range(3) for j in range(3)], dtype=float)
    assert np.isclose(nn_distance(grid), 2.0)
    r = simulate(Params(N=50, pred_time=2, attack_time=1), seed=0)
    assert 0 < r["nn_distance_start"] < 3


def test_iid():
    # three prey on a line, 1 and 2 apart: pairs at 1, 2 and 3
    assert np.isclose(iid(np.array([[0.0, 0.0], [1.0, 0.0], [3.0, 0.0]])), 2.0)
    r = simulate(Params(N=50, pred_time=2, attack_time=1), seed=0)
    assert r["iid_start"] > r["nn_distance_start"]
