"""Tests of the predprey package. Run from the repository root: python3 -m pytest"""

import numpy as np

from predprey.geometry import (
    cluster_labels,
    delaunay_edges,
    largest_cluster,
    norm,
    perpendicular,
    ray_to_alpha_shape,
    turn,
    turn_towards,
    unique_edges,
)
from predprey.params import Params
from predprey.simulation import alignment_strengths, init_school, prey_step, simulate


def test_delaunay_edges_square_with_centre():
    points = np.array([[0, 0], [2, 0], [2, 2], [0, 2], [1, 1]], dtype=float)
    edges = delaunay_edges(points)
    assert (edges[:, 0] < edges[:, 1]).all()
    assert len(np.unique(edges, axis=0)) == len(edges)
    # the centre is a Voronoi neighbour of the four corners
    assert {tuple(e) for e in edges if 4 in e} == {(0, 4), (1, 4), (2, 4), (3, 4)}


def test_delaunay_edges_3d_cube_with_centre():
    corners = np.array([[x, y, z] for x in (0, 2) for y in (0, 2) for z in (0, 2)], dtype=float)
    points = np.vstack([corners, [1, 1, 1]])
    edges = delaunay_edges(points)
    assert (edges[:, 0] < edges[:, 1]).all()
    # the centre is a Voronoi neighbour of the eight corners
    assert {tuple(e) for e in edges if 8 in e} == {(i, 8) for i in range(8)}


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


def test_turn_towards_is_turn_in_2d_and_works_in_3d():
    u = np.array([[1.0, 0.0], [0.6, 0.8]])
    angle = np.array([0.3, -1.2])
    assert np.allclose(turn_towards(u, angle[:, None] * perpendicular(u)), turn(u, angle))
    u3 = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    w = np.array([[0.0, 0.0, np.pi / 2], [0.0, 0.0, 0.0]])  # quarter turn up, no turn
    assert np.allclose(turn_towards(u3, w), [[0, 0, 1], [0, 0, 1]])


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


def test_init_school_3d():
    p = Params(dim=3, N=400)
    pos, u = init_school(p, np.random.default_rng(0))
    assert pos.shape == (400, 3)
    assert np.allclose(u, u[0])
    assert np.isclose(norm(u[0]), 1)
    assert len(largest_cluster(pos, p.cludist)) == 400  # one compact school


def test_prey_step_3d():
    p = Params(dim=3, cpp_neighbour_count=False, Dphi=0)
    # two prey closer than rep_range, side by side: they move apart in 3D too
    pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.0]])
    u = np.tile([1.0, 0.0, 0.0], (2, 1))
    new_pos, new_u = prey_step(pos, u, np.full(2, 2.0), p, np.random.default_rng(0))
    assert new_pos.shape == (2, 3)
    assert norm(new_pos[1] - new_pos[0]) > 0.5
    assert np.allclose(norm(new_u), 1)
    # with noise, headings stay unit vectors
    p.Dphi = 0.5
    _, new_u = prey_step(pos, u, np.full(2, 2.0), p, np.random.default_rng(0))
    assert np.allclose(norm(new_u), 1)


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
