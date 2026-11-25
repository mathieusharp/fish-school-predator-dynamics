"""Tests of the predprey package. Run from the repository root: python3 -m pytest"""

from dataclasses import replace

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
    voronoi_neighbours,
)
from predprey.params import Params
from predprey.simulation import (
    alignment_strengths,
    attack,
    create_predator,
    iid,
    init_school,
    nn_distance,
    prey_step,
    reflect,
    simulate,
)


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


def test_voronoi_neighbours_match_delaunay_edges():
    rng = np.random.default_rng(0)
    for dim in (2, 3):
        for spread in (0.5, 5, 30):  # q inside, at the edge of, and far from the points
            points = rng.normal(0, 4, (200, dim))
            q = rng.normal(0, spread, dim)
            edges = delaunay_edges(np.vstack([points, q]))
            assert np.array_equal(voronoi_neighbours(points, q), edges[edges[:, 1] == 200, 0])
    # q on one of the points, and too few points: falls back to the triangulation
    edges = delaunay_edges(np.vstack([points, points[0]]))
    assert np.array_equal(voronoi_neighbours(points, points[0]), edges[edges[:, 1] == 200, 0])
    assert np.array_equal(voronoi_neighbours(points[:2], np.zeros(3)), [0, 1])


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


def test_ray_to_alpha_shape_3d():
    # a cubic block of prey from x = 10 to x = 14
    xs, ys, zs = np.meshgrid(np.arange(10, 15), np.arange(-2, 3), np.arange(-2, 3))
    points = np.column_stack([xs.ravel(), ys.ravel(), zs.ravel()]).astype(float)
    points += np.random.default_rng(0).normal(0, 1e-3, points.shape)  # avoid degenerate cubes
    origin = np.zeros(3)
    assert np.isclose(ray_to_alpha_shape(points, origin, np.array([1.0, 0, 0]), 3), 10, atol=0.01)
    assert ray_to_alpha_shape(points, origin, np.array([-1.0, 0, 0]), 3) == np.inf


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


def test_simulate_3d():
    p = Params(dim=3, cpp_neighbour_count=False, N=50, pred_time=2, attack_time=2)
    a = simulate(p, seed=3)
    b = simulate(p, seed=3)
    assert np.allclose(a["pred_pos"], b["pred_pos"])
    assert a["prey_pos"].shape == (10, 50, 3)
    assert np.allclose(norm(a["pred_u"]), 1)
    # the predator starts outside the school, about pred_radius away
    assert np.nanmin(norm(a["prey_pos"][0] - a["pred_pos"][0])) > p.pred_radius / 2


def test_reflect_mirrors_below_the_surface():
    surface = (np.array([0.0, 0.0, 1.0]), 1.0)
    pos = np.array([[0.0, 0.0, 1.5], [0.0, 0.0, 0.5]])
    u = np.array([[0.0, 0.6, 0.8], [0.0, 0.6, 0.8]])
    pos2, u2 = reflect(pos, u, surface)
    assert np.allclose(pos2, [[0, 0, 0.5], [0, 0, 0.5]])  # only the first one crossed
    assert np.allclose(u2, [[0, 0.6, -0.8], [0, 0.6, 0.8]])
    x, v = reflect(np.array([0.0, 0.0, 2.0]), np.array([0.0, 0.0, 1.0]), surface)  # one agent
    assert np.allclose(x, [0, 0, 0]) and np.allclose(v, [0, 0, -1])


def test_surface():
    for below in (False, True):
        p = Params(dim=3, N=50, pred_time=3, attack_time=4, surface=True, pred_below=below)
        r = simulate(p, seed=2)
        up, level = r["surface"]
        assert np.isclose(norm(up), 1)
        depth = level - r["prey_pos"] @ up
        assert np.nanmin(depth) >= 0  # no prey above the surface
        pred_depth = level - r["pred_pos"][0] @ up
        school_depth = np.nanmean(depth[0])
        if below:  # starts below the school, heading up
            assert pred_depth > school_depth + p.pred_radius / 2
            assert np.allclose(r["pred_u"][0], up)
        else:  # starts at the depth of the school
            assert abs(pred_depth - school_depth) < p.pred_radius


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


def test_create_predator_attack_angle_and_distance():
    # pred_angle = 0: from the front (heading against the school), pi: from behind
    for dim in (2, 3):
        for angle in (0, np.pi / 2, np.pi):
            p = Params(dim=dim, N=200, pred_angle=angle, pred_angle_noise=0)
            rng = np.random.default_rng(1)
            pos, u = init_school(p, rng)
            x, v = create_predator(pos, np.arange(p.N), u[0], p, rng)
            assert np.isclose(norm(v), 1)
            assert np.isclose(v @ u[0], -np.cos(angle))
            assert (pos.mean(axis=0) - x) @ v > 0  # heading towards the school
            assert abs(norm(pos - x).min() - p.pred_radius) < 1  # pred_radius from its edge


def test_create_predator_with_surface():
    p = Params(dim=3, N=200, pred_angle=np.pi / 2, pred_angle_noise=0)
    rng = np.random.default_rng(1)
    pos, u = init_school(p, rng)
    up = np.array([0.0, 0.0, 1.0])
    heading = np.array([1.0, 0.0, 0.0])
    x, v = create_predator(pos, np.arange(p.N), heading, p, rng, up=up)
    assert np.isclose(v @ up, 0)  # horizontal attack, from the side
    assert np.isclose(v @ heading, 0)
    x, v = create_predator(pos, np.arange(p.N), heading, replace(p, pred_below=True), rng, up=up)
    assert np.allclose(v, up)  # from straight below, heading up
    assert x @ up < pos[:, 2].min()


def one_attack_step(prey, p, seed):
    """Kills after one step of attack() with the predator at the origin heading along x.

    Prey given by their positions; four more prey far behind the predator complete the
    triangulation. The prey only move by speed0 * dt (0.02) before the catch.
    """
    p = replace(p, Dphi=0.0, attack_time=p.dt)
    behind = np.array([[-20.0, 5.0], [-20.0, -5.0], [-25.0, 0.0], [-30.0, 2.0]])
    pos = np.vstack([prey, behind])
    u = np.tile([1.0, 0.0], (len(pos), 1))
    kills, _ = attack(
        pos,
        u,
        np.zeros(len(pos)),
        np.arange(len(pos)),
        np.zeros(2),
        np.array([1.0, 0.0]),
        p,
        np.random.default_rng(seed),
        record=False,
    )
    return kills[:, 1].astype(int)


def test_catch_only_in_front_and_in_reach():
    p = Params(kill_rate=1e6)  # any prey that can be caught is caught
    prey = np.array([[1.0, 0.0], [-1.0, 0.5], [0.0, p.kill_range + 0.5]])
    # in front and in reach / behind / out of reach
    assert list(one_attack_step(prey, p, seed=0)) == [0]


def test_catch_probability_is_shared_between_prey_in_reach():
    # one prey at distance d: probability (kill_range - d) / kill_range * kill_rate * dt
    p = Params()
    p = replace(p, kill_rate=0.75 / p.dt)  # 0.5 for one prey at d = 1
    one = np.array([[1.0, 0.0]])
    two = np.array([[np.cos(0.3), np.sin(0.3)], [np.cos(0.3), -np.sin(0.3)]])  # both at d = 1
    trials = 2000
    rate_one = sum(len(one_attack_step(one, p, s)) for s in range(trials)) / trials
    rate_two = sum(len(one_attack_step(two, p, s)) for s in range(trials)) / (2 * trials)
    assert abs(rate_one - 0.5) < 0.04
    assert abs(rate_two - 0.25) < 0.04  # shared: half the chance each


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
