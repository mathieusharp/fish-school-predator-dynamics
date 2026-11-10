"""Neighbours, clusters and other geometric helpers.

Positions are arrays of shape (n, 2) or (n, 3); directions are unit vectors, never angles.
"""

from itertools import combinations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import Delaunay, QhullError, cKDTree


def delaunay_edges(points):
    """Pairs (i, j), i < j, of Voronoi neighbours (= edges of the Delaunay triangulation)."""
    n, dim = points.shape
    if n <= dim:  # too few points for a triangulation: all are neighbours
        return np.array([(i, j) for i in range(n) for j in range(i + 1, n)], dtype=int).reshape(
            -1, 2
        )
    try:
        tri = Delaunay(points)
    except QhullError:  # degenerate configuration
        try:
            tri = Delaunay(points, qhull_options="QJ")
        except QhullError:  # all points aligned: neighbours are consecutive along the line
            direction = points[-1] - points[0]
            order = np.argsort(points @ direction)
            return np.sort(np.column_stack([order[:-1], order[1:]]), axis=1)
    s = tri.simplices  # triangles in 2D, tetrahedra in 3D
    edges = np.concatenate([s[:, [a, b]] for a, b in combinations(range(dim + 1), 2)])
    edges.sort(axis=1)
    return edges[unique_edges(edges, n)]


def unique_edges(edges, n):
    """Index of the first occurrence of each edge (i, j), i < j, in lexicographic order.

    Same result as np.unique(edges, axis=0, return_index=True)[1], but much faster:
    each edge is encoded as one integer i * n + j (issue #3).
    """
    return np.unique(edges[:, 0].astype(np.int64) * n + edges[:, 1], return_index=True)[1]


def delaunay_edges_counted(points, rng):
    """Delaunay edges (i, j), i < j, and for each edge the end whose neighbour counter the C++ increases.

    C++: CGAL reports each edge once, from one of its two triangles (the one with the
    smaller handle), and IntCalcPrey only increases the counter of f.vertex(cw(i)).
    Handles are modelled by a random rank per triangle. A prey inside the school then
    has a counter between 1 and (number of neighbours - 1), never 0.
    """
    try:
        tri = Delaunay(points)
    except (QhullError, ValueError):  # fewer than 3 points or degenerate: end chosen at random
        edges = delaunay_edges(points)
        return edges, edges[np.arange(len(edges)), rng.integers(0, 2, len(edges))]
    s, nb = tri.simplices.copy(), tri.neighbors.copy()
    p = points[s]
    ab, ac = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]
    clockwise = ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0] < 0  # z of the cross product
    s[clockwise] = s[clockwise][:, [0, 2, 1]]  # CGAL triangles are counter-clockwise
    nb[clockwise] = nb[clockwise][:, [0, 2, 1]]
    rank = rng.permutation(len(s))

    edges, counted = [], []
    for k in range(3):  # edge opposite to vertex k of each triangle
        vi, vj = s[:, (k + 2) % 3], s[:, (k + 1) % 3]  # cw(k), ccw(k)
        other = nb[:, k]
        # on the hull, the other triangle is CGAL's infinite face: reported from either side
        here = np.where(other < 0, rng.random(len(s)) < 0.5, rank < rank[other])
        edges.append(np.column_stack([vi, vj]))
        counted.append(np.where(here, vi, vj))
    edges, counted = np.concatenate(edges), np.concatenate(counted)
    edges.sort(axis=1)
    first = unique_edges(edges, len(points))  # interior edges appear twice
    return edges[first], counted[first]


def norm(v):
    return np.linalg.norm(v, axis=-1)


def unit(v):
    return v / norm(v)[..., None]


def perpendicular(u):
    """u rotated by +90 degrees."""
    return np.stack([-u[..., 1], u[..., 0]], axis=-1)


def turn(u, angle):
    """Rotates the unit vectors u by angle (in radians)."""
    angle = np.asarray(angle)[..., None]
    return unit(u * np.cos(angle) + perpendicular(u) * np.sin(angle))


def perpendicular_part(v, u):
    """Part of the vectors v perpendicular to the unit vectors u."""
    return v - np.sum(v * u, axis=-1, keepdims=True) * u


def turn_towards(u, w):
    """Rotates the unit vectors u by the angle |w| towards w (w perpendicular to u, any dimension).

    In 2D, with w = angle * perpendicular(u), this is turn(u, angle).
    """
    angle = norm(w)[..., None]
    direction = np.divide(w, angle, out=np.zeros_like(w), where=angle > 0)
    return unit(u * np.cos(angle) + direction * np.sin(angle))


def cluster_labels(points, cludist):
    """Label of the cluster of each point (chains of points closer than cludist)."""
    n = len(points)
    pairs = cKDTree(points).query_pairs(cludist, output_type="ndarray")
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    return connected_components(graph, directed=False)[1]


def largest_cluster(points, cludist):
    """Indices of the points in the largest cluster."""
    labels = cluster_labels(points, cludist)
    return np.flatnonzero(labels == np.bincount(labels).argmax())


def ray_to_alpha_shape(points, origin, u, alpha_radius):
    """Distance from origin, along u, to the boundary of the alpha shape of points (C++: DistP2AlphaShape).

    The C++ code uses the exact alpha shape (CGAL); here it is approximated by the
    Delaunay edges shorter than 2 * alpha_radius. Starting from outside the school,
    the first edge hit is on its boundary.
    Returns inf if the ray misses the school.
    In 3D, the boundary is made of the Delaunay triangles whose edges are all shorter than
    2 * alpha_radius.
    """
    if points.shape[1] == 3:
        return ray_to_alpha_shape_3d(points, origin, u, alpha_radius)
    edges = delaunay_edges(points)
    a, b = points[edges[:, 0]], points[edges[:, 1]]
    keep = norm(b - a) <= 2 * alpha_radius
    a, b = a[keep], b[keep]
    # solve origin + t * u = a + s * (b - a), with t > 0 and 0 <= s <= 1
    ab = b - a
    denom = u[0] * ab[:, 1] - u[1] * ab[:, 0]
    ao = a - origin
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (ao[:, 0] * ab[:, 1] - ao[:, 1] * ab[:, 0]) / denom
        s = (ao[:, 0] * u[1] - ao[:, 1] * u[0]) / denom
    hit = (denom != 0) & (t > 0) & (s >= 0) & (s <= 1)
    return t[hit].min() if hit.any() else np.inf


def ray_to_alpha_shape_3d(points, origin, u, alpha_radius):
    """3D version of ray_to_alpha_shape: first Delaunay triangle hit by the ray (Moller-Trumbore)."""
    s = Delaunay(points).simplices
    faces = np.concatenate([s[:, list(f)] for f in combinations(range(4), 3)])
    a, b, c = points[faces[:, 0]], points[faces[:, 1]], points[faces[:, 2]]
    short = (
        (norm(b - a) <= 2 * alpha_radius)
        & (norm(c - a) <= 2 * alpha_radius)
        & (norm(c - b) <= 2 * alpha_radius)
    )
    a, b, c = a[short], b[short], c[short]
    # solve origin + t * u = a + s * (b - a) + v * (c - a), with t > 0, s, v >= 0, s + v <= 1
    ab, ac = b - a, c - a
    pvec = np.cross(u, ac)
    det = np.sum(ab * pvec, axis=1)
    ao = origin - a
    qvec = np.cross(ao, ab)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.sum(ao * pvec, axis=1) / det
        v = (qvec @ u) / det
        t = np.sum(ac * qvec, axis=1) / det
    hit = (det != 0) & (t > 0) & (s >= 0) & (v >= 0) & (s + v <= 1)
    return t[hit].min() if hit.any() else np.inf
