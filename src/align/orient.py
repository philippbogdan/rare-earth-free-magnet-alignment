"""Orientation mathematics: rotations, point groups, lattices, sampling.

Convention: an orientation g (3x3) maps crystal Cartesian coordinates to sample
coordinates, v_sample = g @ v_crystal.
"""
import numpy as np


def rot(axis, angle_deg):
    """Rotation matrix about `axis` by `angle_deg` (right-handed)."""
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    t = np.radians(angle_deg)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K


def quat_to_mat(q):
    """Unit quaternions (n,4) [w,x,y,z] -> rotation matrices (n,3,3)."""
    w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    return np.stack([
        np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1),
        np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
        np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1),
    ], -2)


def random_rotations(n, seed=0):
    """n rotations uniformly distributed on SO(3) (Shoemake), fixed seed."""
    rng = np.random.default_rng(seed)
    u1, u2, u3 = rng.random((3, n))
    q = np.stack([np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
                  np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
                  np.sqrt(u1) * np.sin(2 * np.pi * u3),
                  np.sqrt(u1) * np.cos(2 * np.pi * u3)], -1)
    return quat_to_mat(q)


def fibonacci_hemisphere(n):
    """n nearly uniform unit vectors on the upper hemisphere (z >= 0)."""
    i = np.arange(n) + 0.5
    z = 1 - i / n            # z in (0,1)
    phi = np.pi * (1 + 5 ** 0.5) * i
    r = np.sqrt(1 - z * z)
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], -1)


def _unique_mats(mats, tol=1e-6):
    out = []
    for m in mats:
        if not any(np.allclose(m, o, atol=tol) for o in out):
            out.append(m)
    return out


def _closure(gens):
    G = [np.eye(3)]
    changed = True
    while changed:
        changed = False
        for g in list(G):
            for h in gens:
                p = g @ h
                if not any(np.allclose(p, o, atol=1e-8) for o in G):
                    G.append(p)
                    changed = True
    return G


def point_group(name):
    """Proper rotations of the Laue class (orientation symmetry), Cartesian frame.

    Hexagonal frame: a1 along x, c along z. Tetragonal/cubic/orthorhombic: axes along x,y,z.
    """
    z, x, y = [0, 0, 1], [1, 0, 0], [0, 1, 0]
    if name in ('m-3m', 'cubic', '432'):
        gens = [rot(z, 90), rot([1, 1, 1], 120)]
    elif name in ('6/mmm', 'hexagonal', '622'):
        gens = [rot(z, 60), rot(x, 180)]
    elif name in ('-3m', 'trigonal', '32'):
        gens = [rot(z, 120), rot(x, 180)]
    elif name in ('4/mmm', 'tetragonal', '422'):
        gens = [rot(z, 90), rot(x, 180)]
    elif name in ('mmm', 'orthorhombic', '222'):
        gens = [rot(z, 180), rot(x, 180)]
    elif name in ('2/m', 'monoclinic', '2'):
        gens = [rot(y, 180)]
    elif name in ('1', '-1', 'triclinic'):
        gens = []
    else:
        raise ValueError(name)
    G = _closure(gens)
    G = [np.round(g, 12) for g in G]
    return np.array(G)


class Lattice:
    """Direct/reciprocal lattice in a Cartesian frame (hexagonal: a1 || x, c || z)."""

    def __init__(self, system, a, c=None, b=None):
        self.system = system
        c = a if c is None else c
        b = a if b is None else b
        if system == 'hexagonal':
            A = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]]).T
        elif system in ('cubic', 'tetragonal', 'orthorhombic'):
            A = np.diag([a, b, c])
        else:
            raise ValueError(system)
        self.A = A                      # columns = a1, a2, a3
        self.B = np.linalg.inv(A).T     # columns = reciprocal vectors (no 2*pi)

    def direction(self, uvw):
        """Unit Cartesian vector of direction [uvw] or [uvtw] (4-index, hexagonal)."""
        uvw = list(uvw)
        if len(uvw) == 4:
            u, v, t, w = uvw
            uvw = [u - t, v - t, w]
        v = self.A @ np.asarray(uvw, float)
        return v / np.linalg.norm(v)

    def plane_normal(self, hkl):
        """Unit Cartesian normal of plane (hkl) or (hkil) (4-index, hexagonal)."""
        hkl = list(hkl)
        if len(hkl) == 4:
            h, k, i, l = hkl
            hkl = [h, k, l]
        n = self.B @ np.asarray(hkl, float)
        return n / np.linalg.norm(n)
