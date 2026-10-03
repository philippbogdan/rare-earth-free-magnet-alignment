"""The alignment engine: texture source -> transformation (OR variants, bias, twins)
-> optional free-body rotation -> alignment factor A = <|e.z|>.

All functions are deterministic (fixed seeds / grids)."""
import numpy as np

from .orient import random_rotations, rot, fibonacci_hemisphere

N_RANDOM = 40000          # parent orientations for a random texture
N_AZIMUTH = 720           # rotations about a fibre axis
TIE = 1e-6


# ---------------------------------------------------------------- texture sources
def _align_vec(a, b):
    """Rotation taking unit vector a onto unit vector b."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    s, c = np.linalg.norm(v), a @ b
    if s < 1e-12:
        if c > 0:
            return np.eye(3)
        perp = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
        return rot(np.cross(a, perp), 180)
    return rot(v, np.degrees(np.arctan2(s, c)))


def _spread(g, spread_deg, seed):
    """Misorientation spread: each orientation is rotated by a rotation vector whose three
    Cartesian components are independent N(0, spread_deg). The tilt of any crystal axis is
    then a 2-D Gaussian with that standard deviation, so a rocking curve has
    FWHM = 2.3548 * spread_deg. 0 = none."""
    if not spread_deg:
        return g
    rng = np.random.default_rng(seed)
    rv = rng.normal(scale=spread_deg, size=(len(g), 3))
    ang = np.linalg.norm(rv, axis=1)
    R = np.array([rot(v, a) if a > 0 else np.eye(3) for v, a in zip(rv, ang)])
    return R @ g


def texture_random(n=N_RANDOM, seed=1):
    return random_rotations(n, seed)


def texture_fibre(d_crystal, n_az=N_AZIMUTH, spread_deg=0.0, seed=2, axis=(0, 0, 1)):
    """Crystal direction d_crystal (Cartesian, crystal frame) parallel to the sample
    axis, free rotation about it."""
    f = np.asarray(axis, float)
    g0 = _align_vec(np.asarray(d_crystal, float), f)
    psis = np.arange(n_az) * 360.0 / n_az
    g = np.array([rot(f, p) @ g0 for p in psis])
    return _spread(g, spread_deg, seed)


def texture_single(g=None, n=1, spread_deg=0.0, seed=3):
    g = np.eye(3) if g is None else np.asarray(g, float)
    return _spread(np.repeat(g[None], n, 0), spread_deg, seed)


# ---------------------------------------------------------------- transformation
def product_axes(parents, variants, e_H, twins=()):
    """Easy axes (sample frame) of all variants of all parents.

    parents: (N,3,3); variants: (k,3,3) product->parent rotations; e_H: easy axis in the
    product crystal frame; twins: sequence of (T, fraction) twin rotations (product frame).
    Returns axes (N, k*(1+len(twins)), 3) and per-variant base fractions (same shape[:2]).
    """
    ops, frac = [], []
    f_twin = sum(f for _, f in twins)
    for Rk in variants:
        ops.append(Rk @ e_H)
        frac.append(1.0 - f_twin)
        for T, f in twins:
            ops.append(Rk @ T @ e_H)
            frac.append(f)
    ops = np.array(ops)                         # (m,3) in parent frame
    axes = np.einsum('nij,mj->nmi', parents, ops)
    nvar = len(variants)
    frac = np.array(frac) / nvar
    return axes, np.broadcast_to(frac, axes.shape[:2]).copy(), len(ops) // nvar


def _group_select(score, nvar, per):
    """Complete selection: per parent, keep the variant group(s) (a variant with its twins)
    whose untwinned member has the maximum score; equal split on ties."""
    n = score.shape[0]
    s = score.reshape(n, nvar, per)[:, :, 0]
    best = s.max(1, keepdims=True)
    keep = (s >= best - TIE).astype(float)
    keep /= keep.sum(1, keepdims=True)
    return np.repeat(keep, per, axis=1)


def bias_weights(axes, frac, nvar, per, kind='none', direction=(0, 0, 1),
                 strength='complete', sign=+1):
    """Variant weights. kind: 'none' | 'field' | 'stress'.

    field: favours variants whose easy axis is nearest the field direction.
    stress: sign=+1 favours easy axis parallel to the stress axis (e.g. tension with c/a>1,
            or compression with c/a<1); sign=-1 favours easy axis perpendicular to it.
    strength: 'complete' or a float beta: w ~ exp(beta * cos^2) (anisotropy-energy form).
    """
    if kind == 'none':
        return frac
    h = np.asarray(direction, float)
    h /= np.linalg.norm(h)
    c2 = np.einsum('nmi,i->nm', axes, h) ** 2
    score = c2 if (kind == 'field' or sign > 0) else -c2
    if strength == 'complete':
        g = _group_select(score, nvar, per)
        n = axes.shape[0]
        base = frac.reshape(n, nvar, per)
        base = base / base.sum(2, keepdims=True)            # twin split within a group
        return g * base.reshape(n, -1)
    beta = float(strength)
    n = axes.shape[0]
    s = score.reshape(n, nvar, per)[:, :, 0]
    w = np.exp(beta * (s - s.max(1, keepdims=True)))
    w /= w.sum(1, keepdims=True)
    base = frac.reshape(n, nvar, per)
    base = base / base.sum(2, keepdims=True)
    return (np.repeat(w, per, axis=1) * base.reshape(n, -1))


# ---------------------------------------------------------------- alignment
def alignment(axes, weights, z=(0, 0, 1)):
    z = np.asarray(z, float)
    z /= np.linalg.norm(z)
    a = np.abs(np.einsum('...i,i->...', axes, z))
    return float((a * weights).sum() / weights.sum())


def body_alignment(axes, weights, tol=1e-6, n_sub=64):
    """Free-body rotation: each body (first index) turns the principal axis of
    T = sum w e e^T to the field (second-order anisotropy energy). If the largest eigenvalue
    is degenerate the body has no unique axis: the field direction is uniform over the
    degenerate subspace (a circle, or the sphere for an isotropic body)."""
    T = np.einsum('nm,nmi,nmj->nij', weights, axes, axes)
    vals, vecs = np.linalg.eigh(T)
    tot, wsum = 0.0, 0.0
    D_sphere = fibonacci_hemisphere(400)
    t = np.linspace(0, np.pi, n_sub, endpoint=False)
    for k in range(len(T)):
        v, U = vals[k], vecs[k]
        deg = int(np.sum(v >= v[-1] - tol * max(1.0, abs(v[-1]))))
        if deg == 1:
            dirs = U[:, -1][None]
        elif deg == 2:
            a, b = U[:, -1], U[:, -2]
            dirs = np.cos(t)[:, None] * a + np.sin(t)[:, None] * b
        else:
            dirs = D_sphere
        a_k = (np.abs(axes[k] @ dirs.T) * weights[k][:, None]).sum(0).mean()
        tot += a_k
        wsum += weights[k].sum()
    return float(tot / wsum)


def best_axis(axes, weights, n_dir=2000):
    """Largest A over candidate magnet axes (hemisphere grid), and that axis."""
    D = fibonacci_hemisphere(n_dir)
    flat = axes.reshape(-1, 3)
    w = weights.reshape(-1)
    vals = (np.abs(flat @ D.T) * w[:, None]).sum(0) / w.sum()
    i = int(np.argmax(vals))
    return float(vals[i]), D[i]
