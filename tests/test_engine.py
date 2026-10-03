"""Engine checks against analytic benchmarks (run: .venv/bin/python -m pytest -q or directly)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
import numpy as np
from align.orient import Lattice, point_group
from align.variants import or_rotation, variant_rotations, distinct_axes, twin_rotation, axis_angle_deg
from align import engine as E

def test_groups():
    assert len(point_group('m-3m')) == 24
    assert len(point_group('6/mmm')) == 12
    assert len(point_group('4/mmm')) == 8

def test_random_uniaxial():
    g = E.texture_random()
    axes, frac, per = E.product_axes(g, np.eye(3)[None], np.array([0, 0, 1.0]))
    A = E.alignment(axes, frac)
    assert abs(A - 0.5) < 0.005, A

def cube_variants():
    L = Lattice('cubic', 1.0)
    R = or_rotation(L, L, [0, 0, 1], [0, 0, 1], [1, 0, 0], [1, 0, 0])
    V = variant_rotations(R, 'm-3m', '4/mmm')
    return V

def test_three_variant_selection():
    V = cube_variants()
    assert len(V) == 3
    g = E.texture_random()
    axes, frac, per = E.product_axes(g, V, np.array([0, 0, 1.0]))
    A0 = E.alignment(axes, frac)
    w = E.bias_weights(axes, frac, len(V), per, 'field')
    A1 = E.alignment(axes, w)
    assert abs(A0 - 0.5) < 0.005, A0
    assert abs(A1 - 0.8312) < 0.004, A1          # cubic K1>0 isotropic remanence

def test_fibre_field():
    V = cube_variants()
    g = E.texture_fibre([0, 0, 1])
    axes, frac, per = E.product_axes(g, V, np.array([0, 0, 1.0]))
    assert abs(E.alignment(axes, frac) - 1/3) < 1e-6
    w = E.bias_weights(axes, frac, len(V), per, 'field')
    assert abs(E.alignment(axes, w) - 1.0) < 1e-9


def test_spread_fwhm():
    # easy axis along fibre axis with Gaussian tilt sigma: A = E[cos theta], theta Rayleigh(sigma)
    import numpy as np
    for sig in (5.0, 15.0):
        g = E.texture_fibre([0, 0, 1], n_az=60, spread_deg=sig)
        g = np.concatenate([E.texture_fibre([0, 0, 1], n_az=60, spread_deg=sig, seed=s) for s in range(200)])
        axes, frac, per = E.product_axes(g, np.eye(3)[None], np.array([0, 0, 1.0]))
        A = E.alignment(axes, frac)
        rng = np.random.default_rng(0)
        th = np.hypot(*rng.normal(scale=np.radians(sig), size=(2, 400000)))
        assert abs(A - np.cos(th).mean()) < 0.003, (sig, A, np.cos(th).mean())


if __name__ == '__main__':
    for k, f in list(globals().items()):
        if k.startswith('test_'):
            f(); print('ok', k)
