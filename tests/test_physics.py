import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
import numpy as np
from align import physics as P

def test_ms_ratio_limits():
    assert abs(P.ms_ratio(0.0, 600) - 1) < 1e-12
    assert P.ms_ratio(600.0, 600) == 0
    assert P.ms_ratio(700.0, 600) == 0
    r = [float(P.ms_ratio(t, 600)) for t in np.linspace(0, 599, 50)]
    assert all(a >= b for a, b in zip(r, r[1:]))

def test_sw_bias():
    K, Ms = 1e6, 6e5
    HK = 2 * K / (P.MU0 * Ms)
    assert abs(P.ferro_bias_density(K, Ms, 2 * HK) - K) < 1e-6
    # weak field: K h (2-h) -> 2 K h = mu0 Ms H
    H = 0.001 * HK
    assert abs(P.ferro_bias_density(K, Ms, H) / (P.MU0 * Ms * H) - 1) < 1e-3

def test_selection_fraction():
    assert abs(P.selection_fraction(0.0, 3) - 1/3) < 1e-12
    assert P.selection_fraction(10.0, 3) > 0.999

if __name__ == '__main__':
    for k, f in list(globals().items()):
        if k.startswith('test_'):
            f(); print('ok', k)
