"""Benchmarks of the rod, cubic and scoring geometry added for the spinodal and cubic phases."""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np

from align.account import cubic_A
from align.cases import rods_A
from align.scoring import tilt_A

Z, X = (0, 0, 1), (1, 0, 0)


def close(a, b, tol=3e-3):
    assert abs(a - b) < tol, (a, b)


def test_rods_nearest_random():
    # rods on the <100> nearest the field, random grains: the cubic K1>0 remanence
    close(rods_A('random', Z, Z, 'near')[0], 0.8312)


def test_rods_unbiased():
    ind, cou = rods_A('random', Z, Z, 'none')
    close(ind, 0.5)                       # independent uniaxial rods
    close(cou, 0.8660)                    # coupled network, <111> easy (cubic K1<0 remanence)
    ind, cou = rods_A('fibre', Z, Z, 'none')
    close(ind, 1 / 3)
    close(cou, 1 / math.sqrt(3))


def test_rods_ties():
    d = np.ones(3)
    ind, cou = rods_A('single', d, d, 'near')
    close(ind, 1 / math.sqrt(3))
    close(cou, 1.0)
    d = np.array([1.0, 1.0, 0.0])
    ind, cou = rods_A('single', d, d, 'near')
    close(ind, 1 / math.sqrt(2))
    close(cou, 1.0)


def test_rods_follow_field():
    t = math.radians(30)
    close(rods_A('fibre', (math.sin(t), 0, math.cos(t)), Z, 'field')[0], math.cos(t))


def test_cubic_multiaxial():
    close(cubic_A('<111>'), 0.8660)
    close(cubic_A('<100>'), 0.8312)
    close(cubic_A('<111>', 'f111'), 1.0)
    close(cubic_A('<111>', 'f100'), 1 / math.sqrt(3))


def test_tilt():
    close(tilt_A(0.0), 1.0)
    # small spread: <cos> ~ 1 - sigma^2 (2-D Gaussian tilt)
    s = math.radians(5.0)
    close(tilt_A(5.0), 1 - s * s, 2e-3)


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f()
            print('ok', name)
