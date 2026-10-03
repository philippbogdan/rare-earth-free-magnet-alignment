"""Engine checks against values this work does not produce (sourced records).

Random assemblies of non-interacting grains, each relaxing to the easy direction nearest
the field axis (Stoner-Wohlfarth remanence): the engine's "complete selection among the
easy directions of one crystal" on a random texture must reproduce the published values.
"""
import numpy as np

from . import engine as E
from .orient import point_group

REF_ADAMS = 'uni_phys_random_remanence_table_adams2023'
REF_GESHEV = 'uni_phys_cubic_random_remanence_geshev1998'
REF_SW = 'uni_phys_sw1948_random_uniaxial_remanence'


def _axes_from_group(e, group='m-3m'):
    out = []
    for g in point_group(group):
        v = g @ e
        if not any(abs(abs(v @ o) - 1) < 1e-9 for o in out):
            out.append(v)
    return np.array(out)


def _cone_axes(theta_deg, n=360):
    t, p = np.radians(theta_deg), np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t) * np.ones_like(p)], 1)


def nearest_easy_remanence(easy_dirs, n=40000, seed=7):
    """<max_k |e_k . z|> over random crystal orientations."""
    g = E.texture_random(n, seed)
    ax = np.einsum('nij,mj->nmi', g, easy_dirs)
    return float(np.abs(ax[..., 2]).max(1).mean())


def cases():
    return [
        ('random uniaxial (easy axis)', nearest_easy_remanence(np.array([[0, 0, 1.0]])), 0.5,
         [REF_SW, REF_ADAMS]),
        ('random easy plane', nearest_easy_remanence(_cone_axes(90.0)), 0.785, [REF_ADAMS]),
        ('random easy cone (45 deg)', nearest_easy_remanence(_cone_axes(45.0)), 0.909, [REF_ADAMS]),
        ('random cubic, easy <100>', nearest_easy_remanence(_axes_from_group(np.array([0, 0, 1.0]))),
         0.831, [REF_GESHEV, REF_ADAMS]),
        ('random cubic, easy <111>', nearest_easy_remanence(
            _axes_from_group(np.array([1, 1, 1.0]) / np.sqrt(3))), 0.866, [REF_GESHEV, REF_ADAMS]),
        ('random cubic, easy <110>', nearest_easy_remanence(
            _axes_from_group(np.array([1, 1, 0.0]) / np.sqrt(2))), 0.912, [REF_ADAMS]),
    ]
