"""Route evaluation: the same rules for every phase.

A phase description (data/curated/phases.yaml; every number is a reference to a sourced
record) is turned into, for each route: verdict, ceiling A_max and its magnet axis, the
chain that produces it, the prerequisites checked, any unpublished quantity that decides
the answer (with the ceiling for each of its outcomes) and the single measurement that
tests the prediction.
"""
import numpy as np

from . import engine as E
from . import physics as P
from .orient import Lattice
from .variants import or_rotation, variant_rotations, twin_rotation

ROUTES = ['FP', 'FS', 'DS', 'AM', 'FA', 'DS+FA', 'AM+FA', 'RX+FA', 'SA', 'DS+SA', 'RX+SA',
          'HD', 'EP', 'TGG']

ROUTE_NAMES = {
    'FP': 'field alignment of free particles, then consolidation',
    'FS': 'solidification in a magnetic field',
    'DS': 'directional solidification / columnar casting / melt spinning, no bias',
    'AM': 'additive manufacturing, no bias',
    'FA': 'formation in a magnetic field from an untextured parent',
    'DS+FA': 'columnar parent + formation in a field',
    'AM+FA': 'printed parent + formation in a field',
    'RX+FA': 'rolled/recrystallised parent + formation in a field',
    'SA': 'formation under stress from an untextured parent',
    'DS+SA': 'columnar parent + formation under stress',
    'RX+SA': 'rolled/recrystallised parent + formation under stress',
    'HD': 'hot deformation of the hard phase',
    'EP': 'epitaxial or template-defined growth (films, wire arrays; not bulk)',
    'TGG': 'templated grain growth during sintering',
}

# Best-case sheet textures after rolling/drawing + recrystallisation, by parent lattice
# (sample axes: z = ND, x = RD). Cube {100}<001> for fcc, Goss {110}<001> for bcc, basal
# for hexagonal. Sources: data/curated/records.yaml (texture references).
RX_TEXTURE = {'fcc': ([0, 0, 1], [1, 0, 0]), 'bcc': ([1, 1, 0], [0, 0, 1]),
              'hexagonal': ([0, 0, 0, 1], None)}

FIELD_PRACTICAL_T = 2.0      # electromagnet / permanent-magnet fixtures
FIELD_HIGH_T = 10.0          # superconducting magnets used in magnetic processing
V_STAR = (1e-26, 1e-23)      # nucleus volume bounds (r ~ 1.3-13 nm), an estimate


def L_of(s):
    return Lattice(s['system'], s['a'], s.get('c'))


class Path:
    """Formation path parent -> hard phase with its variants and twins."""

    def __init__(self, ph, p):
        self.ph, self.p = ph, p
        self.mode = p['mode']                       # direct | OR | massive | spinodal
        self.LH = L_of(ph['structure'])
        self.eH = np.asarray(ph['easy_axis_vec'], float)
        self.eH = self.eH / np.linalg.norm(self.eH)
        par = p.get('parent')
        if par and self.mode in ('OR', 'spinodal', 'massive'):
            self.LP = L_of(par['structure'])
            o = p['OR']
            R = or_rotation(self.LP, self.LH, o['plane_P'], o['plane_H'], o['dir_P'], o['dir_H'])
            self.V = variant_rotations(R, par['structure']['laue'], ph['structure']['laue'])
        else:
            self.LP = self.LH
            self.V = np.eye(3)[None]
        self.twins = [(twin_rotation(self.LH, t['axis'], t['angle']), 0.0)
                      for t in p.get('twins', [])]

    @property
    def parent(self):
        return self.p.get('parent') or {}

    def axes(self, parents):
        return E.product_axes(parents, self.V, self.eH)


def _fibre(path, d):
    return E.texture_fibre(path.LP.direction(d))


def _equal_best(path, parents):
    """Equal variants; best magnet axis; for fibres check fibre axis and perpendicular."""
    ax, fr, per = path.axes(parents)
    A_par = E.alignment(ax, fr, (0, 0, 1))
    A_perp = E.alignment(ax, fr, (1, 0, 0))
    return (A_par, 'fibre axis') if A_par >= A_perp else (A_perp, 'perpendicular to fibre axis')


def _select_best(path, parents, kind='field', sign=+1, fibre=True):
    """Complete selection; bias axis scanned over polar angle to the fibre axis
    (fibre parents) or fixed (random parents)."""
    ax, fr, per = path.axes(parents)
    angles = range(0, 91, 1) if fibre else (0,)
    best = (-1.0, None)
    for al in angles:
        h = (np.sin(np.radians(al)), 0.0, np.cos(np.radians(al)))
        w = E.bias_weights(ax, fr, len(path.V), per, kind, h, sign=sign)
        A = E.alignment(ax, w, h) if (kind == 'field' or sign > 0) else E.best_axis(ax, w, 400)[0]
        if A > best[0] + 1e-9:
            best = (A, al)
    return best


def massive_or(path):
    return path.mode == 'massive'


def formation_magnetic(ph, path):
    """Is the forming phase ferromagnetic at (some part of) its formation window?"""
    Tc = ph['Tc']
    lo, hi = path.p['T_f']
    if Tc > hi:
        return 'yes', f'Tc {Tc:.0f} K above the whole formation window {lo:.0f}-{hi:.0f} K'
    if Tc > lo:
        return 'partly', (f'Tc {Tc:.0f} K inside the formation window {lo:.0f}-{hi:.0f} K: '
                          f'acts only if formed below {Tc:.0f} K')
    return 'no', f'Tc {Tc:.0f} K below the formation window {lo:.0f}-{hi:.0f} K'


def para_bound(ph, path, B):
    """Upper bound of the paramagnetic bias (kT) at the lowest formation temperature."""
    lo = max(path.p['T_f'][0], ph['Tc'] + 1)
    return P.para_bias_kT_bound(ph['Js0'], ph['n_mag'], ph['Tc'], lo, B, V_STAR[1])
