"""Route evaluation: the same rules for every phase.

A phase description (from data/curated/phases.yaml, every number tied to a source record)
is turned into, for each route: a verdict, the ceiling A_max with its magnet axis, the
chain that produces it, the prerequisites checked, any unpublished quantity that decides
the answer (with both outcomes), and the single measurement that would test it.
"""
import numpy as np

from . import engine as E
from .orient import Lattice
from .variants import or_rotation, variant_rotations, twin_rotation

ROUTES = ['FP', 'FS', 'DS', 'AM', 'FA', 'DS+FA', 'AM+FA', 'RX+FA', 'SA', 'DS+SA', 'RX+SA',
          'HD', 'EP', 'TGG']

ROUTE_NAMES = {
    'FP': 'field alignment of free particles + consolidation (dry/wet pressing, bonding, colloidal assembly)',
    'FS': 'solidification in a magnetic field (crystals rotate in the melt)',
    'DS': 'directional solidification / columnar casting, then transformation without bias',
    'AM': 'additive manufacturing (texture of the printed phase), then transformation without bias',
    'FA': 'field annealing: formation of the hard phase from an untextured parent in a field',
    'DS+FA': 'columnar (DS) parent + formation in a field',
    'AM+FA': 'printed (AM) parent + formation in a field',
    'RX+FA': 'rolled/drawn and recrystallised parent + formation in a field',
    'SA': 'stress-assisted formation from an untextured parent',
    'DS+SA': 'columnar (DS) parent + formation under stress',
    'RX+SA': 'rolled/drawn and recrystallised parent + formation under stress',
    'HD': 'hot deformation of the hard phase (extrusion, die-upsetting, rolling)',
    'EP': 'epitaxial growth on a single-crystal substrate (films; not a bulk magnet)',
    'TGG': 'templated grain growth on aligned seed crystals during sintering',
}

CUBIC_LAUE = ('m-3m',)


def lattice_of(s):
    return Lattice(s['system'], s['a'], s.get('c'))


def easy_vec(ph):
    return np.asarray(ph['easy_axis_vec'], float)


class Path:
    """A formation path parent -> hard phase."""

    def __init__(self, ph, path):
        self.ph, self.p = ph, path
        self.mode = path['mode']                    # 'direct' | 'OR' | 'massive' | 'spinodal'
        self.LH = lattice_of(ph['structure'])
        self.eH = easy_vec(ph)
        if self.mode in ('OR', 'spinodal'):
            par = path['parent']
            self.LP = lattice_of(par['structure'])
            o = path['OR']
            R = or_rotation(self.LP, self.LH, o['plane_P'], o['plane_H'], o['dir_P'], o['dir_H'])
            self.V = variant_rotations(R, par['structure']['laue'], ph['structure']['laue'])
        else:
            self.LP = self.LH
            self.V = np.eye(3)[None]
        self.twins = [(twin_rotation(self.LH, t['axis'], t['angle']), t.get('fraction', 0.0))
                      for t in path.get('twins', [])]

    def axes(self, parents, with_twins=False):
        tw = self.twins if with_twins else ()
        return E.product_axes(parents, self.V, self.eH, tw)


def fibre_parents(path, growth_dir):
    """Perfect fibre of the parent growth direction (crystal Cartesian)."""
    d = path.LP.direction(growth_dir) if len(growth_dir) in (3, 4) else growth_dir
    return E.texture_fibre(d)


def a_equal(path, parents):
    axes, frac, per = path.axes(parents)
    return axes, frac, per


def best_field_angle(path, parents, kind='field', sign=+1, angles=range(0, 91, 2)):
    """Complete selection with the field/stress axis at polar angle alpha to the fibre
    axis; returns (A, alpha) maximising A along the bias axis."""
    axes, frac, per = path.axes(parents)
    best = (-1, None)
    for al in angles:
        h = (np.sin(np.radians(al)), 0.0, np.cos(np.radians(al)))
        w = E.bias_weights(axes, frac, len(path.V), per, kind, h, sign=sign)
        if kind == 'stress' and sign < 0:
            A, _ = E.best_axis(axes, w, 600)
        else:
            A = E.alignment(axes, w, h)
        if A > best[0] + 1e-9:
            best = (A, al)
    return best
