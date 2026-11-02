"""Export numerical figure data without drawing or changing frozen predictions.

Run from the repository root with PYTHONPATH=src python -m align.figdata.
Orientation populations come from the engine; measurements come from the saved
score files. Outputs are CSVs under figures/data/.
"""
import csv
import json
import os

import numpy as np

from . import engine as E
from .curate import load_phases
from .data import ROOT
from .routes import Path

OUT = os.path.join(ROOT, 'figures', 'data')
Z = np.array([0.0, 0.0, 1.0])
SPREADS = [0, 2, 5, 10, 15, 20, 30, 45, 60, 90]


def _write(name, header, rows):
    with open(os.path.join(OUT, name), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def _axes_rows(ax, fr, n_parents):
    """Fold axial directions onto the upper hemisphere and normalise per parent."""
    rows = []
    for i in range(n_parents):
        for k in range(ax.shape[1]):
            x, y, z = ax[i, k]
            if z < 0:
                x, y, z = -x, -y, -z
            rows.append([i, k, round(float(x), 6), round(float(y), 6), round(float(z), 6),
                         round(float(fr[i, k] / fr[i].sum()), 6)])
    return rows


def mnal(phases):
    """Export inherited MnAl(C) easy axes and alignment along the parent fibre."""
    ph = {p['id']: p for p in phases}['tau-MnAl-C']
    p = ph['paths'][0]
    path = Path(ph, p)
    d = path.LP.direction(p['parent']['am_texture']['dir'])
    head = ['parent', 'variant', 'x', 'y', 'z', 'weight']
    g = E.texture_fibre(d, n_az=720)
    ax, fr, _ = path.axes(g)
    _write('mnal_c_axes_fibre_parent.csv', head, _axes_rows(ax, fr, len(g)))
    gr = E.texture_random(4000)
    ax2, fr2, _ = path.axes(gr)
    _write('mnal_c_axes_random_parent.csv', head, _axes_rows(ax2, fr2, len(gr)))
    rows = []
    for s in SPREADS:
        gs = E.texture_fibre(d, n_az=720, spread_deg=s)
        a, w, _ = path.axes(gs)
        rows.append([s, round(float((np.abs(a @ Z) * w).sum() / w.sum()), 4)])
    a, w, _ = path.axes(E.texture_random(20000))
    rows.append(['random', round(float((np.abs(a @ Z) * w).sum() / w.sum()), 4)])
    _write('mnal_c_A_vs_parent_spread.csv', ['parent_fibre_spread_deg', 'A_along_fibre'], rows)


def heldout():
    """Export scored held-out and known cases with version and status labels."""
    rows = []
    for v in ('v1', 'v2', 'v3', 'v4', 'v5'):
        s = json.load(open(os.path.join(ROOT, 'results', f'scores_{v}.json')))
        for c in s['cases']:
            if c.get('result') != 'scored' or c['status'] not in ('heldout', 'known'):
                continue
            rows.append([v, c['status'], c['case'], c['phase'], c['route'], c['cls'], c['A_hi'], c['A_base'],
                         c['m_lo'], c['m_hi'], c['tol'], c['ceiling']])
    _write('measurements_vs_ceiling.csv', ['version', 'status', 'case', 'phase', 'route', 'predicted_class',
                                           'A_ceiling', 'A_base', 'A_measured_lo', 'A_measured_hi', 'tol',
                                           'ceiling'], rows)


def main():
    """Write the orientation and measurement CSVs to figures/data/."""
    os.makedirs(OUT, exist_ok=True)
    phases, _ = load_phases()
    mnal(phases)
    heldout()
    print('wrote', sorted(os.listdir(OUT)))


if __name__ == '__main__':
    main()
