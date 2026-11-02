"""Write the route account, scores of frozen predictions and numerical summary.

Run from the repository root with PYTHONPATH=src python -m align.cli.
Use python reproduce.py to run checks and compare result hashes as well.
"""
import csv
import json
import os
import sys

from . import account, evaluate, summary
from .curate import load_phases
from .data import ROOT, Registry

RES = os.path.join(ROOT, 'results')

ROUTE_INPUTS = {          # which parts of a phase description each route reads
    'FP': ['Js', 'particles', 'consolidation'], 'MO': ['Js', 'particles'],
    'SC+FA': ['Js', 'Tc', 'paths'], 'SC+SA': ['Js', 'paths'],
    'FS': ['Js', 'melt', 'paths'],
    'DS': ['Js', 'paths'], 'AM': ['Js', 'paths'],
    'FA': ['Js', 'Tc', 'paths'], 'DS+FA': ['Js', 'Tc', 'paths'], 'AM+FA': ['Js', 'Tc', 'paths'],
    'RX+FA': ['Js', 'Tc', 'paths'], 'SA': ['Js', 'paths'], 'DS+SA': ['Js', 'paths'],
    'RX+SA': ['Js', 'paths'], 'HD': ['Js', 'hd'], 'HDDR': ['Js', 'hddr'], 'EP': ['Js', 'epitaxy'],
    'TGG': ['Js', 'tgg'],
}

def cell_refs(ph, route):
    """Collect the source record IDs read by a phase and route combination."""
    from .curate import refs_in
    raw = ph.get('_raw', {})
    out = []
    for k in ROUTE_INPUTS[route] + ['structure', 'easy_axis']:
        out += refs_in(raw.get(k))
    return sorted(set(out))

def matrix(phases):
    """Return account rows with mechanisms, alternatives, source IDs and tests."""
    rows = []
    for ph in phases:
        for c in account.evaluate_with_tests(ph):
            phi = 1.0 if c.route == 'EP' else float((ph.get('phi_routes') or {}).get(
                c.route, ph.get('phi_max', 1.0)))
            Br = ph['Js'] * c.A_max * phi if ph.get('Js') is not None else None
            rows.append({'phase': ph['id'], 'route': c.route, 'verdict': c.verdict,
                         'A_max': round(c.A_max, 4), 'phi_max': round(phi, 4),
                         'Br_max_T': None if Br is None else round(Br, 4),
                         'axis': c.axis, 'chain': c.chain,
                         'prereq': '; '.join(f'{n}: {"yes" if ok else ("no" if ok is False else "?")} ({w})'
                                             for n, ok, w in c.prereq),
                         'decides': json.dumps([[d[0], {k: (None if v is None else round(v, 4))
                                                        for k, v in d[1].items()}, d[2]]
                                                for d in c.decides]),
                         'refs': ' '.join(cell_refs(ph, c.route)),
                         'test': c.test, 'A_alt': json.dumps({k: (None if v is None else round(v, 4))
                                                               for k, v in c.A_alt.items()})})
    return rows

def write_csv(fn, rows):
    if not rows:
        return
    with open(fn, 'w', newline='') as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def main():
    """Regenerate results using current phase data and the frozen case files."""
    os.makedirs(RES, exist_ok=True)
    reg = Registry()
    phases, R = load_phases(reg)
    if R.missing:
        print('MISSING record references:', R.missing)
        sys.exit(1)
    rows = matrix(phases)
    for ph in phases:
        for c in account.evaluate(ph):
            R.used.update(r for r in (c.refs or []) if r in reg.rec)
    write_csv(os.path.join(RES, 'account.csv'), rows)
    scores = {v: evaluate.run(v, phases) for v in evaluate.FILES}
    out = summary.build(phases, rows, scores, reg, R)
    with open(os.path.join(RES, 'summary.json'), 'w') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
        f.write('\n')
    print(f"wrote results/account.csv ({len(rows)} cells, {len(phases)} phases including Nd2Fe14B reference), "
          f"results/scores_v1.json to scores_v5.json and results/summary.json")
    return out

if __name__ == '__main__':
    main()
