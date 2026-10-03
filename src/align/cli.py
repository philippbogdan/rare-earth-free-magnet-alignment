"""The single regeneration command: ./regenerate.sh  (PYTHONPATH=src python -m align.cli)."""
import csv
import json
import os
import sys

import yaml

from . import account, report
from .curate import load_phases
from .data import ROOT, Registry
from .stamp import stamp

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
    from .curate import refs_in
    raw = ph.get('_raw', {})
    out = []
    for k in ROUTE_INPUTS[route] + ['structure', 'easy_axis']:
        out += refs_in(raw.get(k))
    return sorted(set(out))


def matrix(phases):
    rows = []
    for ph in phases:
        for c in account.evaluate_with_tests(ph):
            Br = ph['Js'] * c.A_max if ph.get('Js') is not None else None
            rows.append({'phase': ph['id'], 'route': c.route, 'verdict': c.verdict,
                         'A_max': round(c.A_max, 4), 'Br_max_T': None if Br is None else round(Br, 4),
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
    os.makedirs(RES, exist_ok=True)
    st = stamp()
    reg = Registry()
    phases, R = load_phases(reg)
    if R.missing:
        print('MISSING record references:', R.missing)
        sys.exit(1)
    rows = matrix(phases)
    write_csv(os.path.join(RES, 'account.csv'), rows)
    ctx = report.build(reg, phases, rows, R, st)
    json.dump({'stamp': st, 'n_phases': len(phases), 'n_cells': len(rows)},
              open(os.path.join(RES, 'stamp.json'), 'w'), indent=1)
    report.write_docs(ctx)
    print(f"regenerated {len(rows)} cells for {len(phases)} phases from commit {st['commit']}"
          + ('' if st['clean_inputs'] else ' (UNCOMMITTED INPUTS: ' + ', '.join(st['dirty_files']) + ')'))


if __name__ == '__main__':
    main()
