"""Smoke test of the route rules on a FIXTURE (round numbers, not sourced data)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from align.account import evaluate

FIXTURE = {
    'id': 'fixture_MnAl_like', 'structure': {'system': 'tetragonal', 'a': 3.92, 'c': 3.57, 'laue': '4/mmm'},
    'easy_axis_vec': [0, 0, 1], 'Tc': 650.0, 'Js0': 0.8, 'n_mag': 4.1e28,
    'particles': {'single_variant': 'unknown', 'internal': 'path', 'path': 0, 'reason': 'twinned'},
    'paths': [{'mode': 'massive',
               'parent': {'structure': {'system': 'hexagonal', 'a': 2.69, 'c': 4.38, 'laue': '6/mmm'},
                          'family': 'hexagonal', 'am_texture': {'dir': [0, 0, 0, 1]},
                          'ds_texture': {'dir': [0, 0, 0, 1]}, 'rollable': False},
               'OR': {'plane_P': [0, 0, 0, 1], 'plane_H': [1, 1, 1], 'dir_P': [1, 1, -2, 0], 'dir_H': [1, -1, 0]},
               'T_f': [673.0, 973.0], 'stress_sign': 0}],
    'melt': {'primary': 'parent', 'path': 0},
    'hd': {'fibre': [0, 0, 1], 'reason': 'fixture'},
    'epitaxy': {'demonstrated': True}, 'tgg': {},
}

def test_fixture():
    cells = evaluate(FIXTURE)
    for c in cells:
        print(f'{c.route:6s} {c.verdict:7s} A_max={c.A_max:.3f} {c.axis:40s} | {c.chain[:110]}')
        for d in c.decides:
            print('        decides:', d[0], {k: (round(v, 3) if v is not None else None) for k, v in d[1].items()})

if __name__ == '__main__':
    test_fixture()
