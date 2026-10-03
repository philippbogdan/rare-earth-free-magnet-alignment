"""Frontier bookkeeping: phases covered, route classes covered, predictive accuracy.

results/frontier_history.csv holds one row per account version (written when a version
is frozen); the current row is recomputed by the regeneration command.
"""
import csv
import os

from .data import ROOT

HIST = os.path.join(ROOT, 'results', 'frontier_history.csv')
FIELDS = ['version', 'commit', 'phases_covered', 'phases_universe', 'routes_covered',
          'routes_universe', 'n_cases', 'accuracy', 'sharpness', 'ceiling_violations',
          'pushed_next', 'note']


def read_history():
    if not os.path.exists(HIST):
        return []
    return list(csv.DictReader(open(HIST)))


def append_history(row):
    new = not os.path.exists(HIST)
    with open(HIST, 'a', newline='') as f:
        w = csv.DictWriter(f, FIELDS)
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, '') for k in FIELDS})


def scores(row):
    """Each criterion as a fraction of its limit (1 = at the limit)."""
    p = float(row['phases_covered']) / float(row['phases_universe'])
    r = float(row['routes_covered']) / float(row['routes_universe'])
    a = float(row['accuracy']) if row.get('accuracy') not in (None, '') else 0.0
    return {'phases': p, 'routes': r, 'accuracy': a}


def furthest(row):
    s = scores(row)
    return min(s, key=s.get)


def non_dominated(rows):
    pts = [(r, scores(r)) for r in rows]
    keep = []
    for r, s in pts:
        dominated = any(all(t[k] >= s[k] for k in s) and any(t[k] > s[k] for k in s)
                        for _, t in pts)
        if not dominated:
            keep.append(r)
    return keep
