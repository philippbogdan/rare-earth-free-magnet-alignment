"""Assemble headline counts, geometry, evidence totals and frozen-file hashes.

Called by align.cli after account generation and scoring. The returned mapping
is written to results/summary.json; this module has no command-line entry point.
"""
import collections
import hashlib
import json
import os

import numpy as np
import yaml

from . import account
from .data import ROOT
from .routes import ROUTES, Path, _fibre

REFERENCE_PHASES = ('Nd2Fe14B',)
LOW_CEILING = 0.97


def _sha(path):
    with open(os.path.join(ROOT, path), 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def _account(rows):
    ver = collections.Counter(r['verdict'] for r in rows)
    ref = [r for r in rows if r['phase'] not in REFERENCE_PHASES]
    phases = sorted({r['phase'] for r in rows})
    return {
        'phases': len(phases),
        'rare_earth_free_phases': len([p for p in phases if p not in REFERENCE_PHASES]),
        'reference_phases': list(REFERENCE_PHASES),
        'route_classes': len({r['route'] for r in rows}),
        'cells': len(rows),
        'verdicts': {k: ver[k] for k in ('full', 'capped', 'none', 'open')},
        'rare_earth_free_cells': len(ref),
        'rare_earth_free_verdicts': {k: sum(r['verdict'] == k for r in ref)
                                     for k in ('full', 'capped', 'none', 'open')},
        'full_threshold_A': 0.95,
        'full_by_route': dict(collections.Counter(r['route'] for r in rows if r['verdict'] == 'full').most_common()),
        'full_in_film_or_single_crystal_routes': sum(r['verdict'] == 'full' and r['route'] in ('EP', 'SC+FA', 'SC+SA')
                                                     for r in rows),
        'cells_with_refuting_measurement': sum(bool(r['test'].strip()) for r in rows),
        'open_cells_with_deciding_quantity': sum(r['verdict'] == 'open' and json.loads(r['decides']) != []
                                                 for r in rows),
    }


def _mnal(phases, rows):
    ph = {p['id']: p for p in phases}['tau-MnAl-C']
    cells = {r['route']: r for r in rows if r['phase'] == 'tau-MnAl-C'}
    p = ph['paths'][0]
    path = Path(ph, p)
    d = p['parent']['am_texture']['dir']
    ax, fr, _ = path.axes(_fibre(path, d))
    z = np.array([0.0, 0.0, 1.0])
    ang = np.degrees(np.arccos(np.clip(np.abs(ax.reshape(-1, 3) @ z), 0, 1)))
    return {
        'phase': 'tau-MnAl-C',
        'parent_fibre': d,
        'product_variants': int(ax.shape[1]),
        'easy_axis_angle_to_fibre_deg': round(float(ang.mean()), 1),
        'easy_axis_angle_spread_deg': round(float(ang.max() - ang.min()), 6),
        'A_max_AM': float(cells['AM']['A_max']),
        'A_max_DS': float(cells['DS']['A_max']),
        'verdict_AM': cells['AM']['verdict'],
        'A_random': account.ISO,
        'cos_of_angle': round(float(np.cos(np.radians(ang.mean()))), 4),
    }


def _validation(scores):
    out = {}
    for v, s in scores.items():
        h = s['summary'].get('heldout', {})
        out[v] = {'heldout_measurements': h.get('n_scored', 0), 'within': h.get('ceiling_ok', 0),
                  'straddling': h.get('ceiling_possible', 0), 'violated': h.get('ceiling_violated', 0)}
    s1 = scores['v1']
    ho = [c for c in s1['cases'] if c['status'] == 'heldout' and c.get('result') == 'scored']
    low = [c for c in ho if c['A_hi'] < LOW_CEILING]
    h = s1['summary']['heldout']
    out['v1_detail'] = {
        'within_share': round(h['ceiling_ok'] / h['n_scored'], 4),
        'violated_share': round(h['ceiling_violated'] / h['n_scored'], 4),
        'within_or_straddling_share': round((h['ceiling_ok'] + h['ceiling_possible']) / h['n_scored'], 4),
        'records_in_heldout_batches': sum(c['status'] == 'heldout' for c in s1['cases']),
        'not_scored': dict(collections.Counter(c['result'] for c in s1['cases']
                                               if c['status'] == 'heldout' and c.get('result') != 'scored')),
        'ceiling_below_%.2f' % LOW_CEILING: {
            'measurements': len(low),
            'within': sum(c['ceiling'] == 'ok' for c in low),
            'straddling': sum(c['ceiling'] == 'possible violation' for c in low),
            'violated': sum(c['ceiling'] == 'violated' for c in low)},
        'cannot_align_refuted': [h['none_refuted'], h['none_n']],
        'alignment_shown': [h['acts_shown'], h['acts_n']],
        'pairs_as_expected': [s1['summary']['pairs']['ok'], s1['summary']['pairs']['n']],
        'violations': sorted(c['case'] for c in ho if c['ceiling'] == 'violated'),
        'exposed': s1['summary'].get('exposed'),
        'known': s1['summary'].get('known'),
    }
    return out


def _evidence(reg, R):
    c = collections.Counter()
    for rid in R.used:
        v = reg.ver.get(rid) or {}
        if v.get('quote_ok'):
            c['quote_found_' + (v.get('level') or 'online')] += 1
        elif v:
            c['quote_not_found'] += 1
        else:
            c['not_checked'] += 1
        if v.get('doi_ok'):
            c['doi_resolves'] += 1
    found = sum(n for k, n in c.items() if k.startswith('quote_found_'))
    return {'records_used_by_account': len(R.used), 'quote_found_at_source': found, **dict(sorted(c.items())),
            'labelled_estimates': len(R.estimates),
            'records_in_repository': len(reg.rec)}


def _coverage(phases, reg):
    univ = yaml.safe_load(open(os.path.join(ROOT, 'data', 'curated', 'universe.yaml')))
    names = [p['name'] for p in univ['phases']]
    lit = [r['class'] for r in univ['route_classes']]
    return {'candidate_phases_named_in_reviews': len(names),
            'candidate_phases_covered': len({p.get('universe') for p in phases} & set(names)),
            'reviews': len({(reg.rec[r].get('source') or {}).get('doi') or r
                            for p in univ['phases'] for r in p['named_in']}),
            'literature_route_classes': len(lit),
            'literature_route_classes_covered': len([r for r in lit if r in ROUTES])}


def _frozen():
    man = json.load(open(os.path.join(ROOT, 'data', 'heldout', 'MANIFEST.json')))
    files = {'v1': 'cases.yaml', 'v2': 'cases_v2.yaml', 'v3': 'cases_v3.yaml', 'v4': 'cases_v4.yaml',
             'v5': 'cases_v5.yaml'}
    return {
        'predictions_sha256': {v: _sha(os.path.join('data', 'curated', f)) for v, f in files.items()},
        'heldout_batches': {b: {f: m['n_records'] for f, m in d['files'].items()}
                            for b, d in man['batches'].items()},
    }


def build(phases, rows, scores, reg, R):
    """Collect headline fields from evaluated rows, saved-version scores and evidence."""
    return {
        'account': _account(rows),
        'mnal_c_textured_parent': _mnal(phases, rows),
        'validation': _validation(scores),
        'coverage': _coverage(phases, reg),
        'evidence': _evidence(reg, R),
        'frozen': _frozen(),
    }
