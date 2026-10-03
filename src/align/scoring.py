"""Scoring of the conditional predictions (data/curated/cases.yaml) against measured outcomes.

Fixed before any sealed batch is opened (docs/PROTOCOL_LOG.md). For each case:

  measurement -> A interval [m_lo, m_hi]
      Mr/Ms, A_reported, variant fraction : value -/+ its stated uncertainty
      Br                                  : Br / (Js * phi), with Js from the same-sample record
                                            (js_case), the case's sourced interval (js), or the
                                            phase's Js -/+ 5 %, and phi from the case (default 1)
      rocking FWHM / spread angle (deg)   : <|cos|> of a 2-D Gaussian tilt (sigma = FWHM/2.3548)
      other texture metrics, (BH)max, Hc  : not converted (not scored)
  tolerance  tol = max(DELTA, 2 x stated uncertainty in A)

  ceiling       violated (certain)  if m_lo > A_hi + tol
                violated (possible) if m_hi > A_hi + tol >= m_lo  (conversion range straddles it)
  class 'none'  refuted if m_lo > A_hi + tol   (alignment beyond what no mechanism allows)
  class 'acts'  shown   if m_lo > A_base + tol; otherwise 'not shown' (not a refutation)
  class 'open'  each hypothesis h kept if m_lo <= hyp[h] + tol, else refuted
  interval      (pre-registered) inside if [m_lo, m_hi] overlaps [min(A_base, A_hi) - tol, A_hi + tol];
                sharpness = mean (A_hi - A_base)
  pairs         'higher': value(case) > value(reference) beyond the stated uncertainties;
                'equal' : |A(case) - A(reference)| <= 2 DELTA
"""
import glob
import math
import os

import numpy as np
import yaml

from .account import DELTA
from .data import ROOT

UNIT_T = {'T': 1.0, 'kG': 0.1, 'kGs': 0.1, 'mT': 1e-3, 'G': 1e-4, 'Gs': 1e-4}


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def load_outcomes():
    """Opened outcome records (data/outcomes/*.yaml) and known cases (data/raw/known_cases_*.yaml)."""
    recs = {}
    for f in sorted(glob.glob(os.path.join(ROOT, 'data', 'outcomes', '*.yaml'))):
        for r in yaml.safe_load(open(f)) or []:
            recs[r['id']] = r
    known = {}
    for f in sorted(glob.glob(os.path.join(ROOT, 'data', 'raw', 'known_cases_*.yaml'))):
        for r in yaml.safe_load(open(f)) or []:
            known[r['id']] = r
    return recs, known


def tilt_A(sigma_deg, n=200000, seed=7):
    """<|cos theta|> for a tilt whose two components are N(0, sigma)."""
    if sigma_deg <= 0:
        return 1.0
    rng = np.random.default_rng(seed)
    t = np.radians(rng.normal(scale=sigma_deg, size=(n, 2)))
    th = np.hypot(t[:, 0], t[:, 1])
    return float(np.abs(np.cos(th)).mean())


def measured_A(spec, rec, recs, phases):
    """(m_lo, m_hi, note) or None when the record does not convert to A."""
    kind = spec['meas']['kind']
    v = _num(rec.get('value'))
    if v is None:
        return None
    u = _num(rec.get('uncertainty')) or 0.0
    unit = str(rec.get('unit') or '').strip()
    if kind in ('MrMs', 'A', 'vfrac'):
        if '%' in unit or 'percent' in unit.lower():
            v, u = v / 100, u / 100
        return v - u, v + u, ''
    if kind == 'Br':
        f = UNIT_T.get(unit)
        if f is None:
            return None
        br_lo, br_hi = (v - u) * f, (v + u) * f
        m = spec['meas']
        if m.get('js_case') and m['js_case'] in recs:
            jr = recs[m['js_case']]
            jv = _num(jr.get('value'))
            ju = _num(jr.get('uncertainty')) or 0.0
            jf = UNIT_T.get(str(jr.get('unit') or 'T').strip())
            if jv is None or jf is None:          # e.g. a mass magnetisation: no conversion to T
                return None
            js = [(jv - ju) * jf, (jv + ju) * jf]
            src = 'same-sample Js'
        elif m.get('js'):
            js, src = m['js'], 'Js interval of the grade'
        else:
            J = (phases.get(spec['phase']) or {}).get('Js')
            if J is None:
                return None
            js, src = [0.95 * J, 1.05 * J], 'phase Js -/+5%'
        phi = m.get('phi', [1.0, 1.0])
        return br_lo / (js[1] * phi[1]), br_hi / (js[0] * phi[0]), src
    if kind == 'class':
        if 'fwhm' in unit.lower():
            return tilt_A(v / 2.3548) - 0.0, tilt_A(v / 2.3548), 'FWHM -> 2-D Gaussian tilt'
        if 'spread_deg' in unit.lower():
            return tilt_A(v), tilt_A(v), 'spread -> 2-D Gaussian tilt'
    return None


def score_case(spec, rec, recs, phases):
    p = spec['pred']
    out = {'case': spec['case'], 'phase': spec['phase'], 'route': spec['route'], 'topic': spec.get('topic'),
           'status': spec['status'], 'cls': p['cls'], 'A_hi': p['A_hi'], 'A_base': p['A_base'],
           'why': p['why']}
    if p['cls'] == 'info':
        out['result'] = 'not a test'
        return out
    if not spec['meas'].get('saturating', True):
        out['result'] = 'not saturated'
        return out
    m = measured_A(spec, rec, recs, phases)
    if m is None:
        out['result'] = 'not convertible'
        return out
    lo, hi, note = m
    u = (_num(rec.get('uncertainty')) or 0.0)
    tol = max(DELTA, 2 * u if spec['meas']['kind'] in ('MrMs', 'A', 'vfrac') else DELTA)
    out.update({'m_lo': round(lo, 4), 'm_hi': round(hi, 4), 'conv': note, 'tol': tol})
    out['ceiling'] = 'violated' if lo > p['A_hi'] + tol else (
        'possible violation' if hi > p['A_hi'] + tol else 'ok')
    if p['cls'] == 'none':
        out['class'] = 'refuted' if lo > p['A_hi'] + tol else 'consistent'
    elif p['cls'] == 'acts':
        out['class'] = 'shown' if lo > p['A_base'] + tol else 'not shown'
    elif p['cls'] == 'open':
        out['class'] = 'open'
    if p.get('hyp'):
        out['hyp'] = {h: ('kept' if lo <= a + tol else 'refuted') for h, a in p['hyp'].items()}
    a_lo = min(p['A_base'], p['A_hi']) - tol
    out['interval'] = 'inside' if (hi >= a_lo and lo <= p['A_hi'] + tol) else 'outside'
    out['result'] = 'scored'
    return out


def score_pairs(specs, scored, recs):
    by = {s['case']: s for s in scored}
    res = []
    for sp in specs:
        pr = sp.get('pair')
        if not pr or pr['expect'] == 'open':
            continue
        a, b = sp['case'], pr['with']
        ra, rb = recs.get(a), recs.get(b)
        if ra is None or rb is None:
            continue
        va, vb = _num(ra.get('value')), _num(rb.get('value'))
        if va is None or vb is None or str(ra.get('unit')) != str(rb.get('unit')):
            continue
        ua, ub = _num(ra.get('uncertainty')) or 0.0, _num(rb.get('uncertainty')) or 0.0
        if pr['expect'] == 'higher':
            ok = va - ua > vb + ub
        else:
            sa, sb = by.get(a, {}), by.get(b, {})
            if 'm_lo' in sa and 'm_lo' in sb:
                ok = abs((sa['m_lo'] + sa['m_hi']) / 2 - (sb['m_lo'] + sb['m_hi']) / 2) <= 2 * DELTA
            else:
                ok = abs(va - vb) <= 2 * DELTA * max(abs(vb), 1e-9)
        res.append({'case': a, 'reference': b, 'expect': pr['expect'], 'ok': bool(ok)})
    return res


def summarise(scored, pairs):
    out = {}
    for st in sorted({s['status'] for s in scored}):
        ss = [s for s in scored if s['status'] == st and s.get('result') == 'scored']
        n = len(ss)
        out[st] = {
            'n_scored': n,
            'n_not_test': sum(1 for s in scored if s['status'] == st and s.get('result') != 'scored'),
            'ceiling_ok': sum(s['ceiling'] == 'ok' for s in ss),
            'ceiling_possible': sum(s['ceiling'] == 'possible violation' for s in ss),
            'ceiling_violated': sum(s['ceiling'] == 'violated' for s in ss),
            'none_n': sum(s['cls'] == 'none' for s in ss),
            'none_refuted': sum(s.get('class') == 'refuted' for s in ss),
            'acts_n': sum(s['cls'] == 'acts' for s in ss),
            'acts_shown': sum(s.get('class') == 'shown' for s in ss),
            'open_n': sum(s['cls'] == 'open' for s in ss),
            'interval_inside': sum(s['interval'] == 'inside' for s in ss),
            'sharpness': round(float(np.mean([s['A_hi'] - s['A_base'] for s in ss])), 4) if n else None,
        }
    out['pairs'] = {'n': len(pairs), 'ok': sum(p['ok'] for p in pairs)}
    return out
