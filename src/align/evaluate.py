"""Score a set of case predictions against the opened outcomes and the known cases.

  v1  data/curated/cases.yaml       frozen before opening (docs/PROTOCOL_LOG.md): the held-out test
  v2  data/curated/cases_v2.yaml    revised after opening (post hoc; every change listed in
                                    docs/PROTOCOL_LOG.md); scored for information, never as held-out
Writes results/scores_<version>.json.
"""
import json
import os

import yaml

from . import scoring
from .curate import load_phases
from .data import ROOT

FILES = {'v1': 'cases.yaml', 'v2': 'cases_v2.yaml', 'v3': 'cases_v3.yaml', 'v4': 'cases_v4.yaml',
         'v5': 'cases_v5.yaml'}
# topics whose outcomes were opened only after the version's predictions were frozen
HELDOUT_AFTER = {'v1': None, 'v2': {'mn4n'}, 'v3': {'heldout2_b'}, 'v4': {'heldout2_a'}, 'v5': {'heldout3'}}


def run(version, phases=None):
    path = os.path.join(ROOT, 'data', 'curated', FILES[version])
    if not os.path.exists(path):
        return None
    if phases is None:
        phases, _ = load_phases()
    P = {p['id']: p for p in phases}
    recs, known = scoring.load_outcomes()
    specs = yaml.safe_load(open(path)) or []
    kpath = os.path.join(ROOT, 'data', 'curated', 'cases_known.yaml' if version == 'v1' else f'cases_known_{version}.yaml')
    kspecs = yaml.safe_load(open(kpath)) if os.path.exists(kpath) else []
    rows = [scoring.score_case(s, recs[s['case']], recs, P) for s in specs if s['case'] in recs]
    rows += [scoring.score_case(s, known[s['case']], known, P) for s in kspecs if s['case'] in known]
    ho = HELDOUT_AFTER.get(version)
    if ho is not None:
        for x in rows:
            if x['status'] == 'heldout' and x.get('topic') not in ho:
                x['status'] = 'post hoc'
    pairs = scoring.score_pairs(specs, rows, recs)
    out = {'version': version, 'file': FILES[version], 'summary': scoring.summarise(rows, pairs),
           'pairs': pairs, 'cases': rows,
           'unopened': sorted(s['case'] for s in specs if s['case'] not in recs)}
    if version == 'v3':
        out['sensitivity_js'] = js_sensitivity(specs, recs, P)
    os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
    json.dump(out, open(os.path.join(ROOT, 'results', f'scores_{version}.json'), 'w'), indent=1,
              ensure_ascii=False)
    return out


def js_sensitivity(specs, recs, P):
    """Post-hoc check of a v3 conversion bug: Br cases outside AlNiCo whose same-sample Js was
    not linked in the frozen specs, rescored with it (the frozen result is not changed)."""
    import copy
    from . import cases
    now = {x['case']: x for x in cases.build()}
    out = []
    for sp in specs:
        cid = sp['case']
        nw = now.get(cid)
        if cid not in recs or sp['meas'].get('js_case') or not nw or not nw['meas'].get('js_case'):
            continue
        sp2 = copy.deepcopy(sp)
        sp2['meas']['js_case'] = nw['meas']['js_case']
        a, b = scoring.score_case(sp, recs[cid], recs, P), scoring.score_case(sp2, recs[cid], recs, P)
        out.append({'case': cid, 'status': sp['status'], 'cls': sp['pred']['cls'],
                    'frozen': [a.get('m_lo'), a.get('m_hi'), a.get('ceiling'), a.get('class')],
                    'same_sample_js': [b.get('m_lo'), b.get('m_hi'), b.get('ceiling'), b.get('class')]})
    return out


if __name__ == '__main__':
    import sys
    for v in sys.argv[1:] or ['v1', 'v2']:
        r = run(v)
        if r:
            print(v, json.dumps(r['summary']))
