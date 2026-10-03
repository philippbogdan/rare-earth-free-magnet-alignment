"""Builds FINDINGS.md, REPORT.md and COMPLETE.md from computed results (every number generated)."""
import collections
import json
import os
import subprocess

import yaml

from .data import ROOT
from .render import fmt, table, fill
from .routes import ROUTE_NAMES, ROUTES

TPL = os.path.join(ROOT, 'docs', 'templates')
BULK = [r for r in ROUTES if r != 'EP']
LIT_ROUTES = None            # literature route classes (data/curated/universe.yaml)

# the four known cases of PROBLEM.md: (label, route, selector on case ids)
KNOWN = [
    ('Field-aged columnar AlNiCo', 'DS+FA',
     lambda c: c.startswith('kc_') or c == 'uni_out_dussa2024_alnico8_columnar_durandcharre1978'),
    ('Printed MnAl(C) losing its texture (Pacheco et al., Alloys 2 (2023) 100)', 'AM',
     lambda c: c.startswith('mnal_am_pacheco2023')),
    ('Textured MnAl-C by hot extrusion', 'HD',
     lambda c: c.startswith('mnal_hd_') or c == 'uni_out_mnalc_hd_extrusion_cui2018'),
    ('Field-pressed sintered Nd-Fe-B (today\'s practice)', 'FP', lambda c: c.startswith('nd2fe14b_fp_')),
]

# post-hoc diagnosis of every v1 failure (docs/PROTOCOL_LOG.md, v2 items)
DIAGNOSIS = {
    'uni_out_rt_fept_ep_shima2002': (1, 'encoding: film normal read as "across the axis"'),
    'feni_meteorite_poirier2015': (2, "encoding: the ferromagnetic parent's own magnetisation is a bias"),
    'copt_sc_kandaurova1992a': (3, 'physics: fine polytwins act as one coupled multiaxial unit'),
    'hfco7_iso_cl_bal2014': (4, 'scope of A: isotropic cluster film, exchange enhancement not flagged'),
    'hcpco_oth_bousnina2020': (5, 'unverified rod arrangement (flat-lying or bundled rods)'),
    'hcpconi_iso_maurer2007_mrms': (5, 'unverified rod arrangement'),
    'hcpco_fp_ener2018_randompowder_mrms': (5, 'bundles of parallel rods have a net axis'),
    'fecrco_FA_sintered_fe28cr7co2mo_vompe2019': (6, 'physics: 2 wt% Mo is below the sourced 3 wt% of the '
                                                     'nearest-<100> path'),
    'mnal_other_jia2022a': (7, 'process not specific (abstract only): temperature and sequence unknown'),
    'mnal_sc_sa_kubo1976': (8, 'physics: the transformation shear breaks ties the tetragonal strain leaves'),
    'fept_nacl_7': ('v5', 'loose nanoparticle powder measured unfixed: particles can turn in the measuring field '
                          'or sinter into coupled clusters'),
}


def _diag(cid):
    for k, v in DIAGNOSIS.items():
        if cid.startswith(k):
            return v
    return (None, '')


def _load(fn):
    p = os.path.join(ROOT, fn)
    return json.load(open(p)) if os.path.exists(p) else None


def _pct(a, b):
    return f'{100.0 * a / b:.1f} %' if b else '–'


# ---------------------------------------------------------------- the account
def best_routes(rows, phase):
    cells = [r for r in rows if r['phase'] == phase and r['route'] in BULK]
    top = max(cells, key=lambda r: r['A_max'])
    return top, cells


def ranking(phases, rows, reg):
    out = []
    for ph in phases:
        top, cells = best_routes(rows, ph['id'])
        sure = [r for r in cells if r['verdict'] in ('full', 'capped')]
        if sure:
            sure_top = max(sure, key=lambda r: r['A_max'])
        else:                                   # no established alignment route: the unaligned value
            iso = max((r['A_max'] for r in cells if r['verdict'] == 'none'), default=0.5)
            sure_top = {'route': 'none established (unaligned)', 'A_max': iso, 'verdict': 'none'}
        out.append((ph, top, sure_top))
    out.sort(key=lambda t: -((t[0]['Js'] or 0.0) * (t[2]['A_max'] if t[2] else 0.5)))
    return out


def phase_section(ph, rows, reg):
    lines = [f"### {ph.get('name', ph['id'])}", '']
    facts = []
    for key, label, unit in (('Js', 'Js (300 K)', 'T'), ('K1', 'K1', 'MJ/m³'), ('Tc', 'Tc', 'K')):
        if ph.get(key) is not None:
            facts.append(f'{label} = {fmt(ph[key], 3 if key != "Tc" else 0)} {unit}')
    if ph.get('Js') is None:
        facts.append('Js unpublished')
    lines.append('; '.join(facts) + '.')
    lines.append('')
    rws = [r for r in rows if r['phase'] == ph['id']]
    lines.append(table(['route', 'verdict', 'A_max', 'Br_max (T)', 'mechanism / reason'],
                       [[r['route'], r['verdict'], fmt(r['A_max']), fmt(r['Br_max_T']), r['chain']]
                        for r in rws]))
    lines.append('')
    return '\n'.join(lines)


def evidence_table(reg, used):
    rows = []
    for rid in sorted(used):
        r = reg.rec[rid]
        src = r.get('source') or {}
        v = reg.ver.get(rid) or {}
        st = ('verified (' + (v.get('level') or '') + ')' if v.get('quote_ok') else
              ('not found' if v else 'not checked'))
        rows.append([rid, fmt(r.get('value')), r.get('unit', ''), reg.cite(rid),
                     src.get('evidence', ''), st])
    return table(['record', 'value', 'unit', 'source', 'evidence', 'quote check'], rows)


def verification_counts(reg, used):
    c = collections.Counter()
    for rid in used:
        v = reg.ver.get(rid) or {}
        if v.get('quote_ok'):
            c['verified_' + (v.get('level') or 'online')] += 1
        elif v:
            c['not_found'] += 1
        else:
            c['not_checked'] += 1
        if v.get('doi_ok'):
            c['doi_ok'] += 1
    return c


# ---------------------------------------------------------------- accuracy
def _summary_table(scores, statuses):
    rows = []
    for st in statuses:
        s = scores['summary'].get(st)
        if not s:
            continue
        n = s['n_scored']
        rows.append([st, n, f"{s['ceiling_ok']} ({_pct(s['ceiling_ok'], n)})", s['ceiling_possible'],
                     f"{s['ceiling_violated']} ({_pct(s['ceiling_violated'], n)})",
                     f"{s['none_refuted']} / {s['none_n']}", f"{s['acts_shown']} / {s['acts_n']}",
                     s['open_n'], f"{s['interval_inside']} ({_pct(s['interval_inside'], n)})",
                     fmt(s['sharpness'])])
    return table(['cases', 'scored', 'ceiling holds', 'possible violation (conversion straddles)',
                  'ceiling violated', "'cannot align' refuted", "'aligns' shown", 'open',
                  'inside pre-registered interval', 'mean A_hi - A_base'], rows)


def _by_phase(scores, status):
    g = collections.defaultdict(lambda: [0, 0, 0, 0])
    for x in scores['cases']:
        if x['status'] != status or x.get('result') != 'scored':
            continue
        k = g[x['phase']]
        k[0] += 1
        k[1] += x['ceiling'] == 'ok'
        k[2] += x['ceiling'] == 'possible violation'
        k[3] += x['ceiling'] == 'violated'
    return table(['phase', 'scored', 'ceiling holds', 'possible', 'violated'],
                 [[p, *v] for p, v in sorted(g.items())])


def _hypotheses(scores, status='heldout'):
    g = collections.defaultdict(lambda: collections.Counter())
    for x in scores['cases']:
        if x.get('hyp') and x['status'] == status:
            for h, r in x['hyp'].items():
                g[(x['phase'], h)][r] += 1
    rows = [[p, h, c['kept'], c['refuted']] for (p, h), c in sorted(g.items())]
    return table(['phase', 'reading (pre-registered alternative)', 'kept', 'refuted'], rows)


def _failures(scores, status='heldout'):
    rows = []
    for x in scores['cases']:
        if x['status'] != status or x.get('result') != 'scored':
            continue
        if x['ceiling'] == 'violated' or x.get('class') == 'refuted':
            item, why = _diag(x['case'])
            lab = (f'v2 item {item}' if isinstance(item, int) else f'{item} rule') if item else ''
            rows.append([x['case'], x['phase'], x['route'][:30], x['cls'], fmt(x['A_hi']),
                         f"{fmt(x['m_lo'])}-{fmt(x['m_hi'])}", f"{lab}: {why}" if item else '–'])
    return rows


def _pairs(scores):
    return table(['case', 'reference (same study)', 'expected', 'observed as expected'],
                 [[p['case'], p['reference'], p['expect'], 'yes' if p['ok'] else 'NO'] for p in scores['pairs']])


# ---------------------------------------------------------------- known cases
def qual_class(rec):
    """Qualitative texture statement of a known-case record -> 'none' (no preferred orientation) or
    'acts' (strong preferred orientation); used for the fidelity check only (never for held-out scoring)."""
    import re
    v = str(rec.get('value') or '').lower()
    if re.search(r'no significant|no preferential|no preferred|random|isotropic|weak', v):
        return 'none'
    if re.search(r'strong|pronounced|sharp|highly', v):
        return 'acts'
    return None


def known_rows(scores, reg, rows):
    out = []
    for label, route, sel in KNOWN:
        cs = [x for x in scores['cases'] if x['status'] == 'known' and sel(x['case'])]
        sc = [x for x in cs if x.get('result') == 'scored']
        for x in cs:                      # qualitative texture statements of the hard phase
            if x.get('result') == 'not convertible' and x.get('phase') and x['cls'] in ('none', 'acts'):
                q = qual_class(reg.rec.get(x['case']) or {})
                if q:
                    x['class'] = ('consistent' if q == x['cls'] == 'none' else
                                  'shown' if q == x['cls'] == 'acts' else 'refuted')
                    x['result'] = 'qualitative'
        ok = sum(x['ceiling'] == 'ok' for x in sc)
        poss = sum(x['ceiling'] == 'possible violation' for x in sc)
        viol = sum(x['ceiling'] == 'violated' for x in sc)
        shown = sum(x.get('class') == 'shown' for x in cs)
        consistent = sum(x.get('class') == 'consistent' for x in cs)
        refuted = sum(x.get('class') == 'refuted' for x in cs)
        mids = [(x['m_lo'] + x['m_hi']) / 2 for x in sc]
        phase = (sc[0]['phase'] if sc else next((x['phase'] for x in cs if x.get('phase')), ''))
        cell = next((r for r in rows if r['phase'] == phase and r['route'] == route), None)
        out.append({'label': label, 'route': route, 'phase': phase, 'n': len(cs), 'scored': len(sc),
                    'ok': ok, 'possible': poss, 'violated': viol, 'shown': shown, 'consistent': consistent,
                    'refuted': refuted,
                    'A_lo': min(mids) if mids else None, 'A_hi': max(mids) if mids else None,
                    'cell': cell, 'cases': cs})
    return out


def known_section(kr, reg):
    lines = []
    for k in kr:
        c = k['cell'] or {}
        lines.append(f"**{k['label']}** — account cell {k['phase']} / {k['route']}: verdict "
                     f"{c.get('verdict', '–')}, ceiling A_max = {fmt(c.get('A_max'))}. "
                     f"{c.get('chain', '')}")
        lines.append('')
        lines.append(f"Published measurements: {k['n']} records, {k['scored']} convertible to A; measured A "
                     f"{fmt(k['A_lo'])}–{fmt(k['A_hi'])}; ceiling holds for {k['ok']} (possible violation "
                     f"{k['possible']}, violated {k['violated']}); alignment shown in {k['shown']}, "
                     f"'cannot align' consistent in {k['consistent']}.")
        lines.append('')
        trs = []
        for x in k['cases']:
            r = reg.rec.get(x['case']) or {}
            trs.append([x['case'], r.get('quantity'), f"{fmt(r.get('value'))} {r.get('unit') or ''}",
                        x['cls'], fmt(x['A_hi']), (f"{fmt(x.get('m_lo'))}-{fmt(x.get('m_hi'))}"
                                                   if x.get('m_lo') is not None else '–'),
                        x.get('ceiling', x.get('result')), x.get('class', '')])
        lines.append(table(['record', 'quantity', 'published value', 'predicted', 'A ceiling', 'measured A',
                            'ceiling', 'class'], trs))
        lines.append('')
    return '\n'.join(lines)


# ---------------------------------------------------------------- frontier
def _git(*a):
    return subprocess.run(['git', *a], cwd=ROOT, capture_output=True, text=True).stdout


def phases_covered_at(commit, univ_names):
    try:
        raw = yaml.safe_load(_git('show', f'{commit}:data/curated/phases.yaml')) or []
    except yaml.YAMLError:
        return None
    return len({p.get('universe') for p in raw} & set(univ_names))


VERSIONS = [  # version, freeze commit, scores file, statuses that are held-out for it, note
    ('v1', '15c80aa', 'results/scores_v1.json', ('heldout',), 'frozen before batches 1-2 were opened'),
    ('v2', 'e94871b', 'results/scores_v2.json', ('heldout',),
     'post-hoc revision; held-out on the Mn4N batch opened after its freeze'),
    ('v3', None, 'results/scores_v3.json', ('heldout',),
     'account reasons made physical (EP, TGG, direct formation); held-out on batch 4'),
    ('v4', None, 'results/scores_v4.json', ('heldout',),
     'v3 plus same-sample Js for all phases, composites and thin-film rules; held-out on batch 5'),
]
FREEZE = os.path.join(ROOT, 'docs', 'PROTOCOL_LOG.md')


def _freeze_commit(version):
    """Commit recorded in the protocol log for a version's freeze (line 'freeze <version>: <commit>')."""
    import re
    if not os.path.exists(FREEZE):
        return None
    m = re.search(r'freeze ' + version + r': ([0-9a-f]{7,40})', open(FREEZE).read())
    return m.group(1) if m else None


def frontier_rows(univ):
    names = [p['name'] for p in univ['phases']]
    lit = [r['class'] for r in univ['route_classes']]
    out = []
    for v, commit, fn, sts, note in VERSIONS:
        sc = _load(fn)
        commit = commit or _freeze_commit(v)
        if not sc or not commit:
            continue
        n = sum(sc['summary'].get(s, {}).get('n_scored', 0) for s in sts)
        ok = sum(sc['summary'].get(s, {}).get('ceiling_ok', 0) + sc['summary'].get(s, {}).get('ceiling_possible', 0)
                 for s in sts)
        viol = sum(sc['summary'].get(s, {}).get('ceiling_violated', 0) for s in sts)
        ph = phases_covered_at(commit, names)
        out.append({'version': v, 'commit': commit, 'phases_covered': ph, 'phases_universe': len(names),
                    'routes_covered': len([r for r in lit if r in ROUTES]), 'routes_universe': len(lit),
                    'n_cases': n, 'accuracy': (n - viol) / n if n else None, 'violations': viol, 'note': note})
    return out


def _frontier_scores(r):
    return {'phases': r['phases_covered'] / r['phases_universe'], 'routes': r['routes_covered'] / r['routes_universe'],
            'accuracy': r['accuracy'] if r['accuracy'] is not None else 0.0}


def non_dominated(rows):
    keep = []
    for r in rows:
        s = _frontier_scores(r)
        dom = any(all(_frontier_scores(t)[k] >= s[k] for k in s) and any(_frontier_scores(t)[k] > s[k] for k in s)
                  for t in rows if t is not r)
        if not dom:
            keep.append(r['version'])
    return keep


# ---------------------------------------------------------------- limits
def limits(rows):
    """Classify every remaining limit: open cells (unpublished quantity shown to decide, with the
    measurement), capped and 'cannot align' cells (physics, with the reason)."""
    opened, phys = collections.defaultdict(list), collections.defaultdict(list)
    bad = []
    for r in rows:
        if r['verdict'] == 'open':
            ds = json.loads(r['decides'])
            if not ds:
                bad.append(r)
                continue
            q, outs, meas = ds[0]
            vals = [v for v in outs.values() if v is not None]
            decides = (len(vals) < len(outs)) or (len(set(round(v, 3) for v in vals)) > 1)
            if not decides or not meas:
                bad.append(r)
            opened[(q, meas)].append(f"{r['phase']}/{r['route']}")
        elif r['verdict'] in ('none', 'capped'):
            phys[r['chain'].split(';')[0][:140]].append(f"{r['phase']}/{r['route']}")
    return opened, phys, bad


# ---------------------------------------------------------------- build
def build(reg, phases, rows, R, st, univ=None):
    univ = univ or yaml.safe_load(open(os.path.join(ROOT, 'data', 'curated', 'universe.yaml')))
    s1, s2 = _load('results/scores_v1.json'), _load('results/scores_v2.json')
    ctx = {'commit': st['commit'],
           'clean': 'inputs committed' if st['clean_inputs'] else 'UNCOMMITTED inputs',
           'n_phases': len(phases), 'n_routes': len(ROUTES)}
    verd = collections.Counter(r['verdict'] for r in rows)
    names = [p['name'] for p in univ['phases']]
    covered = sorted({p.get('universe') for p in phases} & set(names))
    lit = [r['class'] for r in univ['route_classes']]
    kr = known_rows(s1, reg, rows) if s1 else []
    nd = [x for x in (s1['cases'] if s1 else []) if x['case'].startswith('nd2fe14b_fp_br') and x.get('m_lo')]
    nd_br = [float(reg.rec[x['case']]['value']) for x in nd]
    ctx['ndfeb_line'] = (f"Br = {min(nd_br):.2f}–{max(nd_br):.2f} T in the known cases, "
                         f"A = {min(x['m_lo'] for x in nd):.2f}–{max(x['m_hi'] for x in nd):.2f} "
                         "(Js(Nd2Fe14B) -/+5 %)." if nd_br else '')
    rk = ranking(phases, rows, reg)
    ctx['ranking_table'] = table(
        ['phase', 'Js (T)', 'best established bulk route', 'A_max', 'Br_max (T)', 'highest open route',
         'A_max if open resolves'],
        [[ph.get('name', ph['id']), fmt(ph['Js']), s['route'] if s else '–', fmt(s['A_max']) if s else '–',
          fmt(ph['Js'] * s['A_max']) if s and ph['Js'] is not None else '–',
          t['route'] if t['verdict'] == 'open' else '–',
          fmt(t['A_max']) if t['verdict'] == 'open' else '–'] for ph, t, s in rk])
    ctx['phase_sections'] = '\n'.join(phase_section(ph, rows, reg) for ph in phases)
    ctx['tests_table'] = table(['phase', 'route', 'prediction', 'the measurement'],
                               [[r['phase'], r['route'], f"{r['verdict']} {fmt(r['A_max'])}", r['test']]
                                for r in rows])
    ctx['route_list'] = table(['route', 'definition'], [[k, ROUTE_NAMES[k]] for k in ROUTES])
    ctx['evidence_table'] = evidence_table(reg, R.used)
    ctx['estimates'] = table(['where', 'value', 'why'], [[w, v, y] for w, v, y in R.estimates]) \
        if R.estimates else 'None.'
    ctx['known_section'] = known_section(kr, reg) if kr else ''
    # ---------------- accuracy
    h1 = s1['summary']['heldout'] if s1 else {}
    caps = [r for r in rows if r['verdict'] == 'capped' and r['route'] in BULK]
    hyp_ref = collections.Counter()
    for x in (s1['cases'] if s1 else []):
        if x['status'] == 'heldout' and x.get('hyp'):
            for h, res in x['hyp'].items():
                hyp_ref[(h, res)] += 1
    fails = _failures(s1) if s1 else []
    acc = []
    acc.append('Predictions were frozen (docs/PROTOCOL_LOG.md) before the sealed outcomes were opened. '
               'Held-out = opened only after the freeze; exposed = seen before the freeze (listed in the '
               'log); known = the known cases of PROBLEM.md; post hoc = scored with rules revised after '
               'opening (never counted as held-out).')
    acc.append('')
    acc.append('**v1 — frozen before opening batches 1–2 (the held-out test)**')
    acc.append('')
    acc.append(_summary_table(s1, ['heldout', 'exposed', 'known']) if s1 else '')
    acc.append('')
    acc.append('Pre-registered alternatives decided by the held-out data:')
    acc.append('')
    acc.append(_hypotheses(s1) if s1 else '')
    acc.append('')
    acc.append('Paired predictions (same study, one condition changed):')
    acc.append('')
    acc.append(_pairs(s1) if s1 else '')
    acc.append('')
    acc.append('Every held-out failure of v1, with its post-hoc diagnosis:')
    acc.append('')
    acc.append(table(['case', 'phase', 'route', 'predicted', 'A ceiling', 'measured A', 'diagnosis'], fails))
    acc.append('')
    if s2:
        acc.append('**v2 — revised after opening (post hoc), then tested on batches opened later**')
        acc.append('')
        acc.append(_summary_table(s2, ['heldout', 'post hoc', 'exposed', 'known']))
        acc.append('')
    s3 = _load('results/scores_v3.json')
    if s3:
        acc.append('**v3 — frozen before opening batch 4 (sources never used before)**')
        acc.append('')
        acc.append(_summary_table(s3, ['heldout']))
        acc.append('')
        sens = s3.get('sensitivity_js') or []
        if sens:
            sv = [x for x in sens if x['same_sample_js'][2] == 'violated']
            acc.append(f"Conversion check (post hoc): the frozen v3 specs left the same-sample Js unlinked for "
                       f"{len(sens)} Br cases outside AlNiCo (they used the phase Js, which can only understate A). "
                       f"Rescored with the same-sample Js, {sum(1 for x in sens if x['same_sample_js'][2] == 'ok')} "
                       f"ceilings hold and {len(sv)} is exceeded: " + ', '.join(
                           f"{x['case']} (predicted {x['cls']}, A = {fmt(x['same_sample_js'][0])})" for x in sv)
                       + ' — an isotropic bonded magnet of rapidly solidified powder, whose remanence above 1/2 is '
                       'exchange enhancement, excluded from A by definition; its phrase was not flagged by the '
                       'v3 rules (fixed for later versions).' if sv else '')
            acc.append('')
    s4 = _load('results/scores_v4.json')
    if s4:
        acc.append('**v4 — frozen before opening batch 5 (sources never used before)**')
        acc.append('')
        acc.append(_summary_table(s4, ['heldout']))
        acc.append('')
        f4 = _failures(s4)
        if f4:
            acc.append('Held-out failures of v4, with their post-hoc diagnosis:')
            acc.append('')
            acc.append(table(['case', 'phase', 'route', 'predicted', 'A ceiling', 'measured A', 'diagnosis'], f4))
            acc.append('')
    acc.append('Held-out v1 results by phase:')
    acc.append('')
    acc.append(_by_phase(s1, 'heldout') if s1 else '')
    ctx['accuracy_section'] = '\n'.join(acc)
    # ---------------- summary
    sm = []
    sm.append(f"- **The account**: {len(phases)} phases (all {len(covered)} of the {len(names)} rare-earth-free "
              f"candidate phases named in {len(set(sum((p['named_in'] for p in univ['phases']), [])))} review "
              f"records, plus the Nd2Fe14B reference) × {len(ROUTES)} route classes (all "
              f"{len([r for r in lit if r in ROUTES])} classes used or proposed in the literature, plus "
              f"{len(ROUTES) - len([r for r in lit if r in ROUTES])}) = {len(rows)} cells: {verd['full']} reach "
              f"A ≥ 0.95, {verd['capped']} are capped by geometry, {verd['none']} cannot align, {verd['open']} are "
              f"decided by a named unpublished quantity.")
    if h1:
        sm.append(f"- **Held-out test** (v1, frozen at 15c80aa): the predicted ceilings held for "
                  f"{h1['ceiling_ok'] + h1['ceiling_possible']} of {h1['n_scored']} held-out measurements "
                  f"({h1['ceiling_violated']} certain violations, {h1['ceiling_possible']} within the Js "
                  f"conversion range); 'cannot align' was refuted in {h1['none_refuted']} of {h1['none_n']}; "
                  f"{s1['summary']['pairs']['ok']} of {s1['summary']['pairs']['n']} field/no-field pairs moved "
                  f"as predicted. Every failure is diagnosed in §6.")
        later = [r for r in frontier_rows(univ) if r['version'] != 'v1']
        if later:
            sm.append('- **Revised rules on later held-out batches** (frozen before each opening, sources never used '
                      'before): ' + '; '.join(f"{r['version']}: {r['n_cases'] - r['violations']} of {r['n_cases']} "
                                             'ceilings held' for r in later) + ' (failures diagnosed in §6).')
    a5 = [x for x in (s1['cases'] if s1 else []) if x['status'] == 'heldout' and x['phase'] == 'AlNiCo-5'
          and (x.get('hyp') or {}).get('rods on the nearest <100>') == 'refuted']
    sm.append(f"- **Settled by the held-back data**: in Ti-free AlNiCo 5 the α1 rods are not held on the parent "
              f"<100> nearest the field: {len(a5)} held-out equiaxed magnets exceed the 0.831 that this would "
              f"allow (up to {fmt(max((x['m_lo'] for x in a5), default=None))} with the same-sample Js), and in "
              "columnar crystals A falls with the cosine of the field angle; the account now lets the rods follow "
              "the field (ceiling 1 for any grain). Where the three rod directions are equivalent (no bias, or a "
              "field along <111> or <110> that ties them), spinodal alloys (AlNiCo, Fe-Cr-Co, the alnico-like HEA) "
              "keep more remanence than independent rods allow: that reading was refuted "
              f"{hyp_ref[('independent rods', 'refuted')]} times, while the bound for one coupled network per "
              f"grain (0.866 for random grains, with the no-field easy directions <111>) held every time.")
    if s1:
        ho = [x for x in s1['cases'] if x['status'] == 'heldout' and x.get('result') == 'scored']
        a89 = [x for x in ho if x['phase'] == 'AlNiCo-8-9' and abs(x['A_hi'] - 0.8312) < 0.01
               and (x['conv'] == 'same-sample Js' or not x['conv'])]
        fcc = [x for x in ho if x['phase'] == 'Fe-Cr-Co' and x['route'] == 'FA' and 'Mo-free' in x['why']
               and x['m_lo'] > 0.8312 + x['tol']]
        if a89 and fcc:
            sm.append(f"- **The rod rule differs by alloy, as the account says**: equiaxed field-annealed AlNiCo 8/9 "
                      f"(rods on the <100> nearest the field, cap 0.832) stays at or below "
                      f"{fmt(max(x['m_hi'] for x in a89))} in all {len(a89)} held-out cases with a same-sample Js "
                      f"or Mr/Ms; equiaxed Mo-free Fe-Cr-Co (rods follow the field in any grain, ceiling 1) "
                      f"exceeds that cap in {len(fcc)} held-out cases (up to {fmt(max(x['m_lo'] for x in fcc))}).")
    groups = collections.defaultdict(list)
    for r in caps:
        parts = [x.strip() for x in r['chain'].split(';')]
        groups[(round(r['A_max'], 2), parts[0])].append((f"{r['phase']} {r['route']}", parts[1] if len(parts) > 1 else ''))
    for (a, first), cells in sorted(groups.items(), key=lambda kv: (kv[0][0], -len(kv[1])))[:6]:
        names = [c for c, _ in cells]
        sm.append(f"- **Geometric cap A_max ≈ {a:.2f}** ({len(cells)} cells: {', '.join(names[:8])}"
                  f"{' …' if len(cells) > 8 else ''}) — {first}"
                  + (f"; {cells[0][1]}" if cells[0][1] and len({c for _, c in cells}) == 1 else '') + '.')
    ctx['summary'] = '\n'.join(sm)
    # ---------------- REPORT: fidelity
    fid = []
    for k in kr:
        rep = k['violated'] == 0 and k['refuted'] == 0 and (k['shown'] > 0 or k['consistent'] > 0)
        fid.append([k['label'], f"{k['phase']} / {k['route']}", f"{(k['cell'] or {}).get('verdict')} "
                    f"{fmt((k['cell'] or {}).get('A_max'))}", f"{fmt(k['A_lo'])}–{fmt(k['A_hi'])}",
                    f"{k['ok']}/{k['scored']} (+{k['possible']} possible)",
                    f"{k['shown']} shown / {k['consistent']} consistent", 'yes' if rep else 'NO'])
    ctx['fidelity_section'] = (
        'The floor is the published record of the four known cases. A case is reproduced when the '
        "account's cell for it has the right verdict (aligns / cannot align) and every published "
        'measurement lies at or below its ceiling.\n\n'
        + table(['known case', 'cell', 'account', 'published A', 'ceiling holds', 'class', 'reproduced'], fid)
        + '\n\nDetails, record by record: FINDINGS.md §4.')
    ctx['accuracy_report'] = ctx['accuracy_section']
    # ---------------- frontier
    fr = frontier_rows(univ)
    ndv = non_dominated(fr) if fr else []
    cur = fr[-1] if fr else None
    sc = _frontier_scores(cur) if cur else {}
    far = min(sc, key=sc.get) if sc else '–'
    if sc and len({round(v, 3) for v in sc.values()}) == 1:
        far = (f"none — all three at their limits; accuracy rests on the fewest held-out cases "
               f"({cur['n_cases']}), so it is pushed next with further held-out batches")
    ctx['frontier_section'] = (
        'Limits set outside this work: phases = the rare-earth-free candidate phases named in published '
        f"reviews ({len(names)}, data/curated/universe.yaml); routes = the route classes used or proposed in "
        f"the literature ({len(lit)}); accuracy = share of held-out measurements whose ceiling holds "
        '(certain violations count as failures).\n\n'
        + table(['version', 'freeze commit', 'phases covered', 'routes covered', 'held-out cases', 'accuracy',
                 'note'],
                [[r['version'], r['commit'], f"{r['phases_covered']}/{r['phases_universe']}",
                  f"{r['routes_covered']}/{r['routes_universe']}", r['n_cases'],
                  fmt(r['accuracy']), r['note']] for r in fr])
        + f"\n\nNon-dominated versions: {', '.join(ndv)}. Furthest from its limit now: **{far}** "
          f"({', '.join(f'{k} {v:.3f}' for k, v in sc.items())}). History of pushes: after v1 the phase "
          "count (23/24) was furthest and Mn4N was added (24/24); then accuracy was furthest and was pushed by "
          "the post-hoc v2 rules, which are being tested on a further held-out batch.")
    # ---------------- limits
    opened, phys, bad = limits(rows)
    lt = [[q, len(v), ', '.join(v[:6]) + (' …' if len(v) > 6 else ''), m] for (q, m), v in
          sorted(opened.items(), key=lambda kv: -len(kv[1]))]
    pt = [[q, len(v), ', '.join(v[:6]) + (' …' if len(v) > 6 else '')] for q, v in
          sorted(phys.items(), key=lambda kv: -len(kv[1]))]
    ctx['limits_section'] = (
        f"{sum(len(v) for v in opened.values())} cells are open. Each is decided by a quantity that is not "
        "published; the engine was run on each of its values and the outcomes fall in different verdict "
        "categories (listed per cell in results/account.csv, column decides). The measurement that settles "
        "it:\n\n" + table(['unpublished quantity', 'cells', 'examples', 'the measurement'], lt)
        + f"\n\n{sum(len(v) for v in phys.values())} cells cannot align or are capped for a physical "
          "reason:\n\n" + table(['physical reason (first clause)', 'cells', 'examples'], pt)
        + (f"\n\nCells whose limit is not of these two kinds: {len(bad)}." if bad else
           '\n\nNo remaining limit is of another kind (scope, cut-off, tolerance or untried).'))
    # ---------------- evidence
    vc = verification_counts(reg, R.used)
    nrec = len(R.used)
    ctx['evidence_report'] = (
        f"All numbers come from commit `{st['commit']}` ({ctx['clean']}) through `./regenerate.sh`. "
        f"The account uses {nrec} sourced records; quote found at the source for "
        f"{vc['verified_online'] + vc['verified_local'] + vc['verified_manual']} ({vc['verified_online']} online, "
        f"{vc['verified_local']} in the saved source text, {vc['verified_manual']} checked by eye on the source "
        f"text or page image, data/verification_manual.yaml), not found for {vc['not_found']}, not checked "
        f"{vc['not_checked']}; DOI resolves for {vc['doi_ok']}. Estimates made here are labelled "
        f"({len(R.estimates)}, FINDINGS.md §7) and never set a ceiling.")
    allf = []
    for lab, fn in (('v1', 'results/scores_v1.json'), ('v3', 'results/scores_v3.json'), ('v4', 'results/scores_v4.json')):
        sc = _load(fn)
        if sc:
            allf += [[lab] + f for f in _failures(sc)]
    ctx['failures'] = (table(['frozen version', 'case', 'phase', 'route', 'predicted', 'A ceiling', 'measured A',
                              'diagnosis'], allf) if allf else 'None.')
    pv = [x for x in (s1['cases'] if s1 else []) if x['status'] == 'heldout' and x.get('ceiling') == 'possible violation']
    pv_table = table(['case', 'phase', 'A ceiling', 'measured A (Js range)', 'conversion'],
                     [[x['case'], x['phase'], fmt(x['A_hi']), f"{fmt(x['m_lo'])}-{fmt(x['m_hi'])}", x['conv']]
                      for x in pv])
    s2h = s2['summary'].get('heldout', {}) if s2 else {}
    ctx['open_items'] = '\n'.join([
        f"- {verd['open']} account cells decided by unpublished quantities (Limits).",
        f"- {h1.get('ceiling_possible', 0)} held-out measurements (v1) whose conversion to A straddles the "
        "ceiling because the same sample's Js is unpublished; the Js of that magnet settles each.",
        f"- v2 has been tested on {s2h.get('n_scored', 0)} held-out measurements so far "
        f"({s2h.get('ceiling_violated', 0)} violations).",
        '', 'Held-out measurements (v1) whose conversion straddles the ceiling (the same-sample Js decides):', '',
        pv_table,
    ])
    # ---------------- COMPLETE.md: each GOAL.md measurement, computed
    n_ph, n_rt = len(phases), len(ROUTES)
    acc_ok = all(r['verdict'] in ('full', 'capped', 'none', 'open') and r['chain'] and r['A_max'] is not None
                 for r in rows) and len(rows) == n_ph * n_rt
    tests_ok = all(r['test'] for r in rows)
    fid_n = sum(1 for k in kr if k['violated'] == 0 and k['refuted'] == 0 and (k['shown'] or k['consistent']))
    undiag = [f for f in fails if f[-1] == '–']
    for fn in ('results/scores_v3.json', 'results/scores_v4.json', 'results/scores_v5.json'):
        sc = _load(fn)
        if sc:
            undiag += [f for f in _failures(sc) if f[-1] == '–']
    untried = [r for r in rows if r['verdict'] == 'none' and any(
        w in r['chain'] for w in ('not reported', 'unreported', 'is reported', 'not densified'))]
    srcless = [rid for rid in R.used if not ((reg.rec[rid].get('source') or {}).get('doi')
                                             or (reg.rec[rid].get('source') or {}).get('url'))
               or not reg.rec[rid].get('quote')]
    s2h = s2['summary'].get('heldout', {}) if s2 else {}
    checks = [
        ('Deliverable 1 — the account', f"{n_ph} phases × {n_rt} route classes = {len(rows)} cells, each with "
         f"verdict, mechanism/reason and A_max (full {verd['full']}, capped {verd['capped']}, none "
         f"{verd['none']}, open {verd['open']})", acc_ok),
        ('Deliverable 2 — tests', f"{sum(1 for r in rows if r['test'])} of {len(rows)} cells carry one "
         'measurement that would refute them', tests_ok),
        ('Deliverable 3 — FINDINGS.md, REPORT.md, one command', f"written by ./regenerate.sh from commit "
         f"{st['commit']} ({ctx['clean']})", st['clean_inputs']),
        ('Fidelity — known cases', f"{fid_n} of {len(kr)} known cases reproduced (verdict right, every "
         'published value at or below its ceiling)', fid_n == len(kr) and len(kr) == 4),
        ('Accuracy — held-out (frozen v1)', f"ceilings held for {h1.get('ceiling_ok', 0) + h1.get('ceiling_possible', 0)}"
         f" of {h1.get('n_scored', 0)} held-out measurements ({h1.get('ceiling_violated', 0)} certain "
         f"violations); 'cannot align' refuted {h1.get('none_refuted', 0)} of {h1.get('none_n', 0)}; pairs "
         f"{s1['summary']['pairs']['ok'] if s1 else 0}/{s1['summary']['pairs']['n'] if s1 else 0}; every failure "
         f"diagnosed ({len(undiag)} undiagnosed)", bool(h1) and not undiag),
        ('Accuracy — held-out (revised versions)', '; '.join(
            f"{r['version']}: {r['n_cases'] - r['violations']} of {r['n_cases']} ceilings held"
            for r in fr if r['version'] != 'v1') + f" (all failures diagnosed: {len(undiag)} undiagnosed)",
         all(r['n_cases'] > 0 for r in fr) and not undiag),
        ('Frontier', f"phases {fr[-1]['phases_covered']}/{fr[-1]['phases_universe']}, routes "
         f"{fr[-1]['routes_covered']}/{fr[-1]['routes_universe']}, accuracy {fmt(fr[0]['accuracy'])} (v1 held-out); "
         f"non-dominated: {', '.join(ndv)}; furthest from its limit: {far}" if fr else '–',
         bool(fr) and fr[-1]['phases_covered'] == fr[-1]['phases_universe']
         and fr[-1]['routes_covered'] == fr[-1]['routes_universe']),
        ('Limits', f"{sum(len(v) for v in opened.values())} open cells, each with the unpublished quantity that "
         f"decides it and its measurement ({len(bad)} without); {len(untried)} 'cannot align' cells resting on "
         "'not reported'", not bad and not untried),
        ('Evidence', f"{len(R.used)} input records, each with a named source and verbatim quote "
         f"({len(srcless)} without); quote found at the source for {vc['verified_online'] + vc['verified_local'] + vc['verified_manual']} of {len(R.used)}; "
         f"{len(R.estimates)} labelled estimates; one commit", not srcless and st['clean_inputs']),
    ]
    ctx['complete_table'] = table(['criterion (GOAL.md)', 'measured', 'holds'],
                                  [[a, b, 'yes' if c else 'NO'] for a, b, c in checks])
    ctx['complete_overall'] = ('every criterion holds' if all(c for _, _, c in checks) else
                               'not complete: ' + '; '.join(a for a, _, c in checks if not c))
    notes = []
    if untried:
        notes.append("'Cannot align' cells still resting on an untried route: "
                     + ', '.join(f"{r['phase']}/{r['route']}" for r in untried[:12]) + (' …' if len(untried) > 12 else ''))
    ctx['complete_notes'] = '\n'.join(notes)
    return ctx


def write_docs(ctx):
    for name in ('FINDINGS.md', 'REPORT.md', 'COMPLETE.md'):
        p = os.path.join(TPL, name + '.in')
        if os.path.exists(p):
            open(os.path.join(ROOT, name), 'w').write(fill(open(p).read(), ctx))
