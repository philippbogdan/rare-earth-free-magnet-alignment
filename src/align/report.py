"""Builds FINDINGS.md and REPORT.md from computed results (every number generated)."""
import json
import os

from .data import ROOT
from .render import fmt, table, fill
from .routes import ROUTE_NAMES, ROUTES

TPL = os.path.join(ROOT, 'docs', 'templates')
BULK = [r for r in ROUTES if r != 'EP']


def best_routes(rows, phase):
    cells = [r for r in rows if r['phase'] == phase and r['route'] in BULK]
    top = max(cells, key=lambda r: r['A_max'])
    return top, cells


def ranking(phases, rows, reg):
    out = []
    for ph in phases:
        top, cells = best_routes(rows, ph['id'])
        sure = [r for r in cells if r['verdict'] in ('full', 'capped')]
        sure_top = max(sure, key=lambda r: r['A_max']) if sure else None
        out.append((ph, top, sure_top))
    out.sort(key=lambda t: -(t[0]['Js'] * (t[2]['A_max'] if t[2] else 0.5)))
    return out


def phase_section(ph, rows, reg):
    lines = [f"### {ph.get('name', ph['id'])}", '']
    facts = []
    for key, label, unit in (('Js', 'Js (300 K)', 'T'), ('K1', 'K1', 'MJ/m³'), ('Tc', 'Tc', 'K')):
        if ph.get(key) is not None:
            facts.append(f'{label} = {fmt(ph[key], 3 if key != "Tc" else 0)} {unit}')
    lines.append('; '.join(facts) + '.')
    lines.append('')
    rws = [r for r in rows if r['phase'] == ph['id']]
    lines.append(table(['route', 'verdict', 'A_max', 'Br_max (T)', 'mechanism'],
                       [[r['route'], r['verdict'], fmt(r['A_max']), fmt(r['Br_max_T']), r['chain']]
                        for r in rws]))
    lines.append('')
    return '\n'.join(lines)


def evidence_table(reg, used):
    rows = []
    for rid in sorted(used):
        r = reg.rec[rid]
        src = r.get('source') or {}
        rows.append([rid, fmt(r.get('value')), r.get('unit', ''), reg.cite(rid),
                     src.get('evidence', ''), 'yes' if reg.verified(rid) else 'not yet'])
    return table(['record', 'value', 'unit', 'source', 'evidence', 'quote verified'], rows)


def build(reg, phases, rows, R, st):
    ctx = {'commit': st['commit'],
           'clean': 'inputs committed' if st['clean_inputs'] else 'UNCOMMITTED inputs',
           'n_phases': len(phases), 'n_routes': len(ROUTES)}
    rk = ranking(phases, rows, reg)
    ctx['ranking_table'] = table(
        ['phase', 'Js (T)', 'best established route', 'A_max', 'Br_max (T)', 'highest open route', 'A_max if open resolves'],
        [[ph.get('name', ph['id']), fmt(ph['Js']), s['route'] if s else '–', fmt(s['A_max']) if s else '–',
          fmt(ph['Js'] * s['A_max']) if s else '–', t['route'] if t['verdict'] == 'open' else '–',
          fmt(t['A_max']) if t['verdict'] == 'open' else '–'] for ph, t, s in rk])
    ctx['phase_sections'] = '\n'.join(phase_section(ph, rows, reg) for ph in phases)
    ctx['tests_table'] = table(['phase', 'route', 'prediction', 'the measurement'],
                               [[r['phase'], r['route'], f"{r['verdict']} {fmt(r['A_max'])}", r['test']]
                                for r in rows])
    ctx['route_list'] = table(['route', 'definition'], [[k, ROUTE_NAMES[k]] for k in ROUTES])
    ctx['evidence_table'] = evidence_table(reg, R.used)
    ctx['estimates'] = table(['where', 'value', 'why'], [[w, v, y] for w, v, y in R.estimates]) \
        if R.estimates else 'None.'
    return ctx


def write_docs(ctx):
    for name in ('FINDINGS.md', 'REPORT.md'):
        p = os.path.join(TPL, name + '.in')
        if os.path.exists(p):
            open(os.path.join(ROOT, name), 'w').write(fill(open(p).read(), ctx))
