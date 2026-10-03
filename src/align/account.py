"""The account: every route evaluated for every phase with the same rules.

Verdicts: 'full' (a mechanism acts and A_max >= 0.95), 'capped' (a mechanism acts but the
geometry caps A_max below 0.95), 'none' (no alignment mechanism acts: A_max is the value
without one), 'open' (decided by an unpublished quantity: the ceiling for each outcome is
given, with the measurement that settles it).
"""
from dataclasses import dataclass, field

import numpy as np

from . import engine as E
from .routes import (Path, _equal_best, _fibre, _select_best, formation_magnetic, para_bound,
                     FIELD_PRACTICAL_T, FIELD_HIGH_T, RX_TEXTURE, ROUTES)

ISO = 0.5


@dataclass
class Cell:
    phase: str
    route: str
    verdict: str
    A_max: float
    axis: str = ''
    chain: str = ''
    prereq: list = field(default_factory=list)       # (name, ok, reason)
    decides: list = field(default_factory=list)      # (quantity, {outcome: A}, measurement)
    test: str = ''
    refs: list = field(default_factory=list)
    A_alt: dict = field(default_factory=dict)

    def label(self):
        return self.verdict


def _verdict(A, acts=True):
    if not acts:
        return 'none'
    return 'full' if A >= 0.95 - 1e-9 else 'capped'


def _best_path_value(cells):
    return max(cells, key=lambda c: c.A_max)


# ---------------------------------------------------------------- FP
def route_FP(ph):
    sp = ph['particles']
    status = sp['single_variant']
    refs = sp.get('refs', [])
    if status == 'yes':
        return Cell(ph['id'], 'FP', 'full', 1.0, 'field direction',
                    chain='single-variant particles each turn their easy axis to the field',
                    prereq=[('single-variant particles', True, sp['reason'])],
                    test=('Mr/Ms along the alignment axis (saturating field, demagnetisation '
                          'corrected) of field-aligned particles fixed in resin; refuted if '
                          '< 0.90 with particles verified single-crystalline by EBSD'), refs=refs)
    A_poly = _particle_internal(ph, sp)
    if status == 'no':
        return Cell(ph['id'], 'FP', _verdict(A_poly, A_poly > ISO + 1e-3), A_poly,
                    'field direction',
                    chain=('particles contain several variants/grains; each particle turns only its '
                           'net anisotropy axis to the field: ' + sp['reason']),
                    prereq=[('single-variant particles', False, sp['reason'])],
                    test=('Mr/Ms of field-aligned particles fixed in resin; the account predicts '
                          f'<= {A_poly:.2f}'), refs=refs)
    return Cell(ph['id'], 'FP', 'open', 1.0, 'field direction',
                chain='depends on whether single-variant particles can be made',
                prereq=[('single-variant particles', None, sp['reason'])],
                decides=[('single-variant particles of this phase', {'yes': 1.0, 'no': A_poly},
                          sp.get('measurement', 'EBSD of particles of the size needed'))],
                test=f'Mr/Ms of aligned particles: 1.0 if single-variant, {A_poly:.2f} if not',
                refs=refs, A_alt={'yes': 1.0, 'no': A_poly})


def _particle_internal(ph, sp):
    """A of a free body whose interior is the variant set of one parent grain (equal
    weights) or a random polycrystal."""
    kind = sp.get('internal', 'random')
    if kind == 'random':
        return ISO
    path = Path(ph, ph['paths'][sp.get('path', 0)])
    ax, fr, per = path.axes(E.texture_single())
    return E.body_alignment(ax, fr)


# ---------------------------------------------------------------- texture sources (DS/AM)
def _texture_route(ph, route, key):
    out = []
    for p in ph['paths']:
        path = Path(ph, p)
        src = (p.get('parent') or ph).get(key) if path.mode != 'direct' else ph.get(key)
        if src is None:
            continue
        d = src['dir']
        if path.mode == 'direct':
            A, ax = _equal_best(path, _fibre(path, d))
            out.append(Cell(ph['id'], route, _verdict(A), A, ax,
                            chain=f'hard phase grows from the melt with {d} along the growth/build '
                                  f'axis; easy axis at fixed angle to it',
                            refs=src.get('refs', [])))
            continue
        A_or, ax = _equal_best(path, _fibre(path, d))
        if path.mode == 'massive':
            out.append(Cell(ph['id'], route, 'capped' if A_or > ISO + 1e-3 else 'none', A_or, ax,
                            chain=(f'parent fibre {d}; massive (incoherent) growth: product keeps at '
                                   f'most the OR of its nucleation grain -> A between {ISO:.2f} '
                                   f'(random) and {A_or:.3f} (OR kept); no bias, equal variants'),
                            refs=src.get('refs', []) + p.get('refs', []),
                            A_alt={'OR kept': A_or, 'random': ISO}))
        else:
            out.append(Cell(ph['id'], route, _verdict(A_or), A_or, ax,
                            chain=(f'parent fibre {d}; {len(path.V)} variants with equal weights'),
                            refs=src.get('refs', []) + p.get('refs', [])))
    if not out:
        return Cell(ph['id'], route, 'open', 1.0, '', chain=f'{route} texture of the parent unpublished',
                    decides=[(f'{route} growth texture of the parent/hard phase', {}, 'EBSD of a '
                              f'{route} sample')])
    return _best_path_value(out)


def route_DS(ph):
    return _texture_route(ph, 'DS', 'ds_texture')


def route_AM(ph):
    return _texture_route(ph, 'AM', 'am_texture')


# ---------------------------------------------------------------- biased formation
def _bias_acts(ph, p, kind):
    path = Path(ph, p)
    if kind == 'field':
        mag, why = formation_magnetic(ph, path)
        if mag in ('yes', 'partly'):
            return True, why, None
        x1 = para_bound(ph, path, FIELD_PRACTICAL_T)
        x10 = para_bound(ph, path, FIELD_HIGH_T)
        if x10 < 0.1:
            return False, why + f'; paramagnetic bias < {x10:.2g} kT even at 10 T', None
        return False, why + (f'; paramagnetic bias <= {x1:.2g} kT at 2 T but up to {x10:.2g} kT at '
                             '10 T if dchi ~ chi'), x10
    st = p.get('stress_sign', 0)
    if p['mode'] == 'massive' and not p.get('stress_selectable', False):
        return False, 'incoherent (massive) growth: no transformation strain to couple to', None
    if not st:
        return None, 'transformation strain (sign/size) unpublished', None
    return True, f'coherent transformation strain (sign {st:+d})', None


def _biased(ph, route, kind, texture_key=None):
    cells = []
    for p in ph['paths']:
        path = Path(ph, p)
        acts, why, x10 = _bias_acts(ph, p, kind)
        sign = p.get('stress_sign', +1) if kind == 'stress' else +1
        if texture_key is None:
            parents, fib = E.texture_random(20000), False
            tex = 'untextured parent'
        elif texture_key == 'rx':
            fam = (p.get('parent') or {}).get('family')
            if fam not in RX_TEXTURE or not (p.get('parent') or {}).get('rollable'):
                cells.append(Cell(ph['id'], route, 'none', ISO, '',
                                  chain='parent cannot be rolled/drawn and recrystallised to a sheet texture'))
                continue
            nd, rd = RX_TEXTURE[fam]
            parents, fib = _fibre(path, nd), True
            tex = f'{fam} rolling/recrystallisation texture ({nd} || ND)'
        else:
            src = (p.get('parent') or {}).get(texture_key)
            if path.mode == 'direct' or src is None:
                continue
            parents, fib = _fibre(path, src['dir']), True
            tex = f"parent fibre {src['dir']}"
        if path.mode == 'direct':
            if acts:
                cells.append(Cell(ph['id'], route, 'full', 1.0, 'field direction',
                                  chain=f'nuclei form without an orientation constraint; {why}'))
            else:
                cells.append(Cell(ph['id'], route, 'none', ISO, '', chain=why))
            continue
        A_eq, _ = _equal_best(path, parents) if fib else (ISO, '')
        if acts is True:
            A, al = _select_best(path, parents, kind, sign, fibre=fib)
            ax = (f'bias axis at {al} deg to the parent fibre axis' if fib else 'bias axis')
            if path.mode == 'massive' and kind == 'field':
                chain = (f'{tex}; {why}; complete selection of the nucleus variant nearest the field '
                         f'({len(path.V)} variants)')
            else:
                chain = f'{tex}; {why}; complete selection among {len(path.V)} variants'
            cells.append(Cell(ph['id'], route, _verdict(A), A, ax, chain=chain,
                              prereq=[(kind + ' bias acts at formation', True, why)],
                              A_alt={'no bias': A_eq}))
        elif acts is None:
            A, al = _select_best(path, parents, kind, +1, fibre=fib)
            A2, al2 = _select_best(path, parents, kind, -1, fibre=fib)
            cells.append(Cell(ph['id'], route, 'open', max(A, A2), '',
                              chain=f'{tex}; {why}',
                              decides=[('transformation strain of the variant', {
                                  'c elongated (tension selects)': A, 'c shortened (compression selects)': A2},
                                  'lattice parameters of parent and product at T_f')]))
        else:
            dec = []
            A_hi = None
            if x10 is not None:
                A_hi, _ = _select_best(path, parents, kind, +1, fibre=fib)
                dec = [('paramagnetic susceptibility anisotropy dchi/chi of the forming phase at T_f',
                        {'dchi/chi >~ 0.1 at 10 T': A_hi, 'small': A_eq},
                        'torque magnetometry above Tc on a single crystal, or formation in 10 T '
                        'followed by a pole figure')]
            cells.append(Cell(ph['id'], route, 'open' if dec else 'none', A_eq,
                              '', chain=f'{tex}; {why}', prereq=[(kind + ' bias acts at formation', False, why)],
                              decides=dec, A_alt={'10 T, large dchi': A_hi} if A_hi else {}))
    if not cells:
        return Cell(ph['id'], route, 'none', ISO, '', chain='no formation path to bias')
    return _best_path_value(cells)


def route_FA(ph):
    return _biased(ph, 'FA', 'field')


def route_DSFA(ph):
    return _biased(ph, 'DS+FA', 'field', 'ds_texture')


def route_AMFA(ph):
    return _biased(ph, 'AM+FA', 'field', 'am_texture')


def route_RXFA(ph):
    return _biased(ph, 'RX+FA', 'field', 'rx')


def route_SA(ph):
    return _biased(ph, 'SA', 'stress')


def route_DSSA(ph):
    return _biased(ph, 'DS+SA', 'stress', 'ds_texture')


def route_RXSA(ph):
    return _biased(ph, 'RX+SA', 'stress', 'rx')


# ---------------------------------------------------------------- FS
def route_FS(ph):
    fs = ph.get('melt', {})
    if fs.get('primary') == 'hard':
        frac = fs.get('fraction_max', 1.0)
        return Cell(ph['id'], 'FS', 'full', 1.0, 'field direction',
                    chain=('hard-phase crystals grow in the melt and turn their easy axis to the field '
                           f'(non-cubic: anisotropic susceptibility or magnetisation); volume fraction '
                           f'of these crystals <= {frac}'), refs=fs.get('refs', []))
    if fs.get('primary') == 'parent':
        p = ph['paths'][fs.get('path', 0)]
        path = Path(ph, p)
        if p['parent'].get('family') in ('fcc', 'bcc') or p['parent']['structure']['laue'] == 'm-3m':
            return Cell(ph['id'], 'FS', 'none', ISO, '',
                        chain='parent crystals are cubic: isotropic susceptibility, no torque in a field')
        out = {}
        for name, d in (('parent c parallel to field', [0, 0, 0, 1] if path.LP.system == 'hexagonal'
                         else [0, 0, 1]),):
            par = _fibre(path, d)
            A_eq, ax = _equal_best(path, par)
            out[name] = A_eq if path.mode != 'massive' else max(A_eq, ISO)
        # c perpendicular to field (negative dchi): parent c-axes in the plane normal to H
        par_perp = E.texture_fibre(path.LP.direction([1, 0, -1, 0] if path.LP.system == 'hexagonal'
                                                     else [1, 0, 0]))
        A_perp, _ = _equal_best(path, par_perp)
        out['parent c perpendicular to field'] = A_perp
        A = max(out.values())
        return Cell(ph['id'], 'FS', 'open', A, '',
                    chain=('non-cubic parent crystals align in the melt (sign of their susceptibility '
                           'anisotropy decides the axis), then transform with equal variant weights'),
                    decides=[('sign of dchi of the parent', out, 'texture of the parent solidified in a '
                              'field (quenched before transforming)')], A_alt=out)
    return Cell(ph['id'], 'FS', 'none', ISO, '',
                chain=fs.get('reason', 'neither the hard phase nor its parent crystallises from a melt '
                                       'as free crystals'))


# ---------------------------------------------------------------- HD, EP, TGG
def route_HD(ph):
    hd = ph.get('hd')
    if not hd:
        return Cell(ph['id'], 'HD', 'open', 1.0, '', chain='deformation texture of the phase unpublished',
                    decides=[('fibre formed by hot deformation', {}, 'pole figure after hot deformation')])
    if hd.get('possible') is False:
        return Cell(ph['id'], 'HD', 'none', ISO, '', chain=hd['reason'], refs=hd.get('refs', []))
    path = Path(ph, {'mode': 'direct'})
    A, ax = _equal_best(path, E.texture_fibre(path.LH.direction(hd['fibre'])))
    return Cell(ph['id'], 'HD', _verdict(A), A, ax + ' (deformation axis)',
                chain=hd['reason'], refs=hd.get('refs', []))


def route_HDDR(ph):
    h = ph.get('hddr')
    if not h:
        return Cell(ph['id'], 'HDDR', 'none', ISO, '',
                    chain='no reversible disproportionation-recombination reaction with an orientation '
                          'memory is reported for this phase')
    return Cell(ph['id'], 'HDDR', _verdict(h['A']), h['A'], 'memory axis', chain=h['reason'],
                refs=h.get('refs', []))


def route_EP(ph):
    ep = ph.get('epitaxy', {})
    if ep.get('demonstrated'):
        return Cell(ph['id'], 'EP', 'full', 1.0, 'film normal or in-plane axis',
                    chain='single-crystal film on a lattice-matched substrate (not a bulk magnet)',
                    refs=ep.get('refs', []))
    return Cell(ph['id'], 'EP', 'open', 1.0, '', chain='epitaxial growth not reported',
                decides=[('epitaxial growth', {}, 'XRD of a film on a matched substrate')])


def route_TGG(ph):
    t = ph.get('tgg', {})
    if t.get('applicable'):
        return Cell(ph['id'], 'TGG', 'full', 1.0, 'template axis', chain=t['reason'],
                    refs=t.get('refs', []))
    return Cell(ph['id'], 'TGG', 'none', ISO, '', chain=t.get('reason', 'not densified by sintering '
                                                                      'with grain growth'))


FUN = {'FP': route_FP, 'FS': route_FS, 'DS': route_DS, 'AM': route_AM, 'FA': route_FA,
       'DS+FA': route_DSFA, 'AM+FA': route_AMFA, 'RX+FA': route_RXFA, 'SA': route_SA,
       'DS+SA': route_DSSA, 'RX+SA': route_RXSA, 'HD': route_HD, 'HDDR': route_HDDR, 'EP': route_EP,
       'TGG': route_TGG}


def _settle(c):
    """'open' only if the unpublished quantity changes the verdict category; otherwise the
    larger outcome is the ceiling (the quantity stays listed for information)."""
    if c.verdict != 'open' or not c.decides:
        return c
    outs = [v for d in c.decides for v in d[1].values() if v is not None]
    if not outs:
        return c
    cats = {_verdict(v, v > ISO + 1e-3) for v in outs}
    if len(cats) == 1:
        c.A_max = max(outs)
        c.verdict = cats.pop()
    return c


def evaluate(ph):
    return [_settle(FUN[r](ph)) for r in ROUTES]
