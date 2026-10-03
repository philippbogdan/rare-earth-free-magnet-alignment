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
    if not acts or A <= ISO + 1e-3:
        return 'none'
    return 'full' if A >= 0.95 - 1e-9 else 'capped'


RANK = {'full': 3, 'capped': 2, 'open': 1, 'none': 0}


def _best_path_value(cells):
    """Highest ceiling; on ties a settled verdict is preferred over an open one."""
    return max(cells, key=lambda c: (round(c.A_max, 3), RANK.get(c.verdict, 0)))


CANDIDATE_FIBRES = {'hexagonal': ([0, 0, 0, 1], [1, 0, -1, 0], [1, 1, -2, 0]),
                    'tetragonal': ([0, 0, 1], [1, 0, 0], [1, 1, 0], [1, 0, 1], [1, 1, 1]),
                    'cubic': ([0, 0, 1], [1, 1, 0], [1, 1, 1]),
                    'orthorhombic': ([1, 0, 0], [0, 1, 0], [0, 0, 1])}


def _fibre_outcomes(path, on_parent=True):
    """Ceiling for each low-index candidate fibre (equal variant weights): used when the
    fibre/growth direction is unpublished, to show whether it decides the answer."""
    L = path.LP if on_parent else path.LH
    out = {}
    for d in CANDIDATE_FIBRES[L.system]:
        par = E.texture_fibre(L.direction(d))
        if on_parent:
            A, _ = _equal_best(path, par)
        else:
            A, _ = _equal_best(Path({'structure': path.ph['structure'],
                                      'easy_axis_vec': path.ph['easy_axis_vec']}, {'mode': 'direct'}), par)
        out['<' + ''.join(str(x) for x in d) + '> fibre'] = A
    return out


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


# ---------------------------------------------------------------- MO
SHAPE_A = {   # ceiling of perfect shape alignment, by particle shape and easy-axis relation
    ('plate', 'normal'): (1.0, 'press direction (plate normals align with it)'),
    ('plate', 'in-plane'): (2 / np.pi, 'any direction in the pressing plane'),
    ('rod', 'axis'): (1.0, 'flow/extrusion direction (rod axes align with it)'),
    ('rod', 'perpendicular'): (2 / np.pi, 'perpendicular to the flow direction'),
    ('equiaxed', None): (ISO, ''),
}


def route_MO(ph):
    h = (ph.get('particles') or {}).get('habit')
    if not h:
        outs = {'plates, easy axis normal': 1.0, 'rods, easy axis along rod': 1.0,
                'equiaxed': ISO}
        return Cell(ph['id'], 'MO', 'open', 1.0, '',
                    chain='shape of single-crystal particles and its relation to the easy axis unpublished',
                    decides=[('particle habit (shape vs easy axis)', outs,
                              'SEM + EBSD of single-crystal particles')], A_alt=outs)
    if h.get('single_crystal') is False:
        return Cell(ph['id'], 'MO', 'none', ISO, '', chain='particles are not single crystals: '
                    + h.get('reason', ''), refs=h.get('refs', []))
    key = (h['shape'], h.get('easy') if h['shape'] != 'equiaxed' else None)
    A, ax = SHAPE_A[key]
    if h.get('single_crystal') is None and A > ISO + 1e-3:
        outs = {f"single-crystal {h['shape']}s": A, 'polycrystalline particles': ISO}
        return Cell(ph['id'], 'MO', 'open', A, ax,
                    chain=f"{h.get('reason', '')}; whether such particles are single crystals is unpublished",
                    decides=[('single-crystal anisometric particles', outs, 'EBSD of crushed particles')],
                    A_alt=outs, refs=h.get('refs', []))
    return Cell(ph['id'], 'MO', _verdict(A, A > ISO + 1e-3), A, ax,
                chain=f"single-crystal {h['shape']}s with the easy axis {h.get('easy', '')} to the shape; "
                      'pressing/shear flow orients the shape', refs=h.get('refs', []))


# ---------------------------------------------------------------- texture sources (DS/AM)
def _texture_route(ph, route, key):
    out = []
    for p in ph['paths']:
        path = Path(ph, p)
        if path.mode == 'random_nucleation':
            out.append(Cell(ph['id'], route, 'none', ISO, '',
                            chain='product nucleates without an orientation relationship: parent texture not inherited'))
            continue
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
        p0 = ph['paths'][0]
        path = Path(ph, p0)
        outs = _fibre_outcomes(path, on_parent=path.mode != 'direct')
        if path.mode == 'massive':
            outs = {k: v for k, v in outs.items()}
        return Cell(ph['id'], route, 'open', max(outs.values()), '',
                    chain=f'{route} growth texture of the ' + ('parent' if path.mode != 'direct' else
                                                               'hard phase') + ' unpublished',
                    decides=[(f'{route} growth direction', outs,
                              f'EBSD pole figure of a {route} sample')], A_alt=outs)
    return _best_path_value(out)


def route_DS(ph):
    return _texture_route(ph, 'DS', 'ds_texture')


def route_AM(ph):
    return _texture_route(ph, 'AM', 'am_texture')


# ---------------------------------------------------------------- biased formation
def _field_status(ph, path):
    """(level1, level2, why, extra_decides) for a field bias.

    level = True (acts), False (cannot), None (decided by an unpublished quantity).
    Single-step paths: both levels equal. Two-step paths: level2 = product choice within an
    intermediate variant (product, or ferromagnetic intermediate, ordered at T_f); level1 =
    choice of the intermediate variant (intermediate ferromagnetic at T_f)."""
    p = path.p
    mag, why = formation_magnetic(ph, path)
    acts2 = mag in ('yes', 'partly')
    extra = []
    if not acts2:
        A_label = 'forms below Tc'
        if p.get('product_ferromagnetic', True):
            extra.append((f'formation of the phase below its Curie point '
                          f'({p.get("Tc_forming", ph["Tc"]):.0f} K) at a practical rate',
                          f'isothermal anneal of the parent ~50 K below Tc: phase fraction vs time (XRD)'))
        x10 = para_bound(ph, path, FIELD_HIGH_T)
        if x10 >= 0.1:
            x1 = para_bound(ph, path, FIELD_PRACTICAL_T)
            why += (f'; paramagnetic bias <= {x1:.2g} kT at 2 T, up to {x10:.2g} kT at 10 T '
                    'only if dchi ~ chi (upper bound)')
            extra.append(('paramagnetic susceptibility anisotropy dchi/chi of the forming phase at T_f',
                          'torque magnetometry above Tc on a single crystal, or formation in 10 T '
                          'followed by a pole figure'))
        else:
            why += f'; paramagnetic bias < {x10:.2g} kT even at 10 T'
    if not path.two_step:
        lv = True if acts2 else (None if extra else False)
        return lv, lv, why, extra
    it = p['intermediate']
    l2 = True if acts2 else (None if extra else False)
    if it.get('ferromagnetic') is False:
        l1, w1 = False, 'intermediate not ferromagnetic: its variant choice cannot be biased'
    elif it.get('Tc') is not None:
        lo, hi = p['T_f']
        l1 = bool(it['Tc'] > lo)
        w1 = f"intermediate Tc {it['Tc']:.0f} K vs formation from {lo:.0f} K"
    else:
        l1, w1 = None, 'intermediate ferromagnetic but its Curie point/easy axis are unpublished'
        extra = extra + [('Curie point and easy axis of the intermediate ' + it.get('name', ''),
                          'M(T) and torque on the intermediate made by low-temperature ordering')]
    return l1, l2, why + '; ' + w1, extra


def _stress_status(ph, path):
    p = path.p
    if p['mode'] == 'massive':
        return False, False, 'incoherent (massive) growth: no coherent transformation strain', []
    if p['mode'] == 'direct':
        return False, False, 'no parent lattice: no transformation strain', []
    if p.get('strain_known', True) is False:
        return None, None, 'transformation strain of the variants unpublished', [
            ('transformation strain of the variants', 'lattice parameters of parent and product at T_f')]
    if path.two_step:
        it = p['intermediate']
        l1 = it.get('strain', True)
        return l1, True, ('coherent shear/ordering strains at both steps (stress sign chosen to suit)'
                          if l1 else 'coherent strain at the second step only'), []
    return True, True, 'coherent transformation with a tetragonal strain (stress sign chosen to suit)', []


def _parents_for(path, source):
    """List of (label, parents, kind) for a texture source; kind: 'random'|'fibre'|'single'.
    Several entries mean the source texture is unpublished (candidates)."""
    par = path.parent
    if source == 'random':
        return [('untextured parent', E.texture_random(20000), 'random')], True
    if source == 'sc':
        return [('single-crystal parent', E.texture_single(), 'single')], True
    key = {'ds': 'ds_texture', 'am': 'am_texture', 'rx': 'rx_texture'}[source]
    if source == 'rx' and par.get('rollable') is False:
        return [], False
    t = par.get(key)
    if t:
        if t.get('kind', 'fibre') == 'sheet':
            return [(f"sheet texture {t['nd']}||ND {t['rd']}||RD", E.texture_single(_sheet_g(path, t)),
                     'single')], True
        return [(f"parent fibre {t['dir']}", _fibre(path, t['dir']), 'fibre')], True
    out = []
    for d in CANDIDATE_FIBRES[path.LP.system]:
        out.append((f"<{''.join(str(x) for x in d)}> parent fibre", _fibre(path, d), 'fibre'))
    return out, None


def _sheet_g(path, t):
    nd = path.LP.direction(t['nd'])
    rd = path.LP.direction(t['rd'])
    rd = rd - (rd @ nd) * nd
    rd /= np.linalg.norm(rd)
    td = np.cross(nd, rd)
    return np.stack([rd, td, nd], 0)        # rows: crystal directions along sample x, y, z


def _A_bias(path, parents, kind, tex_kind, l1, l2):
    """Ceiling for given level statuses (True/False) and texture."""
    if not l1 and not l2:
        if tex_kind == 'random':
            return ISO, 'untextured: no alignment'
        if tex_kind == 'single':
            ax, fr, per = path.axes(parents)
            A, _ = E.best_axis(ax, fr, 1500)
            return A, 'best axis'
        return _equal_best(path, parents)
    if tex_kind == 'single':
        ax, fr, per = path.axes(parents)
        best = (-1.0, None)
        for h in E.fibonacci_hemisphere(600):
            w = path.weights(ax, fr, kind, h, l1, l2)
            A = E.alignment(ax, w, h)
            if A > best[0]:
                best = (A, h)
        return best[0], 'bias along the best crystal direction'
    A, al = _select_best(path, parents, kind, +1, fibre=(tex_kind == 'fibre'), level1=l1, level2=l2)
    return A, (f'bias axis at {al} deg to the parent fibre axis' if tex_kind == 'fibre' else 'bias axis')


def _bias_route(ph, route, kind, source):
    cells = []
    for p in ph['paths']:
        path = Path(ph, p)
        if path.mode == 'random_nucleation' and kind == 'stress':
            cells.append(Cell(ph['id'], route, 'none', ISO, '',
                              chain='nucleation without an orientation relationship: no coherent strain to bias'))
            continue
        if path.mode in ('direct', 'random_nucleation'):
            if source != 'random' and path.mode == 'direct':
                continue
            if kind == 'field':
                mag, why = formation_magnetic(ph, path)
                acts = mag in ('yes', 'partly')
                cells.append(Cell(ph['id'], route, 'full' if acts else 'none', 1.0 if acts else ISO,
                                  'field direction' if acts else '',
                                  chain=(f'forms without a parent crystal; {why}' +
                                         ('; nuclei take the orientation the field favours' if acts else ''))))
            continue
        l1, l2, why, extra = (_field_status if kind == 'field' else _stress_status)(ph, path)
        texs, avail = _parents_for(path, source)
        if avail is False:
            cells.append(Cell(ph['id'], route, 'none', ISO, '', chain='parent cannot be rolled or drawn'))
            continue
        # scenarios over unknown levels
        if path.two_step:
            lv1 = [l1] if l1 is not None else [True, False]
            lv2 = [l2] if l2 is not None else [True, False]
            scen = [(a, b) for a in lv1 for b in lv2]
        else:
            scen = [(l1, l2)] if l1 is not None else [(True, True), (False, False)]
        outcomes, axis_of = {}, {}
        for lab, par, tk in texs:
            for a, b in scen:
                if True:
                    if path.mode == 'massive' and not (a or b):
                        A_or, ax_ = _A_bias(path, par, kind, tk, False, False)
                        A, ax_ = max(A_or, ISO), ax_
                    else:
                        A, ax_ = _A_bias(path, par, kind, tk, a, b)
                    key = (lab if len(texs) > 1 else '') + (
                        '' if (l1 is not None and l2 is not None) else
                        f"{' ; ' if len(texs) > 1 else ''}bias {'acts' if (a or b) else 'cannot act'}"
                        + (' (both steps)' if (a and b and path.two_step) else
                           (' (product step only)' if (b and not a and path.two_step) else
                            (' (intermediate step only)' if (a and not b and path.two_step) else ''))))
                    outcomes[key or 'value'] = A
                    axis_of[key or 'value'] = ax_
        tex_label = texs[0][0] if len(texs) == 1 else 'parent texture unpublished'
        if len(outcomes) == 1:
            A = list(outcomes.values())[0]
            acts = bool(l1 or l2)
            cells.append(Cell(ph['id'], route, _verdict(A, A > ISO + 1e-3), A, list(axis_of.values())[0],
                              chain=f'{tex_label}; {why}',
                              prereq=[(kind + ' bias acts at formation', acts, why)]))
            continue
        dec = []
        if l1 is None or l2 is None:
            for q, meas in extra:
                dec.append((q, outcomes, meas))
        if len(texs) > 1:
            dec.append((f'{source.upper()} texture of the parent', outcomes, 'EBSD pole figure of the parent as made'))
        cells.append(Cell(ph['id'], route, 'open', max(outcomes.values()), '',
                          chain=f'{tex_label}; {why}', decides=dec, A_alt=outcomes))
    if not cells:
        return Cell(ph['id'], route, 'none', ISO, '', chain='no formation path to bias with this texture source')
    return _best_path_value(cells)


def route_FA(ph):
    return _bias_route(ph, 'FA', 'field', 'random')


def route_DSFA(ph):
    return _bias_route(ph, 'DS+FA', 'field', 'ds')


def route_AMFA(ph):
    return _bias_route(ph, 'AM+FA', 'field', 'am')


def route_RXFA(ph):
    return _bias_route(ph, 'RX+FA', 'field', 'rx')


def route_SCFA(ph):
    return _bias_route(ph, 'SC+FA', 'field', 'sc')


def route_SA(ph):
    return _bias_route(ph, 'SA', 'stress', 'random')


def route_DSSA(ph):
    return _bias_route(ph, 'DS+SA', 'stress', 'ds')


def route_RXSA(ph):
    return _bias_route(ph, 'RX+SA', 'stress', 'rx')


def route_SCSA(ph):
    return _bias_route(ph, 'SC+SA', 'stress', 'sc')


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
    if fs.get('primary') == 'unknown':
        A_hard = 1.0 if ph['structure']['laue'] != 'm-3m' else ISO
        outs = {'hard phase crystallises as free primary crystals': A_hard, 'it does not': ISO}
        return Cell(ph['id'], 'FS', 'open', A_hard, '', chain=fs.get('reason', ''),
                    decides=[('primary crystallisation of the hard phase from a melt (liquidus)', outs,
                              'DTA/quench series of an off-stoichiometric melt; then solidify in >= 10 T and '
                              'measure the pole figure')], A_alt=outs, refs=fs.get('refs', []))
    return Cell(ph['id'], 'FS', 'none', ISO, '',
                chain=fs.get('reason', 'neither the hard phase nor its parent crystallises from a melt '
                                       'as free crystals'))


# ---------------------------------------------------------------- HD, EP, TGG
def route_HD(ph):
    hd = ph.get('hd')
    if not hd:
        path = Path(ph, {'mode': 'direct'})
        outs = _fibre_outcomes(path, on_parent=False)
        return Cell(ph['id'], 'HD', 'open', max(outs.values()), '',
                    chain='deformation texture of the phase unpublished',
                    decides=[('fibre formed by hot deformation', outs,
                              'pole figure after hot deformation')], A_alt=outs)
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
        n = int(ep.get('in_plane_variants') or 1)
        if n > 1:      # n easy axes in the film plane, 180/n degrees apart: best in-plane direction
            ang = np.radians(np.arange(n) * 180.0 / n)
            A = max(float(np.mean(np.abs(np.cos(ang - t)))) for t in np.radians(np.linspace(0, 180, 721)))
            return Cell(ph['id'], 'EP', _verdict(A), A, 'best in-plane direction',
                        chain=f'epitaxial film with {n} in-plane easy-axis variants (not a bulk magnet)',
                        refs=ep.get('refs', []))
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


FUN = {'FP': route_FP, 'MO': route_MO, 'SC+FA': route_SCFA, 'SC+SA': route_SCSA, 'FS': route_FS, 'DS': route_DS, 'AM': route_AM, 'FA': route_FA,
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
    c.A_max = max(outs)            # for an open cell: the ceiling if the quantity falls its way
    if len(cats) == 1:
        c.verdict = cats.pop()
    return c


def _hypotheses(ph):
    """Expand formation paths whose orientation relationship is unpublished into two
    hypotheses: random nucleation (nothing inherited) and the stated analogue OR."""
    import copy
    unk = [i for i, p in enumerate(ph['paths']) if p.get('OR') == 'unknown']
    if not unk:
        return None
    a, b = copy.deepcopy(ph), copy.deepcopy(ph)
    for i in unk:
        a['paths'][i]['mode'] = 'random_nucleation'
        b['paths'][i]['OR'] = b['paths'][i]['OR_analogue']
        b['paths'][i]['mode'] = b['paths'][i].get('mode_analogue', 'OR')
        if 'intermediate_analogue' in b['paths'][i]:
            b['paths'][i]['intermediate'] = b['paths'][i]['intermediate_analogue']
    return a, b, [ph['paths'][i] for i in unk]


def evaluate(ph):
    hyp = _hypotheses(ph)
    if not hyp:
        return [_settle(FUN[r](ph)) for r in ROUTES]
    a, b, unk = hyp
    out = []
    lab_a = 'no orientation relationship (random nucleation)'
    lab_b = 'orientation relationship as in ' + unk[0].get('OR_analogue_name', 'the analogue')
    for r in ROUTES:
        ca, cb = _settle(FUN[r](a)), _settle(FUN[r](b))
        if abs(ca.A_max - cb.A_max) < 0.005 and ca.verdict == cb.verdict:
            ca.chain += ' (same under either orientation relationship hypothesis)'
            out.append(ca)
            continue
        outs = {lab_a: ca.A_max, lab_b: cb.A_max}
        c = Cell(ph['id'], r, 'open', max(outs.values()), '',
                 chain=f'orientation relationship of {unk[0]["parent"].get("name", "parent")} -> product '
                       f'unpublished. If none: {ca.chain}. If as analogue: {cb.chain}',
                 decides=[('orientation relationship of the transformation', outs,
                           'EBSD/TEM of a partially transformed parent crystal')] + ca.decides + cb.decides,
                 A_alt=outs)
        out.append(_settle(c))
    return out


# ---------------------------------------------------------------- the single measurement
DELTA = 0.03     # stated margin of a measurement of A (Mr/Ms or pole-figure average)


def test_for(c, ph):
    """One measurement per prediction, with the outcome that would refute it."""
    A = c.A_max
    name = ph.get('name', ph['id'])
    if c.verdict == 'open':
        d = c.decides[0] if c.decides else None
        if d:
            outs = ', '.join(f'{k}: A_max = {v:.3f}' for k, v in d[1].items() if v is not None)
            return f'Measure {d[0]} by {d[2]}; it selects the ceiling ({outs}).'
        return 'Unpublished input: see limits.'
    r = c.route
    if r == 'FP':
        if c.verdict == 'full':
            return (f'Mr/Ms of field-aligned, single-variant {name} particles fixed in resin, measured '
                    f'along the alignment axis at 300 K after saturation (demagnetisation-corrected): '
                    f'predicted to reach >= {A - DELTA - 0.02:.2f}; refuted if, with particles shown '
                    'single-variant by EBSD, it stays below that.')
        return (f'Mr/Ms of field-aligned {name} particles of the stated internal structure: refuted if '
                f'> {A + DELTA:.2f}.')
    if r in ('DS', 'AM'):
        return (f'Easy-axis pole figure (EBSD or XRD) of {name} formed from a {r} parent of measured '
                f'fibre texture, converted to A along the best axis: refuted if A > {A + DELTA:.2f}'
                + ('' if c.verdict == 'full' else ' (the cap is geometric: no bias acts)') + '.')
    if r in ('FA', 'DS+FA', 'AM+FA', 'RX+FA'):
        if c.verdict == 'none':
            return (f'Form {name} from the same parent at the same temperature with and without a '
                    f'saturating field and compare A (Mr/Ms or pole figure): refuted if the field '
                    f'raises A by more than {DELTA:.2f}.')
        return (f'Form {name} in a saturating field below its Curie point and measure the variant '
                f'fractions (XRD/EBSD) or Mr/Ms along the field: predicted ceiling {A:.3f}; refuted if '
                'the favoured-variant fraction stays at the no-field value.')
    if r in ('SA', 'DS+SA', 'RX+SA'):
        if c.verdict == 'none':
            return (f'Form {name} under a uniaxial stress near the yield stress and without: refuted if '
                    f'the variant fractions (A) differ by more than {DELTA:.2f}.')
        return (f'Form {name} under uniaxial stress of the favoured sign and measure variant fractions: '
                f'predicted ceiling {A:.3f}; refuted if no selection occurs.')
    if r == 'FS':
        if c.verdict == 'none':
            return (f'Solidify in 0 and >= 10 T and compare the pole figures of {name}: refuted if A '
                    f'differs by more than {DELTA:.2f}.')
        return (f'Solidify in a field (ferromagnetic crystals: >= 1 T; paramagnetic: >= 10 T) and '
                f'measure the pole figure of the {name} crystals: predicted A up to {A:.3f}.')
    if r == 'HD':
        return (f'Pole figure of {name} after hot deformation along the deformation axis: predicted '
                f'ceiling {A:.3f}' + (f'; refuted if A > {A + DELTA:.2f}' if c.verdict != 'full' else '') + '.')
    if r == 'HDDR':
        return ('A of powder after a disproportionation-recombination cycle: refuted if A exceeds '
                f'{A + DELTA:.2f}.' if c.verdict == 'none' else f'Predicted ceiling {A:.3f}.')
    if r == 'EP':
        return f'XRD pole figure / Mr/Ms of an epitaxial {name} film: predicted A >= {A - DELTA:.2f}.'
    if r == 'TGG':
        return (f'Pole figure after templated grain growth: predicted up to {A:.3f}.' if c.verdict != 'none'
                else f'Refuted if templated growth gives A > {0.5 + DELTA:.2f} in {name}.')
    return ''


def evaluate_with_tests(ph):
    cells = evaluate(ph)
    for c in cells:
        c.test = test_for(c, ph)
    return cells
