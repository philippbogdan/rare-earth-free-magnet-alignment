"""Conditional predictions for the sealed published cases, built only from their
value-free conditions (data/sealed/*.conditions.yaml). Writes data/curated/cases.yaml.

Each case gets:
  phase, route          the account phase and route class it tests
  A_hi                  the ceiling of A under the case's own conditions (engine)
  A_base                A without the route's alignment mechanism (engine)
  cls                   'acts' if A_hi > A_base + DELTA, 'none' otherwise; 'open' when an
                        unpublished quantity decides; 'info' when not a test of alignment
  meas                  how the measured value converts to A (MrMs | Br | vfrac | class | skip)
  status                'heldout' or 'exposed' (docs/PROTOCOL_LOG.md)
"""
import os
import re

import numpy as np
import yaml

from . import engine as E
from .account import evaluate, ISO, DELTA
from .curate import load_phases
from .data import ROOT
from .routes import Path, _fibre

PLANAR = 2 / np.pi                       # easy axes random in a plane, measured in that plane

PHASE_MAP = {                            # collector phase label -> account phase id
    'tau-MnAl': 'tau-MnAl', 'tau-MnAl-C': 'tau-MnAl-C', 'tau-MnAl(C)': 'tau-MnAl-C',
    'LTP-MnBi': 'LTP-MnBi', 'L10-FeNi': 'L10-FeNi', 'L10-FePt': 'L10-FePt',
    'L10-CoPt': 'L10-CoPt', 'L10-FePd': 'L10-FePd', "α''-Fe16N2": 'Fe16N2',
    "α'-Fe8N": None, 'D022-Mn3Ga': 'D022-Mn3Ga', 'D022-Mn3Ge': 'D022-Mn3Ge',
    'L10-MnGa': 'L10-MnGa', 'Zr2Co11': 'Zr2Co11', 'HfCo7': 'HfCo7',
    'hcp-Co': 'Co-nanorods', 'hcp-CoNi': 'Co-nanorods', 'Co nanorods (shape anisotropy)': 'Co-nanorods',
    'Co nanorods/nanowires (shape anisotropy)': 'Co-nanorods', 'epsilon-Fe2O3': 'eps-Fe2O3',
    'Co3C': 'Co3C', 'bct-FeCo': 'bct-FeCo', 'Fe5SiB2': 'Fe5SiB2', 'Fe5PB2': 'Fe5PB2',
    '(Fe,Co)5PB2': 'Fe5PB2', 'BaFe12O19': 'SrM-hexaferrite', 'SrFe12O19': 'SrM-hexaferrite',
    'M-type ferrite (unspecified)': 'SrM-hexaferrite', 'M-type hexaferrite': 'SrM-hexaferrite',
    '(Sr,La)(Fe,Co)12O19': 'SrM-hexaferrite', '(Ca,La)(Fe,Co)12O19': 'SrM-hexaferrite',
    '(La,Na)(Fe,Co)12O19': 'SrM-hexaferrite', 'Nd2Fe14B': 'Nd2Fe14B',
}

EXPOSED = {   # docs/PROTOCOL_LOG.md
    'prefix': ['fepd_scfa_', 'uni_out_rt_fepd_scfa', 'uni_out_fepd_l10_fa', 'mnal_hd_', 'mnal_sc_sa_kubo1976'],
    'ids': ['mnal_am_tiden2024_ref_mrms', 'mnal_ep_takata2018_sto', 'mnal_ep_takata2018_mgo_mn4n',
            'mnal_ep_takata2018_mgo', 'mnal_ep_takata2018_rocking', 'mnbi_other_tang2022_fieldsinter300',
            'mnbi_other_tang2022_fieldsinter250'],
}


def exposed(cid):
    return cid in EXPOSED['ids'] or any(cid.startswith(p) for p in EXPOSED['prefix'])


def along_kind(s):
    s = (s or '').lower()
    if any(w in s for w in ('perpendicular', 'perp', 'hard axis', 'hard direction', 'axial (compression')):
        return 'perp'
    if any(w in s for w in ('parallel', 'alignment', 'easy axis', 'field direction', 'preferred',
                            'orientation direction', 'needle long', 'wire axis', 'c-axis', 'c axis',
                            'ribbon length', 'radial', 'diameter', 'rolling direction', 'pressing',
                            'growth', 'film normal', 'out-of-plane', 'out of plane', 'extrusion',
                            'compression axis', 'texture plane')):
        return 'par'
    if any(w in s for w in ('any', 'random', 'n/a', 'not applicable', 'arbitrary')):
        return 'any'
    return 'unknown'


def meas_kind(c):
    q, u = str(c.get('quantity')), str(c.get('unit') or '')
    if q == 'Mr_over_Ms':
        return 'MrMs'
    if q in ('[redacted]', 'Br') and u.strip() in ('T', 'kG'):
        return 'Br'
    if q == 'texture_metric' and u.startswith('variant_fraction'):
        return 'vfrac'
    if q == 'texture_metric' and u.strip() == 'A_reported':
        return 'A'
    if q == 'texture_metric':
        return 'class'
    return 'skip'


def saturating(c):
    m = (str(c.get('measurement') or '') + ' ' + str(c.get('process') or '')).lower()
    return not ('not saturat' in m or 'non-saturat' in m or 'far below saturation' in m)


class Ctx:
    def __init__(self):
        self.phases, _ = load_phases()
        self.ph = {p['id']: p for p in self.phases}
        self.cells = {}

    def cell(self, pid, route):
        if pid not in self.cells:
            self.cells[pid] = {c.route: c for c in evaluate(self.ph[pid])}
        return self.cells[pid][route]


def pred(A_hi, A_base, why, cls=None):
    if cls is None:
        cls = 'acts' if A_hi > A_base + DELTA else 'none'
    return {'A_hi': round(float(A_hi), 4), 'A_base': round(float(A_base), 4), 'cls': cls, 'why': why}


def ceiling_for(ctx, pid, route, along):
    """Generic: the account cell's ceiling along the alignment axis; a measurement across
    the alignment axis of an aligned sample is bounded by the isotropic value."""
    c = ctx.cell(pid, route)
    base = ISO
    if along == 'perp':
        return pred(ISO, ISO, f'{route}: measured across the alignment axis (bounded by the random value)',
                    cls='none')
    cls = 'open' if c.verdict == 'open' else None
    return pred(c.A_max, base, f'{route} cell: {c.verdict}; {c.chain[:160]}', cls=cls)


# ---------------------------------------------------------------- MnAl specific geometry
def mnal_sc_stress(ctx, pid, theta1, theta2, axis='best'):
    """Single-crystal ε compressed along (theta1 from [0001], theta2 azimuth from [1-100]):
    complete selection of the variant(s) with c most perpendicular to the compression axis
    (shear mechanism: compression along c destabilises a variant), A along the best axis."""
    ph = ctx.ph[pid]
    p = [q for q in ph['paths'] if q['mode'] == 'OR'][0]
    path = Path(ph, p)
    t1, t2 = np.radians(theta1), np.radians(theta2)
    c_e = path.LP.direction([0, 0, 0, 1])
    m_e = path.LP.direction([1, -1, 0, 0])
    o_e = np.cross(c_e, m_e)
    s = np.cos(t1) * c_e + np.sin(t1) * (np.cos(t2) * m_e + np.sin(t2) * o_e)
    ax, fr, per = path.axes(E.texture_single())
    c2 = (ax[0] @ s) ** 2
    keep = (c2 <= c2.min() + 1e-6).astype(float)
    w = (keep / keep.sum())[None]
    A, _ = E.best_axis(ax, w, 3000)
    return A, int(keep.sum())


def mnal_ds_stress(ctx, pid, measure):
    """Columnar ε (<0001> fibre along y) compressed along x (perpendicular to growth):
    variants with c most perpendicular to x selected; A along x, y (growth) or z."""
    ph = ctx.ph[pid]
    p = [q for q in ph['paths'] if q['mode'] == 'OR'][0]
    path = Path(ph, p)
    par = E.texture_fibre(path.LP.direction([0, 0, 0, 1]), axis=(0, 1, 0))
    ax, fr, per = path.axes(par)
    x = np.array([1.0, 0, 0])
    c2 = np.einsum('nmi,i->nm', ax, x) ** 2
    keep = (c2 <= c2.min(1, keepdims=True) + 1e-6).astype(float)
    w = keep / keep.sum(1, keepdims=True)
    vec = {'x': (1, 0, 0), 'y': (0, 1, 0), 'z': (0, 0, 1)}[measure]
    return E.alignment(ax, w, vec)


def build():
    ctx = Ctx()
    out = []
    files = sorted(f for f in os.listdir(os.path.join(ROOT, 'data', 'sealed')) if f.endswith('.conditions.yaml'))
    for fn in files:
        topic = fn.replace('_outcomes.conditions.yaml', '')
        for c in yaml.safe_load(open(os.path.join(ROOT, 'data', 'sealed', fn))):
            spec = encode(ctx, topic, c)
            if spec:
                spec.update({'case': c['id'], 'topic': topic,
                             'status': 'exposed' if exposed(c['id']) else 'heldout'})
                out.append(spec)
    return out


def encode(ctx, topic, c):
    """One case -> spec; rules by (phase, route) group with per-case conditions."""
    cid, route = c['id'], str(c.get('route'))
    pid = PHASE_MAP.get(c.get('phase'))
    mk = meas_kind(c)
    al = along_kind(c.get('measured_along'))
    cond = c.get('conditions') or {}
    proc = str(c.get('process') or '').lower()
    base = {'phase': pid, 'route': route, 'meas': {'kind': mk, 'along': al, 'saturating': saturating(c)}}
    if pid is None or pid not in ctx.ph:
        return {**base, 'pred': pred(1, ISO, 'phase not in the account', cls='info')}
    if mk == 'skip':
        return {**base, 'pred': pred(1, ISO, 'quantity is not remanence or texture', cls='info')}
    r = route.lower()
    # ------------------------------------------------ isotropic samples (no alignment step)
    if r.startswith('isotropic') or 'no alignment' in r or 'without field' in r and 'cold rolling' not in r:
        if any(w in proc for w in ('film', 'cast', 'spin-coated')) and pid == 'Co-nanorods' and al != 'any':
            return {**base, 'pred': pred(PLANAR, PLANAR, 'rods lie in the film plane at random azimuth: '
                                                         'planar random, in-plane value 2/pi', cls='none')}
        exch = any(w in proc for w in ('milled', 'melt-spun', 'ribbon', 'nanocrystal', 'quenched', 'mechanically alloyed'))
        return {**base, 'pred': pred(ISO, ISO, 'no alignment step: random grains' +
                                     ('; fine-grained, exchange enhancement possible (not alignment)' if exch else ''),
                                     cls='info' if exch else 'none')}
    # ------------------------------------------------ phase-specific groups
    f = globals().get('enc_' + topic)
    if f:
        res = f(ctx, pid, c, route, r, proc, cond, al, mk)
        if res is not None:
            return {**base, **res}
    # ------------------------------------------------ generic: route class of the account
    rc = route.split(' ')[0].replace(':', '')
    if rc in ctx.cells.get(pid, {}) or rc in ('FP', 'FS', 'DS', 'AM', 'FA', 'SA', 'HD', 'EP', 'TGG', 'MO',
                                               'SC+FA', 'SC+SA', 'DS+FA', 'RX+FA', 'DS+SA'):
        try:
            return {**base, 'pred': ceiling_for(ctx, pid, rc, al)}
        except KeyError:
            pass
    return {**base, 'pred': pred(1, ISO, 'route not specific enough to test', cls='info')}


# ---------------------------------------------------------------- topic rules
def enc_mnal(ctx, pid, c, route, r, proc, cond, al, mk):
    cid = c['id']
    if r == 'hd':
        if pid == 'tau-MnAl' or 'below the tau solubility' in proc or 'c 0.47' in proc or 'binary' in proc:
            return {'pred': pred(ISO, ISO, 'carbon below its solubility / binary MnAl: warm working gives no '
                                           'texture (decomposition)', cls='none')}
        if any(w in proc for w in ('upset', 'compressed', 'free hot deformation')):
            if al == 'perp':
                return {'pred': pred(ISO, ISO, 'upsetting: easy axes leave the compression axis; measured '
                                               'along it the remanence is below the random value', cls='none')}
            return {'pred': pred(PLANAR, ISO, 'upsetting gives a planar texture (easy axes in the plane '
                                              'normal to the load): in-plane ceiling 2/pi')}
        if any(w in proc for w in ('rolled', 'swag', 'extru')):
            return {'pred': pred(1.0, ISO, 'radial compression (extrusion, swaging, groove rolling): '
                                           '[001] fibre along the axis')}
        return {'pred': pred(1.0, ISO, 'hot deformation, mode not stated', cls='info')}
    if r.startswith('sc+sa'):
        if 'binary' in proc or pid == 'tau-MnAl':
            return {'pred': pred(0.6133, 0.6133, 'binary MnAl: stress during the transformation gives no '
                                                 'single variant (no carbon); variants equal', cls='none')}
        t1, t2 = cond.get('theta1_deg'), cond.get('theta2_deg')
        if t1 is not None and t2 is not None and mk in ('Br', 'MrMs'):
            A, n = mnal_sc_stress(ctx, pid, float(t1), float(t2))
            return {'pred': pred(A, 0.6133, f'compression at theta1={t1}, theta2={t2} deg in the ε frame '
                                             f'selects {n} variant(s) with c most perpendicular to the load')}
        if 're-compress' in proc or 're-pressed' in proc:
            return {'pred': pred(1.0, 0.6133, 're-compression of a τ crystal re-selects a single variant')}
        return {'pred': pred(1.0, 0.6133, 'single-crystal ε compressed during the transformation (530-830 C)')}
    if r.startswith('ds+sa'):
        meas = 'x' if 'pressing direction' in str(c.get('measured_along')) and 'perpendicular' not in str(c.get('measured_along')) else (
            'y' if 'growt' in str(c.get('measured_along')) else 'z')
        A = mnal_ds_stress(ctx, pid, meas)
        return {'pred': pred(A, 0.6133 if meas == 'y' else ISO,
                             f'columnar ε (<0001> fibre) compressed across the growth axis; A along {meas}')}
    if r.startswith('am'):
        return {'pred': pred(0.6133, ISO, 'printed ε <0001> fibre; τ formed at 833 K (massive mode): at most '
                                          'the OR value 0.613 along the build direction', cls='none')}
    if r.startswith('ds') or 'melt-spun ribbon' in r or r.startswith('other: ds'):
        return {'pred': pred(0.6133, ISO, 'textured ε parent without bias: capped at 0.613', cls='none')}
    if 'sc (single-crystal epsilon parent, no bias)' in r:
        return {'pred': pred(0.6133, 0.6133, 'single-crystal ε without bias: six equal variants', cls='none')}
    if r.startswith('fa'):
        T = cond.get('T_process_K')
        try:
            T = float(str(T).split('-')[-1])
        except (TypeError, ValueError):
            T = None
        Tc = ctx.ph[pid]['Tc']
        if T is not None and T > Tc:
            return {'pred': pred(0.6133 if 'textured' in str(cond.get('parent_texture')) else ISO, ISO,
                                 f'formed at {T:.0f} K above Tc {Tc:.0f} K: the field cannot select', cls='none')}
        return {'pred': pred(0.899, ISO, f'formed below Tc {Tc:.0f} K: selection possible; the ε\' step decides '
                                         'how much (open)', cls='acts')}
    if r.startswith('sa'):
        T = cond.get('T_process_K')
        if T is not None and float(str(T).split('-')[0]) > 766:
            return {'pred': pred(0.6133, ISO, 'stress during ε->τ above the displacive window (massive mode): '
                                              'no selection', cls='none')}
        return {'pred': pred(0.899, ISO, 'stress during the displacive transformation selects variants')}
    if 'hot compaction of epsilon powder' in r:
        return {'pred': pred(ISO, ISO, 'transformation under load at 973 K (massive): no selection', cls='none')}
    if r.startswith('fp'):
        if 'twin-free' in proc or 'single-crystal' in proc:
            return {'pred': pred(1.0, ISO, 'single-variant particles aligned in a field')}
        if 'milled' in proc or 'twinned' in proc or 'polycrystalline' in proc:
            return {'pred': pred(0.6133, ISO, 'twinned/polycrystalline particles: each turns only its net '
                                              'axis (<= 0.613)', cls='none')}
        return {'pred': pred(1.0, ISO, 'field-aligned particles (crystallinity not stated)', cls='info')}
    if r.startswith('ep'):
        return {'pred': pred(1.0, ISO, 'epitaxial film (c normal)')}
    if 'ds + hot deformation' in r:
        return {'pred': pred(0.6133, ISO, 'binary MnAl: DS capped at 0.613, warm working adds no texture', cls='none')}
    return None


def enc_l10(ctx, pid, c, route, r, proc, cond, al, mk):
    ph = ctx.ph[pid]
    p = ph['paths'][0]
    Tc_prod, Tc_par = ph['Tc'], (p.get('parent') or {}).get('Tc')
    Tmag = max(x for x in (Tc_prod, Tc_par) if x is not None)
    def T_of():
        T = cond.get('T_process_K')
        try:
            return float(str(T).split(' ')[0].split('-')[0])
        except (TypeError, ValueError):
            return None
    if 'sc+fa' in r or (r.startswith('fa') and pid in ('L10-FePd', 'L10-CoPt')):
        T = T_of()
        if 'zero-field' in r or 'no field' in r or str(cond.get('field_T')) in ('0', '0.0'):
            return {'pred': pred(1 / 3, 1 / 3, 'no field: three equal variants (fraction 1/3)', cls='none')}
        if T is not None and T < Tmag:
            return {'pred': pred(1.0, 1 / 3, f'ordering at {T:.0f} K below the Curie point of the product or '
                                             f'ferromagnetic parent ({Tmag:.0f} K): field selects')}
        if T is not None:
            return {'pred': pred(1.0, 1 / 3, f'ordering at {T:.0f} K above both Curie points: only a '
                                             'paramagnetic bias (open)', cls='open')}
        return {'pred': pred(1.0, 1 / 3, 'ordering in a field, temperature not given', cls='info')}
    if 'sc, no field' in r or 'sc, no stress' in r or 'control' in r and 'sc' in r:
        return {'pred': pred(1 / 3, 1 / 3, 'unbiased single crystal: equal variants', cls='none')}
    if 'sc+sa' in r:
        return {'pred': pred(1.0, 1 / 3, 'ordering under uniaxial stress selects the variant')}
    if r.startswith('sa') or 'sa in films' in r or 'biaxial' in r:
        return {'pred': pred(0.8312, ISO, 'polycrystal ordered under stress: complete selection of the nearest '
                                          'variant (0.831 along the stress axis)')}
    if 'polytwinned' in r:
        return {'pred': pred(0.5774, 0.5774, 'polytwinned single crystal, three equal variants: best axis [111]',
                             cls='none')}
    if 'meteorite' in r:
        if 'interaction' in r:
            return {'pred': pred(1.0, 1 / 3, 'islands interact: selection by local magnetisation (not an applied '
                                             'bias); outside the single-island account', cls='info')}
        return {'pred': pred(1 / 3, 1 / 3, 'no external bias: on average equal variant fractions', cls='none')}
    if 'crystallisation of melt-spun amorphous' in r:
        return {'pred': pred(ISO, ISO, 'nucleation from an amorphous matrix without field: random', cls='none')}
    if 'rx+fa+sa' in r:
        return {'pred': pred(1.0, ISO, 'field and tension below Tc during ordering of a rolled strip: both act')}
    if r.startswith('fp'):
        if pid == 'L10-FeNi':
            return {'pred': pred(ISO + DELTA, ISO, 'NITE particles are polycrystalline: field alignment gives '
                                                   '~random', cls='none')}
        if al == 'perp':
            return {'pred': pred(ISO, ISO, 'aligned single-crystal particles measured across the axis', cls='none')}
        return {'pred': pred(1.0, ISO, 'single-crystal L1_0 particles aligned in a liquid')}
    if 'field re-anneal' in r or 'zero-field re-anneal' in r:
        return {'pred': pred(ISO, ISO, 'film already ordered before the field step: no formation event, the '
                                       'field cannot select', cls='none')}
    if 'field (0.6 t) during wet-chemical synthesis' in r or 'zero-field control of the xu' in r:
        return {'pred': pred(ISO, ISO, 'field during A1 synthesis; ordering later at 973 K without field: '
                                       'no selection', cls='none')}
    if r.startswith('ep'):
        return {'pred': pred(1.0, ISO, 'epitaxial film', cls='info')}
    return None


def enc_universe(ctx, pid, c, route, r, proc, cond, al, mk):
    if pid in ('L10-FePd', 'L10-CoPt'):
        return enc_l10(ctx, pid, c, route, r, proc, cond, al, mk)
    if pid == 'SrM-hexaferrite':
        return {'pred': pred(1.0, ISO, 'platelets with c normal: field or pressure orients them')}
    if pid == 'LTP-MnBi' or pid == 'Co-nanorods':
        return {'pred': pred(1.0, ISO, 'route not stated: only the absolute bound A <= 1', cls='info')}
    return None


def enc_fe16n2(ctx, pid, c, route, r, proc, cond, al, mk):
    if r.startswith('sa'):
        return {'pred': pred(1.0, ISO, 'stress during the α\'->α\'\' ordering anneal: acts only if the N '
                                       'sublattice re-orients at 423 K (open)', cls='open')}
    if 'ion implantation' in r:
        return {'pred': pred(1.0, ISO, 'implanted single-crystal foil + strained anneal: open as above', cls='open')}
    if r.startswith('fp') or 'electrospinning' in r:
        if al == 'perp':
            return {'pred': pred(ISO, ISO, 'aligned particles measured across the alignment axis', cls='none')}
        return {'pred': pred(1.0, ISO, 'single-variant α\'\' particles turned in a liquid')}
    if r.startswith('ep'):
        return {'pred': pred(1.0, ISO, 'epitaxial film', cls='info')}
    if 'niron' in r or "rx parent" in r:
        return {'pred': pred(1.0, ISO, 'process not disclosed / parent-only texture', cls='info')}
    return None


def enc_mnbi(ctx, pid, c, route, r, proc, cond, al, mk):
    if r.startswith('fp') or 'fp +' in r:
        if 'random crystallites' in str(cond.get('parent_texture')) or 'sb1.5' in proc:
            return {'pred': pred(ISO + DELTA, ISO, 'flakes contain randomly oriented crystallites: no net axis',
                                 cls='none')}
        if al == 'perp':
            return {'pred': pred(ISO, ISO, 'measured across the alignment axis', cls='none')}
        return {'pred': pred(1.0, ISO, 'single-crystal LTP particles aligned in a field')}
    if 'through htp' in r or 'alpha->beta->alpha' in r:
        return {'pred': pred(1.0, ISO, 'LTP re-forms from paramagnetic HTP on cooling in the field; whether c '
                                       'is inherited from HTP is unpublished', cls='open')}
    if r.startswith('fa') or 'reactive annealing' in r or 'field-assisted' in r:
        return {'pred': pred(1.0, ISO, 'LTP forms below its Curie point in the field: nuclei orient')}
    if 'above the ltp transition' in r:
        return {'pred': pred(1.0, ISO, 'paramagnetic HTP crystals turn in the melt by dchi; LTP inherits c '
                                       'only if the HTP->LTP OR keeps it (open)', cls='open')}
    if r.startswith('fs') or 'semi-solid' in r or 'high field gradient' in r:
        return {'pred': pred(1.0, ISO, 'ferromagnetic LTP crystals turn in liquid Bi')}
    if r.startswith('ds') or 'ds (zone' in r or 'ds +' in r:
        return {'pred': pred(1.0, ISO, 'eutectic rods grow with c along the growth axis (per MnBi volume)')}
    if 'film' in r or r.startswith('ep'):
        return {'pred': pred(1.0, ISO, 'film with c normal', cls='info')}
    return None


def enc_mnga_co(ctx, pid, c, route, r, proc, cond, al, mk):
    if 'cluster deposition' in r:
        if 'soft' in r:
            return {'pred': pred(1.0, ISO, 'aligned clusters + soft phase (exchange composite)', cls='info')}
        if al == 'perp':
            return {'pred': pred(ISO, ISO, 'aligned single-crystal clusters measured across the axis', cls='none')}
        return {'pred': pred(1.0, ISO, 'single-crystal clusters turned in flight')}
    if 'melt-spun' in r:
        return {'pred': pred(1.0, ISO, 'ribbon texture not reported', cls='info')}
    if pid in ('D022-Mn3Ga', 'D022-Mn3Ge') and ('fa' in r):
        T = cond.get('T_process_K')
        return {'pred': pred(1.0, ISO, f'D0_22 forms below its Curie point in the field (T {T} K)')}
    if pid == 'HfCo7' and r.startswith('fa'):
        return {'pred': pred(ISO, ISO, 'crystallisation at 873 K above Tc (751 K): no ferromagnetic bias',
                             cls='none')}
    if 'cold rolling + annealing without field' in r:
        return {'pred': pred(ISO, ISO, 'no field during formation', cls='none')}
    if r.startswith('ep'):
        return {'pred': pred(1.0, ISO, 'epitaxial film', cls='info')}
    return None


def enc_nano_other(ctx, pid, c, route, r, proc, cond, al, mk):
    medium = str(cond.get('alignment_medium', '')).lower() + ' ' + proc
    if pid == 'Co-nanorods':
        if 'field during growth' in r:
            return {'pred': pred(ISO, ISO, 'rods grown in a field but measured randomly embedded', cls='none')}
        if 'dry powder in die' in medium or 'dry agglomerated' in medium:
            return {'pred': pred(ISO + DELTA, ISO, 'dried rods agglomerate irreversibly: agglomerates have no '
                                                   'net axis', cls='none')}
        if al == 'perp':
            return {'pred': pred(ISO, ISO, 'measured across the alignment axis', cls='none')}
        if 'rescaled' in medium:
            return {'pred': pred(1.0, ISO, 'derived value (rescaled loop)', cls='info')}
        return {'pred': pred(1.0, ISO, 'single-crystal rods turned in a liquid (c along the rod)')}
    if pid == 'eps-Fe2O3' and r.startswith('ep'):
        return {'pred': pred(2 / 3, ISO, 'epitaxial film with three in-plane easy-axis variants: 0.667 ceiling')}
    if r.startswith('fp'):
        if al == 'perp':
            return {'pred': pred(ISO, ISO, 'measured across the alignment axis', cls='none')}
        return {'pred': pred(1.0, ISO, 'single-crystal particles turned in a field')}
    if r.startswith('ep'):
        return {'pred': pred(1.0, ISO, 'epitaxial film', cls='info')}
    return None


def enc_ferrite_ndfeb(ctx, pid, c, route, r, proc, cond, al, mk):
    if pid == 'Nd2Fe14B':
        if 'hddr' in r:
            return {'pred': pred(1.0, ISO, 'HDDR texture memory: ceiling depends on the memory phase (open)',
                                 cls='open')}
        if r.startswith('hd'):
            return {'pred': pred(1.0, ISO, 'die-upsetting: c along the press axis')}
        if 'hot pressing' in r:
            return {'pred': pred(1.0, ISO, 'hot pressing without upsetting: little strain', cls='info')}
        return None
    if 'calendering' in r or 'spss' in r or 'sps of' in r or 'shape alignment' in r or 'tgg' in r.lower() or 'shape-induced' in r:
        return {'pred': pred(1.0, ISO, 'platelets (c normal) oriented by pressure/flow or templates')}
    if r.startswith('fp') or 'rubber-bonded' in r or 'sintered' in r or 'conventional sintering' in r or 'not stated' in r:
        if al == 'perp':
            return {'pred': pred(ISO, ISO, 'measured across the alignment axis', cls='none')}
        return {'pred': pred(1.0, ISO, 'single-crystal particles aligned in a slurry/powder; sintering keeps it',
                             cls='acts' if r.startswith('fp') else 'info')}
    return None


def write():
    specs = build()
    hdr = ('# Conditional predictions for the sealed published cases (generated by align.cases from the\n'
           '# value-free condition files; committed before any outcome is opened).\n')
    open(os.path.join(ROOT, 'data', 'curated', 'cases.yaml'), 'w').write(
        hdr + yaml.safe_dump(specs, sort_keys=False, allow_unicode=True, width=110))
    return specs


if __name__ == '__main__':
    s = write()
    import collections
    print(len(s), 'cases')
    print(collections.Counter(x['pred']['cls'] for x in s))
    print(collections.Counter(x['status'] for x in s))
