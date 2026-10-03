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
import json
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
    'AlNiCo': 'AlNiCo-8-9', 'Fe-Cr-Co': 'Fe-Cr-Co',      # grade decides AlNiCo-5 vs 8/9 (enc_alnico_fecrco)
    'FeCoNiAlCu0.4Ti0.4 (alloy A)': 'HEA-FeCoNiAlCuTi', 'Fe2CoNiAlCu0.4Ti0.4 (alloy D)': 'HEA-FeCoNiAlCuTi',
    '(FeCoNiMnCu)Pt (L10 film)': 'L10-HEA-Pt', 'Fe3Se4 (monoclinic)': 'Fe3Se4',
    'Fe3-xMnxSe4 nanorods (composition for low-T value unspecified)': 'Fe3Se4',
    'Fe2+δTa-C14-meltspun': 'Fe2Ta-C14',
    'Mn4N': 'Mn4N', 'Mn-Sn-Co-N multiphase alloy': 'Mn4N',
    'Ba0.8La0.2Fe11.8Cu0.2O19': 'SrM-hexaferrite', 'Ba1-xCoxDyyFe12-yO19': 'SrM-hexaferrite',
    'Ba1-xCoxFe12O19': 'SrM-hexaferrite', 'RE-substituted SrFe12O19': 'SrM-hexaferrite',
    'Sr0.95Gd0.05Fe11.4Cu0.6O19': 'SrM-hexaferrite', 'Sr0.95Gd0.05Fe11.7Cu0.3O19': 'SrM-hexaferrite',
    'epsilon-Fe2O3': 'eps-Fe2O3', 'LTP-MnBi/FeSiB': 'LTP-MnBi', 'tau-MnAl-C/alpha-Fe': 'tau-MnAl-C',
}

PAIRS = [   # (case, reference case, expectation): same study, same alloy, one condition changed
    ('alnico8H_Br_FA_800_zhou2017', 'alnico8H_Br_nofield_800_zhou2017', 'higher'),
    ('alnico8H_Br_FA_840_zhou2017', 'alnico8H_Br_nofield_840_zhou2017', 'higher'),
    ('alnico8H_Br_FA_860_zhou2017', 'alnico8H_Br_nofield_860_zhou2017', 'higher'),
    ('alnico5_Br_equiaxed_field_makino1963', 'alnico5_Br_equiaxed_nofield_makino1963', 'higher'),
    ('alnico5_Br_columnar_field_par_makino1963', 'alnico5_Br_columnar_nofield_par_makino1963', 'higher'),
    ('alnico5_Br_columnar_field_perp_makino1963', 'alnico5_Br_columnar_nofield_perp_makino1963', 'higher'),
    ('alnico9_MrMs_Ydoped_field30kOe_zhang2018', 'alnico9_MrMs_Ydoped_field0kOe_zhang2018', 'higher'),
    ('alnico9_MrMs_Ydoped_textured_FA45kOe_zhang2018', 'alnico9_MrMs_Ydoped_field0kOe_zhang2018', 'higher'),
    ('fecrco_FA_stepaged_okadathesis', 'fecrco_iso_stepaged_okadathesis', 'higher'),
    ('hea_a_fa_mrms_na2021', 'hea_a_ta_mrms_na2021', 'open'),
    ('ndfeb_nguyen2013_field_mr_over_ms', 'ndfeb_nguyen2013_zero_mr_over_ms', 'equal'),
]

EXPOSED = {   # docs/PROTOCOL_LOG.md
    'prefix': ['fepd_scfa_', 'uni_out_rt_fepd_scfa', 'uni_out_fepd_l10_fa', 'mnal_hd_', 'mnal_sc_sa_kubo1976'],
    'ids': ['mnal_am_tiden2024_ref_mrms', 'mnal_ep_takata2018_sto', 'mnal_ep_takata2018_mgo_mn4n',
            'mnal_ep_takata2018_mgo', 'mnal_ep_takata2018_rocking', 'mnbi_other_tang2022_fieldsinter300',
            'mnbi_other_tang2022_fieldsinter250'],
}


KNOWN_TOPICS = ('alnico_fecrco', 'fese_hea', 'ndfeb_routes', 'mn4n', 'mnal', 'l10', 'universe', 'fe16n2',
                'mnbi', 'mnga_co', 'nano_other', 'ferrite_ndfeb', 'borides_feta')
PHASE_TOPIC = {   # rule set for each phase when a collection mixes phases (held-out batch 4)
    'tau-MnAl': 'mnal', 'tau-MnAl-C': 'mnal', 'LTP-MnBi': 'mnbi', 'L10-FeNi': 'l10', 'L10-FePt': 'l10',
    'L10-CoPt': 'l10', 'L10-FePd': 'l10', 'Fe16N2': 'fe16n2', 'D022-Mn3Ga': 'mnga_co', 'D022-Mn3Ge': 'mnga_co',
    'L10-MnGa': 'mnga_co', 'Zr2Co11': 'mnga_co', 'HfCo7': 'mnga_co', 'Co-nanorods': 'nano_other',
    'eps-Fe2O3': 'nano_other', 'Co3C': 'nano_other', 'bct-FeCo': 'nano_other', 'Fe5SiB2': 'nano_other',
    'Fe5PB2': 'nano_other', 'Fe2P-based': 'nano_other', 'SrM-hexaferrite': 'ferrite_ndfeb',
    'Nd2Fe14B': 'ndfeb_all', 'HEA-FeCoNiAlCuTi': 'fese_hea', 'L10-HEA-Pt': 'fese_hea', 'Fe3Se4': 'fese_hea',
    'Mn4N': 'mn4n', '(Fe,Co)2B': 'nano_other', 'Fe2Ta-C14': 'borides_feta',
}
TOPIC_FIRST = ('alnico_fecrco', 'fese_hea', 'ndfeb_routes', 'mn4n', 'ndfeb_all')


def exposed(cid):
    return cid in EXPOSED['ids'] or any(cid.startswith(p) for p in EXPOSED['prefix'])


def along_kind(s):
    s = (s or '').lower()
    if any(w in s for w in ('film normal', 'out-of-plane', 'out of plane', 'perpendicularly magnetized',
                            'perpendicular to cylinder plane', 'perpendicular to the cylinder plane')):
        return 'par'                          # the film/disc normal is the alignment axis
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
    if q in ('[redacted]', 'Br') and u.strip() in ('T', 'kG', 'kGs', 'mT', 'G', 'Gs'):
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
        self.js_same = {}

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


# ---------------------------------------------------------------- AlNiCo / Fe-Cr-Co rods
_D111 = np.array([[1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1]]) / np.sqrt(3)


def rods_A(tex, h, m, rule, spread=0.0):
    """A along m for α1 rods on the parent <100> (crystal x, y, z), per parent grain.

    tex   'random' | 'fibre' (<100> along z) | 'single' (parent axes = sample axes)
    h     field direction (sample frame); m measuring direction
    rule  'field'  rods follow the field in every grain (A = |h.m|)
          'near'   complete selection of the <100> nearest h (ties split equally)
          'none'   no bias: all three <100>
    Returns (A_indep, A_coupled): rods as independent units, or the selected rod sets of a
    grain as one coupled unit (one set: uniaxial; two: easy plane; three: cubic with the
    no-field easy directions <111>, Nesbitt 1952)."""
    h, m = np.asarray(h, float), np.asarray(m, float)
    h, m = h / np.linalg.norm(h), m / np.linalg.norm(m)
    if rule == 'field':
        a = abs(float(h @ m))
        return a, a
    g = {'random': lambda: E.texture_random(),
         'fibre': lambda: E.texture_fibre((0, 0, 1), n_az=720, spread_deg=spread),
         'single': lambda: E.texture_single(n=400 if spread else 1, spread_deg=spread)}[tex]()
    ax = np.transpose(g, (0, 2, 1))                 # (n, 3 rods, 3): columns of g
    if rule == 'near':
        sc = np.abs(ax @ h)
        sel = sc >= sc.max(1, keepdims=True) - 1e-6
    else:
        sel = np.ones(ax.shape[:2], bool)
    k = sel.sum(1)
    cm = np.abs(ax @ m)                             # (n, 3)
    indep = float(((cm * sel).sum(1) / k).mean())
    coup = np.where(k == 1, (cm * sel).sum(1), 0.0)
    two = k == 2
    if two.any():
        n_ax = (~sel[two]).argmax(1)                # the unselected axis is the plane normal
        nm = cm[two][np.arange(two.sum()), n_ax]
        coup[two] = np.sqrt(np.clip(1 - nm ** 2, 0, 1))
    three = k == 3
    if three.any():
        mc = np.einsum('nij,j->ni', np.transpose(ax[three], (0, 2, 1)), m)   # m in crystal frame
        coup[three] = np.abs(mc @ _D111.T).max(1)
    return indep, float(coup.mean())


def build():
    ctx = Ctx()
    out = []
    files = sorted(f for f in os.listdir(os.path.join(ROOT, 'data', 'sealed')) if f.endswith('.conditions.yaml'))
    man = json.load(open(os.path.join(ROOT, 'data', 'sealed', 'MANIFEST.json')))
    real = {f: [c['id'] for c in m['cases']] for b in man['batches'].values() for f, m in b['files'].items()}
    for fn in files:                         # same-sample Js records for Br -> A
        recs = yaml.safe_load(open(os.path.join(ROOT, 'data', 'sealed', fn))) or []
        ids = real.get(fn.replace('.conditions.yaml', '.yaml'), [])
        if len(ids) == len(recs):
            for c, i in zip(recs, ids):
                c['id'] = i
        js_ids = [x['id'] for x in recs if x.get('quantity') == 'Js_same_sample']
        for x in recs:
            if x.get('quantity') in ('Br', '[redacted]'):
                exact = x['id'].replace('_Br_', '_Js_same_')
                exact2 = x['id'].replace('_br_', '_js_same_sample_')
                paper = x['id'].rsplit('_', 1)[-1]
                same = [j for j in js_ids if j.endswith('_' + paper)]
                if exact in js_ids:
                    ctx.js_same[x['id']] = exact
                elif exact2 in js_ids:
                    ctx.js_same[x['id']] = exact2
                elif len(same) == 1:
                    ctx.js_same[x['id']] = same[0]
    man = json.load(open(os.path.join(ROOT, 'data', 'sealed', 'MANIFEST.json')))
    real = {f: [c['id'] for c in m['cases']] for b in man['batches'].values() for f, m in b['files'].items()}
    for fn in files:
        topic = fn.replace('_outcomes.conditions.yaml', '')
        recs = yaml.safe_load(open(os.path.join(ROOT, 'data', 'sealed', fn)))
        ids = real.get(fn.replace('.conditions.yaml', '.yaml'), [])
        if len(ids) == len(recs):              # value-free manifest ids (the redaction also hit ids)
            for c, i in zip(recs, ids):
                c['id'] = i
        for c in recs:
            spec = encode(ctx, topic, c)
            if spec:
                spec.update({'case': c['id'], 'topic': topic,
                             'status': 'exposed' if exposed(c['id']) else 'heldout'})
                out.append(spec)
    ids = {x['case'] for x in out}
    by = {x['case']: x for x in out}
    for a, b, exp in PAIRS:
        if a in ids and b in ids:
            by[a]['pair'] = {'with': b, 'expect': exp}
    return out


def encode(ctx, topic, c):
    """One case -> spec; rules by (phase, route) group with per-case conditions."""
    cid, route = c['id'], str(c.get('route'))
    pid = PHASE_MAP.get(c.get('phase'))
    if pid == 'tau-MnAl' and re.search(r'mn-?al-?c\b|mnal\(c\)|mn-al-c|mnal-c', str(c.get('process') or '').lower()):
        pid = 'tau-MnAl-C'                    # carbon-containing alloy labelled with the binary phase
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
    if topic not in KNOWN_TOPICS and pid is not None:      # mixed-phase collections: rules by phase
        topic = PHASE_TOPIC.get(pid, topic)
    t_first = 'alnico_fecrco' if (c.get('phase') in ('AlNiCo', 'Fe-Cr-Co')
                                  or pid in ('AlNiCo-5', 'AlNiCo-8-9', 'Fe-Cr-Co')) else topic
    if t_first in TOPIC_FIRST:                # topic rules that also cover unbiased samples
        res = globals()['enc_' + t_first](ctx, pid, c, route, r, proc, cond, al, mk)
        if res is not None:
            extra = res.pop('meas_extra', {})
            return {**base, **res, 'meas': {**base['meas'], **extra}}
    # ------------------------------------------------ isotropic samples (no alignment step)
    if r.startswith('isotropic') or 'no alignment' in r or 'without field' in r and 'cold rolling' not in r:
        if pid == 'Co-nanorods' and not any(w in proc for w in ('film', 'cast', 'spin-coated')):
            return {**base, 'pred': _rod_arrangement('rod powder or compact without an alignment step')}
        if any(w in proc for w in ('film', 'cast', 'spin-coated')) and pid == 'Co-nanorods' and al != 'any':
            return {**base, 'pred': pred(PLANAR, PLANAR, 'rods lie in the film plane at random azimuth: '
                                                         'planar random, in-plane value 2/pi', cls='none')}
        exch = any(w in proc for w in ('milled', 'melt-spun', 'ribbon', 'nanocrystal', 'quenched', 'mechanically alloyed',
                                       'cluster', 'mqp', 'bonded nd-fe-b', 'fesib', 'alpha-fe'))
        if '/' in str(c.get('phase') or ''):          # composite with a second (soft) phase
            exch = True
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
            if float(t1) == 0:
                return {'pred': pred(0.6133, 0.6133, 'compression along c_ε: all six variants equivalent', cls='none')}
            return {'pred': pred(1.0, 0.6133, f'v2: compression at theta1={t1}, theta2={t2} deg; the tetragonal strain '
                                              f'alone leaves {n} variant(s) tied (A {A:.3f}), the transformation shear '
                                              'can select one (SC+SA cell ceiling)')}
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
    if 'ds + hot deformation' in r:
        return {'pred': _hyp({'deformation after the transformation (massive)': 0.6133,
                              'stress during a displacive transformation': 0.899},
                             'v2: DS + hot deformation, temperature and sequence not stated', cls='open')}
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
        return {'pred': _hyp({'independent variants': 0.5774, 'coupled fine polytwins (multiaxial)': 1.0},
                             'v2: polytwinned crystal with three equal variants; fine exchange-coupled twins act as one '
                             'multiaxial unit (as the spinodal networks)', cls='none', base=0.5774)}
    if 'meteorite' in r:
        if 'interaction' in r:
            return {'pred': pred(1.0, 1 / 3, 'islands interact: selection by local magnetisation (not an applied '
                                             'bias); outside the single-island account', cls='info')}
        return {'pred': pred(1.0, 1 / 3, 'v2: ordering below the Curie point of the ferromagnetic parent: the '
                                         "parent's own domain magnetisation biases the variant (rule 2 of METHOD "
                                         'section 8), so one variant can dominate')}
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
        if 'melt-spun' in proc and 'mill' in proc:
            return {'pred': _hyp({'particles milled below the ribbon grain size (single crystals)': 1.0,
                                  'particles holding several randomly oriented ribbon grains': ISO + DELTA},
                                 'powder milled from melt-spun ribbon: particle size relative to the ribbon grain '
                                 'size unpublished', cls='open')}
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


def _rod_arrangement(why):
    """v2: rods in an assembly made without a deliberate alignment step: their arrangement is unverified."""
    return _hyp({'3-D random': ISO, 'lying flat at random azimuth (in-plane)': PLANAR,
                 'parallel bundles (dipolar self-assembly)': 1.0},
                why + ': rod arrangement in the measured sample not verified', cls='open')


def enc_nano_other(ctx, pid, c, route, r, proc, cond, al, mk):
    medium = str(cond.get('alignment_medium', '')).lower() + ' ' + proc
    if pid == 'Co-nanorods':
        if 'field during growth' in r:
            return {'pred': _rod_arrangement('rods grown in a field, then embedded without deliberate alignment')}
        if 'dry powder in die' in medium or 'dry agglomerated' in medium:
            return {'pred': pred(1.0, ISO, 'v2: dried rods agglomerate into bundles of parallel rods, which have a net '
                                           'axis and turn in the die field (incompletely)')}
        if al == 'perp':
            return {'pred': pred(ISO, ISO, 'measured across the alignment axis', cls='none')}
        if 'rescaled' in medium:
            return {'pred': pred(1.0, ISO, 'derived value (rescaled loop)', cls='info')}
        return {'pred': pred(1.0, ISO, 'single-crystal rods turned in a liquid (c along the rod)')}
    if pid == 'eps-Fe2O3' and r.startswith('ep') and mk == 'vfrac':
        return {'pred': pred(1.0, ISO, 'fraction of one in-plane variant family of a multi-variant film (not A)',
                             cls='info')}
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


# ---------------------------------------------------------------- AlNiCo / Fe-Cr-Co
_Z, _X = (0, 0, 1), (1, 0, 0)
_JS = {   # sourced saturation intervals (T) by grade, data/raw/alnico_fecrco_inputs.yaml
    'alnico1-4': ([0.87, 1.46], 'grade Js unpublished: whole sourced alnico range (estimate)'),
    'alnico5': ([1.43, 1.46], 'alnico5_Js_iwama1976, alnico_Js_sergeyev_via_takeuchi1976_5 (Ti-free grades)'),
    'alnico6': ([1.30, 1.43], 'between alnico_Js_takeuchi1976_ti1p3 (1.36) and the Ti-free value (estimate)'),
    'alnico8n': ([1.15, 1.23], 'alnico_Js_sergeyev_via_takeuchi1976_8 (5 wt% Ti, MMPA alnico 8 and 9) -/+3%'),
    'alnico8H': ([1.03, 1.10], 'alnico_Js_takeuchi1976_ti7p5, alnico_Js_zou2016_DS_7p5Ti (7-8.7 wt% Ti)'),
    'alnico8': ([1.03, 1.23], 'Ti content not stated: from the 7.5 wt% Ti to the 5 wt% Ti value'),
}


def _alnico_grade(cid, proc):
    k = cid.lower()
    if re.search(r'alnico[1-4]_|alnico_br_alnico[1-4]_|ln10|lng12', k):
        return 'AlNiCo-5', 'alnico1-4'
    if 'alnico6' in k:
        return 'AlNiCo-5', 'alnico6'
    if re.search(r'alnico5|makino1963|no1_takeuchi|alcomax|5dg|_57_|arkomax800|lng52|lng60', k):
        return 'AlNiCo-5', 'alnico5'
    if re.search(r'alnico9|_9_|lngt60|lngt72|lngt80', k) or ('_mmpa' in k and 'alnico8_' in k):
        return 'AlNiCo-8-9', 'alnico8n'
    if re.search(r'alnico8h|8hc|8he|8h_|kassen2018|white20|ames2023|rinko2023|zou2016|no4_takeuchi', k):
        return 'AlNiCo-8-9', 'alnico8H'
    return 'AlNiCo-8-9', 'alnico8'


def _fecrco_js(cid, proc):
    if 'fe31cr23co' in cid or 'okada' in cid or 'fe-31cr-23co' in proc:
        return [1.2, 1.4], 'fecrco_Js_fe31cr23co_okadathesis'
    if re.search(r'\d+\s*-\s*\d+\s*cr', proc) and 'arnokrome' not in cid:
        return [1.2, 1.7], 'composition range only: sourced Fe-Cr-Co range (estimate)'
    if 'arnokrome' in cid:
        proc = 'fe-28cr-8co'
    m = re.search(r'(\d+(?:\.\d+)?)\s*(?:wt%\s*)?cr\b', proc) or re.search(r'fe-?\(?(\d+)', proc)
    cr = None
    if m:
        try:
            cr = float(m.group(1))
        except ValueError:
            cr = None
    if cr is None or not 10 <= cr <= 40:
        return [1.2, 1.7], 'composition not resolved: sourced Fe-Cr-Co range (estimate)'
    js = float(np.interp(cr, [18.4, 26.8, 35.2], [1.68, 1.37, 1.21]))
    return [round(js - 0.12, 2), round(js + 0.12, 2)], \
        f'interpolated in Cr ({cr:g} wt%) between fecrco_Js_he2024_No1/No4/No7, +-0.12 T for Co (estimate)'


def _hyp(vals, why, cls=None, base=None):
    """Prediction carrying alternative readings: vals = {name: A}; A_hi = max."""
    hi = max(vals.values())
    b = base if base is not None else ISO
    if cls is None:
        cl = {('acts' if v > b + DELTA else 'none') for v in vals.values()}
        cls = cl.pop() if len(cl) == 1 else 'open'
    out = pred(hi, b, why, cls=cls)
    out['hyp'] = {k: round(float(v), 4) for k, v in vals.items()}
    return out


def _unbiased(tex, m, why, spread=0.0):
    ind, cou = rods_A(tex, _Z, m, 'none', spread)
    return _hyp({'independent rods': ind, 'coupled network (<111> easy)': cou},
                'no bias: α1 along all three <100>; ' + why, cls='none', base=ind)


def enc_alnico_fecrco(ctx, pid, c, route, r, proc, cond, al, mk):
    cid = c['id']
    alongs = str(c.get('measured_along') or '').lower()
    ptex = str(cond.get('parent_texture') or '').lower()
    js_case = ctx.js_same.get(cid)
    meas = {}
    if c.get('phase') == 'AlNiCo':
        pid, grade = _alnico_grade(cid, proc)
        js, js_src = _JS[grade]
    else:
        pid, grade = 'Fe-Cr-Co', 'fecrco'
        js, js_src = _fecrco_js(cid, proc)
    if mk == 'Br':
        meas = {'js_case': js_case} if js_case else {'js': js, 'js_src': js_src}
        if 'bonded' in proc:
            meas['phi'] = [0.5, 1.0]
        elif 'sinter' in proc or 'sint_' in cid or '_fl' in cid:
            meas['phi'] = [0.95, 1.0]
    out = lambda p: {'phase': pid, 'pred': p, 'meas_extra': meas}
    no_field = str(cond.get('field_T')) in ('0', '0.0') or 'without field' in proc or 'isotropic' in r \
        or re.search(r'no field(?! template)', r + ' ' + proc) is not None
    colum = ('columnar' in ptex and not re.search(r'not (described as )?columnar|non-?columnar', ptex)
             and 'columnar-appearing' not in ptex and 'not measured' not in ptex) or r.startswith('ds')
    fe = pid == 'Fe-Cr-Co'
    mo = fe and bool(re.search(r'\d\s*mo\b|-mo\b|mo-|\dmo', proc + ' ' + cid))
    mo_wt = None
    if mo:
        mm = re.search(r'(\d+(?:\.\d+)?)\s*mo', proc + ' ' + cid.lower())
        mo_wt = float(mm.group(1)) if mm else None
    mo_low = mo and mo_wt is not None and mo_wt < 3   # v2: path 1 is sourced at 3 wt% Mo only
    # ---------------------------------------------------- not tests of alignment
    if mk == 'class' or c.get('quantity') == 'texture_metric':
        return out(pred(1.0, ISO, 'parent texture metric from the same work that supplied the account input '
                                  '(not independent)', cls='info'))
    if 'no coercivity' in proc or 'as-quenched' in cid or 'asquenched' in cid:
        return out(pred(1.0, ISO, 'no α1 rods formed (no coercivity): not a Stoner-Wohlfarth remanence', cls='info'))
    if 'whether a field was applied is not stated' in r or 'post-treatment not give' in proc:
        return out(pred(1.0, ISO, 'process not specific enough to test', cls='info'))
    if 'cold roll' in proc and no_field and fe:
        return out(pred(1.0, ISO, 'cold-rolled strip aged without field: deformation state and texture not '
                                  'reported', cls='info'))
    # ---------------------------------------------------- deformation ageing (Fe-Cr-Co)
    if re.search(r'deformation[- ]ag', r):
        cell = ctx.cell('Fe-Cr-Co', 'HD')
        return out(pred(cell.A_max, ISO, 'deformation ageing elongates α1 along the deformation axis (HD cell)'))
    # ---------------------------------------------------- stress anneal without field (AlNiCo)
    if r.startswith('sa'):
        cell = ctx.cell(pid, 'SA')
        return out(pred(cell.A_max, ISO, 'stress during the coherent spinodal on a random parent: complete '
                                         'selection of the <100> nearest the tensile axis (SA cell)'))
    # ---------------------------------------------------- unbiased (no field, no stress)
    if no_field:
        if colum:
            m = _X if al == 'perp' else _Z
            return out(_unbiased('fibre', m, 'columnar <100> fibre', spread=2.5))
        if 'as-built columnar' in ptex:
            return out(_unbiased('random', _Z, 'AM parent: random to <001> fibre (account input); random '
                                               'gives the larger coupled value'))
        return out(_unbiased('random', _Z, 'random or weakly textured parent'))
    # ---------------------------------------------------- field-treated
    if fe:
        mode = 'Mo-bearing: rods on the nearest <100> (path 1)' if mo else 'Mo-free: rods follow the field (path 0)'
        if 'sc+fa' in r or 'single crystal' in ptex:
            fd = re.search(r'field along <(\d)(\d)(\d)>', ptex)
            if fd:
                d = np.array([float(x) for x in fd.groups()])
                if not mo:
                    b, _ = rods_A('single', d, d, 'none')
                    return out(pred(1.0, b, f'{mode}: A = 1 along the field <{"".join(fd.groups())}>'))
                ind, cou = rods_A('single', d, d, 'near')
                b, _ = rods_A('single', d, d, 'none')
                return out(_hyp({'independent rods': ind, 'coupled sets (tie)': cou},
                                f'{mode}; field along <{"".join(fd.groups())}>: nearest-<100> selection '
                                f'(ties split)', base=b))
            return out(pred(1.0, ISO, f'{mode}; single crystal with field along <100>'))
        if 'rx+fa' in r or 'ds+fa' in r or 'am+fa' in r:
            return out(pred(1.0, ISO, f'{mode}; parent <100> texture along the field (or rods follow the '
                                      'field): ceiling 1'))
        # equiaxed / sintered / unstated polycrystal
        if mo_low and not ('rx+fa' in r or 'ds+fa' in r or 'am+fa' in r or 'single crystal' in ptex):
            return out(_hyp({'rods follow the field (path 0)': 1.0, 'rods on the nearest <100> (path 1)': 0.8312},
                            f'v2: {mo_wt:g} wt% Mo, below the 3 wt% of the path-1 sources; random parent'))
        if mo and ptex.startswith('not stated'):
            return out(pred(1.0, ISO, f'{mode}; parent texture not stated: ceiling of a <100> fibre along the field'))
        if mo:
            return out(pred(0.8312, ISO, f'{mode}; random parent: nearest <100> (0.831)'))
        if 'mo' not in proc and re.search(r'15-35 cr', proc):
            return out(_hyp({'Mo-free (path 0)': 1.0, 'Mo-bearing (path 1)': 0.8312},
                            'commercial grade, Mo content not stated'))
        return out(pred(1.0, ISO, f'{mode}; any grain orientation'))
    # AlNiCo, field-treated
    five = pid == 'AlNiCo-5'
    if al == 'perp' and not colum:
        near_perp = rods_A('random', _Z, _X, 'near')[0]
        _, unb = rods_A('random', _Z, _X, 'none')
        vals = {'rods on the nearest <100>': max(near_perp, unb)}
        if five:
            vals = {'rods follow the field': unb}
        return out(_hyp(vals, f'measured across the anneal field: complete selection gives {near_perp:.3f} '
                              f'(nearest <100>) or 0 (rods along the field); partial selection up to the '
                              f'unbiased {unb:.3f}', cls='none', base=ISO))
    if 'grain structure not described' in r:
        return out(pred(1.0, ISO, 'grain structure not stated', cls='info'))
    if 'agg under uniaxial load' in ptex:
        return out(pred(1.0, ISO, 'parent texture from stress-biased grain growth not quantified', cls='info'))
    if 'sc+fa' in r:
        return out(pred(1.0, 1 / 3, 'single crystal, field along [100]: one rod set along the field'))
    if colum and ('angle' in cid):
        th = float(re.search(r'angle(\d+)', cid).group(1))
        t = np.radians(th)
        if 'perp' in cid:
            h, m = (np.cos(t), 0, np.sin(t)), _X
        else:
            h, m = (np.sin(t), 0, np.cos(t)), _Z
        vals = [rods_A('fibre', h, m, 'near', s) for s in (0.0, 2.5)]
        near = max(v[0] for v in vals)
        b, unb = rods_A('fibre', _Z, m, 'none', 2.5)
        # off the field direction partial selection can exceed complete selection: the ceiling is
        # the larger of the complete-selection and the unbiased (coupled) values
        full = [float(np.cos(t)), near]                 # complete selection decides the class
        cl = {('acts' if v > b + DELTA else 'none') for v in full}
        cl = {('acts' if float(np.cos(t)) > b + DELTA else 'none')}
        return out(_hyp({'rods follow the field': max(float(np.cos(t)), unb)}, cls=cl.pop(),
                        why=
                        f'columnar <100> fibre, field {th:g} deg from the measuring direction; '
                        f'complete selection gives {np.cos(t):.3f} (follow) / {near:.3f} (nearest <100>), '
                        f'unbiased up to {unb:.3f}', base=b))
    if 'columnar-appearing' in ptex or ('as above' in ptex and 'l37' in cid.lower()):
        return out(pred(1.0, ISO, 'printed parent texture not measured: ceiling of a perfect <100> fibre',
                        cls='acts'))
    if colum:
        if 'field_perp' in cid or 'perpendicular to the columnar' in alongs or 'perpendicul' in str(cond.get('magnet_axis')):
            near = rods_A('fibre', _X, _X, 'near', 2.5)[0]
            b, _ = rods_A('fibre', _Z, _X, 'none', 2.5)
            vals = {'rods follow the field': 1.0}
            return out(_hyp(vals if five else {'rods on the nearest <100>': near},
                            'columnar <100> fibre, field and measurement across the columns', base=b))
        return out(pred(1.0, 1 / 3, 'columnar <100> fibre with the field along it: one rod set along the field'))
    if '(200)' in ptex or '<200>' in ptex:
        return out(pred(1.0, ISO, '(200)-textured button, field along <200>: between 0.831 (random) and 1 '
                                  '(perfect fibre); texture degree not quantified'))
    if 'near single crystal with [111]' in ptex:
        d = np.ones(3)
        ind, cou = rods_A('single', d, d, 'near')
        return out(_hyp({'independent rods': ind, 'coupled sets (tie)': cou},
                        'near single crystal with [111] along the field: the three <100> tie', base=ind))
    if 'field-templated' in ptex or 'agg under uniaxial load' in ptex:
        return out(pred(1.0, ISO, 'parent texture from field templating / stress-biased growth not quantified',
                        cls='info'))
    if 'grain structure not described' in r:
        return out(pred(1.0, ISO, 'grain structure not stated', cls='info'))
    weak = 'weak texture' in ptex or 'fairly isotropic' in ptex
    if five:
        return out(pred(1.0, ISO, 'v2: Ti-free alnico, rods follow the field in any grain (settled by the held-back '
                                  'data: equiaxed Makino 1963 0.895 > 0.831)'))
    if weak:
        return out(pred(0.8312 + 0.1 * (1 - 0.8312), ISO, 'nearly random parent (very weak texture, <=10% '
                                                         '<100> fibre assumed): nearest <100>'))
    return out(pred(0.8312, ISO, 'random parent: rods on the <100> nearest the field (0.831)'))


# ---------------------------------------------------------------- HEA, Fe3Se4
def enc_fese_hea(ctx, pid, c, route, r, proc, cond, al, mk):
    cid = c['id']
    if pid == 'HEA-FeCoNiAlCuTi':
        js = [0.49, 0.54] if 'alloy a' in str(c.get('phase')).lower() else [0.86, 0.95]
        meas = {'js': js, 'js_src': 'hea_a_js_na2021 / hea_d_js_na2021 (+-5%)'} if mk == 'Br' else {}
        if r.startswith('isotropic'):
            res = _unbiased('random', _Z, 'random as-cast parent (alnico-like spinodal)')
        else:
            res = _hyp({'precipitates follow the field': 1.0, 'precipitates on the nearest <100>': 0.8312,
                        'no elongated precipitates': 0.8656},
                       'field anneal at 923 K below Tc 1093 K: shape anisotropy and its axis unpublished',
                       cls='open')
        return {'pred': res, 'meas_extra': meas}
    if pid == 'L10-HEA-Pt':
        return {'pred': pred(1.0, ISO, 'rocking width of the (001)-oriented grains does not give the oriented '
                                       'fraction', cls='info')}
    if pid == 'Fe3Se4' and 'alignment unreported' in r:
        return {'pred': pred(1.0, ISO, 'alignment not reported', cls='info')}
    return None


# ---------------------------------------------------------------- Nd2Fe14B, other routes
def enc_ndfeb_routes(ctx, pid, c, route, r, proc, cond, al, mk):
    cid = c['id']
    if mk == 'skip':
        return None
    if 'injection cast' in proc:
        meas = {'phi': [0.5, 1.0]} if mk == 'Br' else {}
        if str(cond.get('field_T')) in ('0', '0.0'):
            return {'pred': pred(ISO, ISO, 'injection cast without field: random grains', cls='none'),
                    'meas_extra': meas}
        return {'pred': pred(1.0, ISO, 'crystals solidifying in 0.37 T far above Tc: paramagnetic torque on '
                                       '>~1 um crystals exceeds kT, but rotation (~ms) competes with the '
                                       'solidification time (unpublished)', cls='open'), 'meas_extra': meas}
    if 'semisolid' in proc or r.startswith('fs'):
        meas = {'phi': [0.25, 0.35], 'phi_src': '30 wt% Nd2Fe14B in Nd70Cu30 (volume fraction estimate)'} \
            if mk == 'Br' else {}
        return {'pred': pred(1.0, ISO, 'Nd2Fe14B crystals free in a liquid at 973 K (paramagnetic) turn in 10 T '
                                       '(FS cell)'), 'meas_extra': meas}
    if 'melt-spun' in proc and 'powdered' in proc:
        if mk == 'MrMs':
            return {'pred': pred(1.0, ISO, 'powdered nanocomposite ribbon measured at random orientation: '
                                           'exchange-enhanced, not an alignment value; the 0.32 T field acts on '
                                           'nanocrystals far above Tc (bias << kT): no field effect', cls='info')}
        return {'pred': pred(1.0, ISO, 'powdered nanocomposite ribbon', cls='info')}
    if 'melt-spun' in proc:
        return {'pred': pred(1.0, ISO, 'melt-spun ribbon: growth texture unpublished, measuring direction relative '
                                       'to the field not stated', cls='info')}
    if r.startswith('am'):
        meas = {'phi': [0.5, 1.0]} if mk == 'Br' else {}
        if 'bonded' in proc or 'isotropic' in str(cond.get('parent_texture')).lower():
            return {'pred': pred(ISO, ISO, 'printed bonded magnet of isotropic rapidly solidified powder: no alignment '
                                           'step; exchange enhancement possible (not alignment)', cls='info'),
                    'meas_extra': meas}
        return {'pred': pred(1.0, ISO, 'printed Nd-Fe-B: crystal orientation not measured (AM cell open)',
                             cls='open'), 'meas_extra': meas}
    return None


def enc_ndfeb_all(ctx, pid, c, route, r, proc, cond, al, mk):
    """Nd2Fe14B in a mixed collection: the sintered/HD/HDDR rules, then the other-route rules."""
    res = enc_ndfeb_routes(ctx, pid, c, route, r, proc, cond, al, mk)
    if res is None and not r.startswith('isotropic'):
        res = enc_ferrite_ndfeb(ctx, pid, c, route, r, proc, cond, al, mk)
    return res


def enc_mn4n(ctx, pid, c, route, r, proc, cond, al, mk):
    iso = ctx.cell('Mn4N', 'FA').A_max            # random cubic multiaxial value
    if 'multiphase' in str(c.get('phase')).lower() or 'mn-sn-co-n' in proc.lower() or 'mn82.5' in proc.lower():
        return {'pred': pred(1.0, iso, 'multiphase Mn-Sn-Co-N alloy: fraction and Js of the Mn4N-like phase '
                                       'unpublished', cls='info')}
    if r.startswith('ep'):
        if al == 'perp':
            return {'pred': pred(iso, iso, 'film measured in plane', cls='none')}
        return {'pred': pred(1.0, iso, 'epitaxial film on MgO(001), tetragonally distorted: uniaxial easy axis '
                                       'along the film normal (EP cell)')}
    return None


def build_known():
    """Specs for the known (fidelity) cases of PROBLEM.md, data/raw/known_cases_*.yaml, from the same
    value-free fields and the same rules as the sealed cases."""
    ctx = Ctx()
    keep = ('id', 'phase', 'route', 'process', 'conditions', 'quantity', 'unit', 'measured_along', 'measurement')
    out = []
    for fn in sorted(os.listdir(os.path.join(ROOT, 'data', 'raw'))):
        if not fn.startswith('known_cases_'):
            continue
        topic = fn[len('known_cases_'):-5]
        for r in yaml.safe_load(open(os.path.join(ROOT, 'data', 'raw', fn))) or []:
            c = {k: r.get(k) for k in keep}
            if not isinstance(c.get('conditions'), dict):
                c['conditions'] = {}
            spec = encode(ctx, topic, c)
            if spec:
                spec.update({'case': c['id'], 'topic': topic, 'status': 'known'})
                out.append(spec)
    return out


def write(version='v2'):
    """Write the current rules' predictions. The v1 files (cases.yaml, cases_known.yaml) are the
    frozen pre-registration (docs/PROTOCOL_LOG.md) and are never rewritten."""
    assert version != 'v1', 'v1 is frozen'
    known = build_known()
    open(os.path.join(ROOT, 'data', 'curated', f'cases_known_{version}.yaml'), 'w').write(
        f'# Known-case predictions, rules {version} (post hoc; see docs/PROTOCOL_LOG.md).\n'
        + yaml.safe_dump(known, sort_keys=False, allow_unicode=True, width=110))
    specs = build()
    open(os.path.join(ROOT, 'data', 'curated', f'cases_{version}.yaml'), 'w').write(
        f'# Case predictions, rules {version}: revised AFTER the outcomes were opened (post hoc, not a\n'
        '# held-out test); every change is listed in docs/PROTOCOL_LOG.md.\n'
        + yaml.safe_dump(specs, sort_keys=False, allow_unicode=True, width=110))
    return specs


if __name__ == '__main__':
    s = write()
    import collections
    print(len(s), 'cases')
    print(collections.Counter(x['pred']['cls'] for x in s))
    print(collections.Counter(x['status'] for x in s))
