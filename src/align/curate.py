"""Load data/curated/phases.yaml, resolving every {ref: <record id>} to the sourced value.

A field may be:
  {ref: id}                -> value of the record (number or text)
  {ref: id, index: k}      -> k-th element of a list-valued record
  {est: value, why: text}  -> an estimate made here (labelled as such in all outputs)
  plain values             -> crystallographic conventions (Miller indices, Laue class),
                              never physical magnitudes.
"""
import os

import yaml

from .data import ROOT, Registry


def _numeric(v):
    """Numeric text -> float; 'lo-hi' text -> [lo, hi]; other values unchanged."""
    import re
    if isinstance(v, str):
        t = v.strip().replace('−', '-').replace('–', '-')
        m = re.match(r'^([0-9]+(?:\.[0-9]+)?)\s*-\s*([0-9]+(?:\.[0-9]+)?)$', t)
        if m:
            return [float(m.group(1)), float(m.group(2))]
        m = re.match(r'^[<>~≈]?\s*([0-9]+(?:\.[0-9]+)?)$', t)
        if m:
            return float(m.group(1))
    if isinstance(v, list):
        return [_numeric(x) for x in v]
    return v


def _key(v, k, where):
    """Component k of a record value given as a mapping {'a': ..} or a text 'a = 5.72, c = 6.29'."""
    import re
    if isinstance(v, dict):
        x = v[k]
    else:
        m = re.search(r'(?<![A-Za-z])' + re.escape(k) + r'\s*(?:\(primitive\))?\s*=\s*([-0-9.]+)', str(v))
        if not m:
            raise KeyError(f'{where}: no {k} in {v!r}')
        x = m.group(1)
    if isinstance(x, str):
        x = x.split('-')[0] if re.match(r'^[0-9.]+-[0-9.]+$', x) else x
    return float(x)


class Resolver:
    def __init__(self, reg):
        self.reg = reg
        self.used = set()
        self.estimates = []
        self.missing = []

    def __call__(self, obj, where=''):
        if isinstance(obj, dict):
            if 'ref' in obj and set(obj) <= {'ref', 'index', 'key', 'scale', 'note', 'pick'}:
                rid = obj['ref']
                if rid not in self.reg.rec:
                    self.missing.append((where, rid))
                    return None
                self.used.add(rid)
                v = _numeric(self.reg.value(rid))
                if 'key' in obj:
                    v = _key(v, obj['key'], where)
                if 'index' in obj:
                    v = v[obj['index']]
                if 'pick' in obj:                       # 'lo' | 'hi' of a [lo, hi] range
                    v = v[0] if obj['pick'] == 'lo' else v[-1]
                if 'scale' in obj:
                    v = float(v) * obj['scale']
                return v
            if 'est' in obj and set(obj) <= {'est', 'why'}:
                self.estimates.append((where, obj['est'], obj.get('why', '')))
                return obj['est']
            return {k: self(v, f'{where}.{k}') for k, v in obj.items()}
        if isinstance(obj, list):
            return [self(v, f'{where}[{i}]') for i, v in enumerate(obj)]
        return obj


def refs_in(obj):
    """All record ids referenced anywhere in a raw curated structure ({ref: id}, refs: [ids],
    *_ref: id, *_refs: [ids])."""
    out = []
    if isinstance(obj, dict):
        if 'ref' in obj and isinstance(obj['ref'], str):
            out.append(obj['ref'])
        for k, v in obj.items():
            if (k == 'refs' or k.endswith('_refs')) and isinstance(v, list):
                out += [x for x in v if isinstance(x, str)]
            elif k.endswith('_ref') and isinstance(v, str):
                out.append(v)
            else:
                out += refs_in(v)
    elif isinstance(obj, list):
        for v in obj:
            out += refs_in(v)
    return out


def load_phases(reg=None, path=None):
    reg = reg or Registry()
    path = path or os.path.join(ROOT, 'data', 'curated', 'phases.yaml')
    raw = yaml.safe_load(open(path)) or []
    R = Resolver(reg)
    phases = []
    for p in raw:
        ph = R(p, p['id'])
        ph['_refs'] = sorted(set(refs_in(p)))
        for rid in ph['_refs']:
            if rid not in reg.rec:
                R.missing.append((p['id'], rid))
            else:
                R.used.add(rid)
        ph['_raw'] = p
        derive(ph, R)
        phases.append(ph)
    return phases, R


def cell_volume(st):
    a, c, b = st['a'], st.get('c', st['a']), st.get('b', st['a'])
    if st['system'] == 'hexagonal':
        return (3 ** 0.5 / 2) * a * a * c * 1e-30
    if st['system'] == 'cubic':
        return a ** 3 * 1e-30
    return a * b * c * 1e-30


def derive(ph, R):
    """Quantities derived here (estimates, labelled): magnetic-atom density and Js at 0 K."""
    from .physics import ms_ratio
    cell = ph.get('cell') or {}
    if cell.get('atoms') and ph['structure'].get('a') not in (None, 1.0):
        ph['n_mag'] = cell['atoms'] * float(cell.get('magnetic_fraction', 1.0)) / cell_volume(ph['structure'])
        R.estimates.append((ph['id'] + '.n_mag', f"{ph['n_mag']:.3g} m^-3",
                            'cell content x magnetic fraction / cell volume (bounds the paramagnetic bias only)'))
    if ph.get('Js') is not None and ph.get('Tc') is not None:
        ph['Js0'] = float(ph['Js']) / float(ms_ratio(300.0, float(ph['Tc'])))
        R.estimates.append((ph['id'] + '.Js0', f"{ph['Js0']:.3f} T",
                            "Js(300 K) / Kuz'min Ms(T) shape (s = 0.5, p = 5/2); bounds the paramagnetic bias only"))
