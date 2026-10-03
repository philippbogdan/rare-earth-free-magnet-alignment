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


class Resolver:
    def __init__(self, reg):
        self.reg = reg
        self.used = set()
        self.estimates = []
        self.missing = []

    def __call__(self, obj, where=''):
        if isinstance(obj, dict):
            if 'ref' in obj and set(obj) <= {'ref', 'index', 'scale', 'note'}:
                rid = obj['ref']
                if rid not in self.reg.rec:
                    self.missing.append((where, rid))
                    return None
                self.used.add(rid)
                v = self.reg.value(rid)
                if 'index' in obj:
                    v = v[obj['index']]
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
    """All record ids referenced anywhere in a raw curated structure."""
    out = []
    if isinstance(obj, dict):
        if 'ref' in obj and isinstance(obj['ref'], str):
            out.append(obj['ref'])
        for v in obj.values():
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
        ph['_raw'] = p
        phases.append(ph)
    return phases, R
