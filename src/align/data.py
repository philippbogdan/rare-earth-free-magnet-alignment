"""Record registry: every number used by the account is a sourced record (id -> value,
unit, source, quote, evidence) from data/raw/*.yaml, data/curated/records.yaml or opened
outcome files; curated phase descriptions refer to records by id."""
import glob
import json
import os

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _load(fn):
    recs = yaml.safe_load(open(fn)) or []
    for r in recs:
        r['_file'] = os.path.relpath(fn, ROOT)
    return recs


class Registry:
    def __init__(self, include_outcomes=True):
        self.rec = {}
        files = sorted(glob.glob(os.path.join(ROOT, 'data', 'raw', '*.yaml')))
        files += sorted(glob.glob(os.path.join(ROOT, 'data', 'curated', 'records*.yaml')))
        if include_outcomes:
            files += sorted(glob.glob(os.path.join(ROOT, 'data', 'outcomes', '*.yaml')))
        self.dups = []
        for fn in files:
            for r in _load(fn):
                rid = r.get('id')
                if not rid:
                    continue
                if rid in self.rec:
                    self.dups.append(rid)
                self.rec[rid] = r
        vf = os.path.join(ROOT, 'data', 'verification.json')
        self.ver = json.load(open(vf)) if os.path.exists(vf) else {}
        mf = os.path.join(ROOT, 'data', 'verification_manual.yaml')
        for m in (yaml.safe_load(open(mf)) or []) if os.path.exists(mf) else []:
            if m.get('verdict') == 'verified':          # checked by eye where the automatic match failed
                v = dict(self.ver.get(m['id']) or {})
                v.update({'quote_ok': True, 'level': 'manual', 'manual_how': m.get('how')})
                self.ver[m['id']] = v

    def __getitem__(self, rid):
        return self.rec[rid]

    def value(self, rid):
        return self.rec[rid]['value']

    def cite(self, rid):
        r = self.rec[rid]
        s = r.get('source') or {}
        au = s.get('first_author') or ''
        return f"{au} {s.get('year', '')}".strip() + (f", doi:{s['doi']}" if s.get('doi') else
                                                      (f", {s.get('url')}" if s.get('url') else ''))

    def evidence(self, rid):
        return (self.rec[rid].get('source') or {}).get('evidence')

    def verified(self, rid):
        v = self.ver.get(rid) or {}
        return bool(v.get('quote_ok'))

    def by(self, **kw):
        out = []
        for r in self.rec.values():
            if all(r.get(k) == v for k, v in kw.items()):
                out.append(r)
        return out
