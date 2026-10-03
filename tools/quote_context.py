"""Print, for records whose quote was not found automatically, the best-matching passage of the
saved source text, for checking by eye (results go to data/verification_manual.yaml by hand).

    .venv/bin/python tools/quote_context.py <record id> ...
"""
import difflib
import hashlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, os.path.dirname(__file__))
from align.data import Registry                       # noqa: E402
import verify_sources as V                             # noqa: E402


def source_texts(r):
    src = r.get('source') or {}
    out = []
    if src.get('url'):
        key = hashlib.sha1(src['url'].encode()).hexdigest()
        for ext in ('.pdftxt', '.txt'):
            p = os.path.join(V.CACHE, key + ext)
            if os.path.exists(p):
                out.append(('cache' + ext, open(p, encoding='utf-8', errors='ignore').read()))
    doi = (src.get('doi') or '').lower().replace('/', '_')
    au = V.norm((src.get('first_author') or '').split(',')[0].split(' ')[-1])
    yr = str(src.get('year') or '')
    for rel, raw, nt in V.local_texts():
        name = rel.lower()
        if (doi and doi in name) or (au and len(au) > 2 and au in V.norm(name) and yr in name):
            out.append((rel, raw))
    return out


def best_window(quote, text, pad=60):
    t = ' '.join(text.split())
    q = ' '.join((quote or '').split())
    sm = difflib.SequenceMatcher(None, q.lower(), t.lower(), autojunk=False)
    b = sm.find_longest_match(0, len(q), 0, len(t))
    lo = max(0, b.b - b.a - pad)
    return t[lo: lo + len(q) + 2 * pad]


if __name__ == '__main__':
    reg = Registry()
    for rid in sys.argv[1:]:
        r = reg.rec[rid]
        print('=' * 100)
        print(rid, '|', (r.get('source') or {}).get('url'))
        print('QUOTE :', ' '.join(str(r.get('quote')).split())[:700])
        for name, text in source_texts(r)[:2]:
            print('SOURCE', name, ':', best_window(r.get('quote'), text)[:900])
