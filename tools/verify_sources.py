"""Verify sources and quotes of every record in data/raw/*.yaml, data/outcomes/*.yaml and
data/curated/*.yaml.

For each record: (1) the DOI resolves on Crossref and its title/year/first author match;
(2) the verbatim quote is found in the abstract (Crossref / OpenAlex / Semantic Scholar)
or in the text at `url` (HTML or PDF via pdftotext). Matching ignores case, whitespace,
punctuation and line-break hyphenation. Results go to data/verification.json.
Fetched texts are cached in .tmp/cache (not committed).

    .venv/bin/python tools/verify_sources.py [files...]
"""
import difflib
import glob
import hashlib
import html
import json
import os
import re
import subprocess
import sys
import time
import unicodedata
import urllib.parse
import urllib.request

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, '.tmp', 'cache')
OUT = os.path.join(ROOT, 'data', 'verification.json')
UA = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 '
                    '(KHTML, like Gecko) Chrome/124 Safari/537.36'}
os.makedirs(CACHE, exist_ok=True)


def _get(url, binary=False, timeout=40):
    key = hashlib.sha1(url.encode()).hexdigest()
    p = os.path.join(CACHE, key + ('.bin' if binary else '.txt'))
    if os.path.exists(p):
        return open(p, 'rb').read() if binary else open(p, encoding='utf-8').read()
    try:
        req = urllib.request.Request(url, headers=UA)
        data = urllib.request.urlopen(req, timeout=timeout).read()
    except Exception as e:  # noqa: BLE001
        return None
    time.sleep(0.2)
    if binary:
        open(p, 'wb').write(data)
        return data
    s = data.decode('utf-8', errors='ignore')
    open(p, 'w', encoding='utf-8').write(s)
    return s


def norm(s):
    s = unicodedata.normalize('NFKC', s or '')
    s = s.replace('­', '').replace('ﬁ', 'fi').replace('ﬂ', 'fl')
    s = re.sub(r'-\s*\n\s*', '', s)
    s = s.lower()
    return re.sub(r'[\W_]+', '', s)


def strip_tags(s):
    s = re.sub(r'<script.*?</script>|<style.*?</style>', ' ', s or '', flags=re.S)
    return html.unescape(re.sub(r'<[^>]+>', ' ', s))


def crossref(doi):
    s = _get('https://api.crossref.org/works/' + urllib.parse.quote(doi))
    try:
        return json.loads(s)['message']
    except Exception:  # noqa: BLE001
        return None


def abstracts(doi):
    texts = []
    m = crossref(doi)
    if m and m.get('abstract'):
        texts.append(strip_tags(m['abstract']))
    s = _get('https://api.openalex.org/works/doi:' + urllib.parse.quote(doi))
    try:
        inv = json.loads(s).get('abstract_inverted_index')
        if inv:
            pos = {}
            for w, ix in inv.items():
                for i in ix:
                    pos[i] = w
            texts.append(' '.join(pos[i] for i in sorted(pos)))
    except Exception:  # noqa: BLE001
        pass
    s = _get('https://api.semanticscholar.org/graph/v1/paper/DOI:' + urllib.parse.quote(doi)
             + '?fields=title,abstract')
    try:
        a = json.loads(s).get('abstract')
        if a:
            texts.append(a)
    except Exception:  # noqa: BLE001
        pass
    return m, texts


def url_text(url):
    """Text at a URL; PDFs are recognised by their content, not their URL."""
    if not url:
        return ''
    key = hashlib.sha1(url.encode()).hexdigest()
    pdftxt = os.path.join(CACHE, key + '.pdftxt')
    if os.path.exists(pdftxt):                  # PDFs are kept as extracted text only (disk)
        return open(pdftxt, encoding='utf-8').read()
    data = _get(url, binary=True)
    if not data:
        return ''
    if data[:5] == b'%PDF-':
        pdf = os.path.join(CACHE, key + '.bin')
        r = subprocess.run(['pdftotext', '-q', pdf, '-'], capture_output=True, text=True)
        open(pdftxt, 'w', encoding='utf-8').write(r.stdout)
        os.remove(pdf)
        return r.stdout
    return strip_tags(data.decode('utf-8', 'ignore'))


_LOCAL = None


def local_texts():
    """Saved source texts from the collection work (.tmp/**/*.txt), normalised once."""
    global _LOCAL
    if _LOCAL is None:
        _LOCAL = []
        for f in glob.glob(os.path.join(ROOT, '.tmp', '**', '*.txt'), recursive=True):
            if os.path.dirname(f) == os.path.join(ROOT, '.tmp'):
                continue                      # our own listings and printouts, not sources
            try:
                if os.path.getsize(f) > 30_000_000:
                    continue
                t = open(f, encoding='utf-8', errors='ignore').read()
                _LOCAL.append((os.path.relpath(f, ROOT), t, norm(t)))
            except OSError:
                pass
    return _LOCAL


def find_quote(quote, text):
    """Score in [0,1]; quotes with '...' elisions are checked fragment by fragment
    (fragments shorter than 12 normalised characters are ignored); the score is the
    length-weighted mean of the fragment scores."""
    frags = [f for f in re.split(r'\.\.\.|…|\[\.\.\.\]', quote or '') if len(norm(f)) >= 12]
    if len(frags) > 1:
        t = norm(text)
        tot = sum(len(norm(f)) for f in frags)
        return sum(_find_one(f, t, normed=True) * len(norm(f)) for f in frags) / tot
    return _find_one(quote, text)


def _find_one(quote, text, normed=False):
    q, t = norm(quote), (text if normed else norm(text))
    if not q or not t:
        return 0.0
    if q in t:
        return 1.0
    sm = difflib.SequenceMatcher(None, q, t, autojunk=False)
    b = sm.find_longest_match(0, len(q), 0, len(t))
    lo = max(0, b.b - b.a - 20)
    win = t[lo: lo + len(q) + 40]
    return difflib.SequenceMatcher(None, q, win, autojunk=False).ratio()


NUM = re.compile(r'(?<![\w.])[-−]?\d+(?:[.,]\d+)?')


def numbers(s):
    return [n.replace('−', '-').replace(',', '.') for n in NUM.findall(unicodedata.normalize('NFKC', s or ''))]


def numbers_found(quote, text):
    """Fraction of the numbers in the quote that occur verbatim (as digit strings) in the
    source text; None if the quote has no numbers."""
    ns = [n for n in numbers(quote) if len(n.strip('-').replace('.', '')) >= 1]
    if not ns:
        return None
    t = unicodedata.normalize('NFKC', text or '').replace('−', '-')
    t = re.sub(r'\s+', ' ', t)
    hit = 0
    for n in ns:
        variants = {n, n.replace('.', ','), n.replace('.', '·')}
        if any(re.search(r'(?<![\d.])' + re.escape(v) + r'(?![\d])', t) for v in variants):
            hit += 1
    return hit / len(ns)


def value_in_quote(value, quote):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    qs = numbers(quote)
    for q in qs:
        try:
            if abs(float(q) - float(value)) <= 1e-9 * max(1, abs(float(value))):
                return True
        except ValueError:
            pass
    return False


def title_score(a, b):
    return difflib.SequenceMatcher(None, norm(a), norm(b)).ratio()


def verify_record(r):
    src = r.get('source') or {}
    doi = (src.get('doi') or '').strip()
    out = {'file': r.get('_file'), 'doi': doi or None}
    quote = r.get('quote') or ''
    texts = []
    if doi:
        m, texts = abstracts(doi)
        if m:
            t = (m.get('title') or [''])[0]
            out['crossref_title'] = t
            out['title_score'] = round(title_score(t, src.get('title', '')), 2)
            y = (m.get('issued') or {}).get('date-parts', [[None]])[0][0]
            out['crossref_year'] = y
            out['year_ok'] = (y is None or src.get('year') is None
                              or abs(int(y) - int(src.get('year'))) <= 1)
            fam = [(a.get('family') or '') for a in m.get('author', [])]
            fa = (src.get('first_author') or '')
            out['author_ok'] = bool(fa) and any(norm(f) and norm(f) in norm(fa) or
                                                 norm(fa) in norm(f) for f in fam if f)
            out['doi_ok'] = True
        else:
            out['doi_ok'] = False
    best, where, btext = 0.0, None, ''
    for t in texts:
        s = find_quote(quote, t)
        if s > best:
            best, where, btext = s, 'abstract', t
    if best < 0.95 and src.get('url'):
        ut = url_text(src['url'])
        s = find_quote(quote, ut)
        if s > best:
            best, where, btext = s, 'url', ut
    if best < 0.9:
        q = norm(quote)
        frags = [norm(f) for f in re.split(r'\.\.\.|…', quote or '') if len(norm(f)) >= 12]
        for rel, raw, nt in local_texts():
            if frags and all(f in nt for f in frags):
                best, where, btext = 1.0, 'local:' + rel, raw
                break
    if best < 0.9:                    # fuzzy match in the saved text of this record's own source
        slug = doi.lower().replace('/', '_') if doi else None
        au = norm((src.get('first_author') or '').split(',')[0].split(' ')[-1]) if src.get('first_author') else ''
        yr = str(src.get('year') or '')
        for rel, raw, nt in local_texts():
            name = rel.lower()
            if (slug and slug in name) or (au and len(au) > 2 and au in norm(name) and yr and yr in name):
                sc = find_quote(quote, raw)
                if sc > best:
                    best, where, btext = sc, 'local:' + rel, raw
    out['quote_score'] = round(best, 3)
    out['quote_where'] = where
    out['level'] = (None if best < 0.9 else ('local' if (where or '').startswith('local:') else 'online'))
    nf = numbers_found(quote, btext) if btext else None
    out['numbers_found'] = None if nf is None else round(nf, 3)
    out['value_in_quote'] = value_in_quote(r.get('value'), quote)
    out['quote_ok'] = best >= 0.9 and (nf is None or nf >= 0.999)
    return out


def main(files):
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for fn in files:
        recs = yaml.safe_load(open(fn)) or []
        for r in recs:
            r['_file'] = os.path.relpath(fn, ROOT)
            rid = r.get('id')
            if not rid:
                continue
            res[rid] = verify_record(r)
            st = res[rid]
            print(f"{rid:50s} doi={st.get('doi_ok')} title={st.get('title_score')} "
                  f"quote={st['quote_score']} ({st['quote_where']})")
        json.dump(res, open(OUT, 'w'), indent=1, sort_keys=True)


def recheck_failed(files):
    """Re-verify only the records whose quote was not found (e.g. after adding saved texts)."""
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for fn in files:
        for r in yaml.safe_load(open(fn)) or []:
            rid = r.get('id')
            if not rid or (res.get(rid) or {}).get('quote_ok'):
                continue
            r['_file'] = os.path.relpath(fn, ROOT)
            res[rid] = verify_record(r)
            print(f"{rid:50s} quote={res[rid]['quote_score']} ({res[rid]['quote_where']})", flush=True)
    json.dump(res, open(OUT, 'w'), indent=1, sort_keys=True)


if __name__ == '__main__':
    args = sys.argv[1:]
    if args and args[0] == '--recheck-failed':
        recheck_failed(args[1:])
    else:
        main(args or sorted(glob.glob(os.path.join(ROOT, 'data', 'raw', '*.yaml'))))
