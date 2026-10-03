"""Seal outcome files: encrypt .tmp/sealed/<name>.yaml into data/sealed/<name>.yaml.enc and
record SHA-256 of the plaintext plus a value-free case list in data/sealed/MANIFEST.json.

    .venv/bin/python tools/seal.py seal <batch>     # encrypt every unsealed file in .tmp/sealed
    .venv/bin/python tools/seal.py open <batch>     # decrypt a batch into data/outcomes/, check hashes

The key lives in .tmp/sealed.key (not committed) until the batch is opened; opening
commits the key-free plaintext and the hash check proves it is unchanged.
Values are never printed by this tool.
"""
import hashlib
import json
import os
import secrets
import subprocess
import sys

import re

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, '.tmp', 'sealed')
DST = os.path.join(ROOT, 'data', 'sealed')
OPEN = os.path.join(ROOT, 'data', 'outcomes')
KEY = os.path.join(ROOT, '.tmp', 'sealed.key')
MAN = os.path.join(DST, 'MANIFEST.json')


def _key():
    if not os.path.exists(KEY):
        with open(KEY, 'w') as f:
            f.write(secrets.token_hex(32))
        os.chmod(KEY, 0o600)
    return KEY


def _manifest():
    if os.path.exists(MAN):
        return json.load(open(MAN))
    return {'batches': {}}


REDACT = re.compile(r'(Br\b|B_r|Mr\b|M_r|Mr/Ms|remanen\w*|squareness|\(BH\)\s*max|BHmax|energy product'
                    r'|pole[- ]figure max\w*|MRD|mrd|Lotgering|orientation factor|alignment degree'
                    r'|texture coefficient|variant fraction)[^;\n]*', re.I)
COND_KEYS = ('id', 'phase', 'route', 'process', 'conditions', 'quantity', 'unit',
             'measured_along', 'measurement', 'n_samples')


def _redact(v):
    if isinstance(v, str):
        return REDACT.sub('[redacted]', v)
    if isinstance(v, dict):
        return {k: _redact(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_redact(x) for x in v]
    return v


COND_ALLOW = re.compile(r'^(field_T|field_Oe|field_kOe|field_direction|T_process_K|T_[a-z_]*K|time_.*|'
                        r'stress_MPa|stress_axis|stress_sign|theta1_deg|theta2_deg|strain_percent|'
                        r'parent_texture|parent_phase.*|magnet_axis|composition.*|particle_size.*|substrate|'
                        r'Tc_relation|la_co|route_detail|alignment.*|consolidation.*|measured_direction|'
                        r'density_rel|phase_content|packing_fraction)$')


def _allowed_conditions(c):
    """Process parameters only: any key not on the allow-list is dropped (outcome-related
    numbers such as Mr, M, Br, (BH)max, Hc, Ku placed in conditions never pass)."""
    if not isinstance(c, dict):
        return _redact(c)
    return {k: _redact(v) for k, v in c.items() if COND_ALLOW.match(str(k))}


def conditions_only(recs):
    """Value-free view of outcome records (no value, uncertainty, quote or note)."""
    out = []
    for r in recs:
        c = {k: _redact(r.get(k)) for k in COND_KEYS if k in r and k != 'conditions'}
        c['conditions'] = _allowed_conditions(r.get('conditions'))
        src = r.get('source') or {}
        c['source'] = {k: src.get(k) for k in ('doi', 'title', 'first_author', 'year',
                                               'evidence', 'locator', 'url') if k in src}
        out.append(c)
    return out


def seal(batch, only=None):
    """Seal unsealed files of .tmp/sealed (only the named ones if `only` is given)."""
    os.makedirs(DST, exist_ok=True)
    man = _manifest()
    sealed_files = {f for b in man['batches'].values() for f in b['files']}
    entry = man['batches'].setdefault(batch, {'files': {}, 'opened': False})
    for fn in sorted(os.listdir(SRC)):
        if not fn.endswith('.yaml') or fn in sealed_files or (only and fn not in only):
            continue
        p = os.path.join(SRC, fn)
        raw = open(p, 'rb').read()
        recs = yaml.safe_load(raw) or []
        cases = [{'id': r.get('id'), 'phase': r.get('phase'), 'route': r.get('route'),
                  'quantity': r.get('quantity'),
                  'doi': (r.get('source') or {}).get('doi')} for r in recs]
        out = os.path.join(DST, fn + '.enc')
        subprocess.run(['openssl', 'enc', '-aes-256-cbc', '-pbkdf2', '-salt', '-in', p,
                        '-out', out, '-pass', 'file:' + _key()], check=True)
        cond = os.path.join(DST, fn.replace('.yaml', '.conditions.yaml'))
        yaml.safe_dump(conditions_only(recs), open(cond, 'w'), sort_keys=False,
                       allow_unicode=True, width=100)
        entry['files'][fn] = {'sha256': hashlib.sha256(raw).hexdigest(),
                              'n_records': len(recs), 'cases': cases}
        print(f'sealed {fn}: {len(recs)} records')
    json.dump(man, open(MAN, 'w'), indent=1)


def export_conditions():
    """Regenerate every value-free conditions file from the sealed plaintexts."""
    for fn in sorted(os.listdir(SRC)):
        if fn.endswith('.yaml'):
            recs = yaml.safe_load(open(os.path.join(SRC, fn))) or []
            cond = os.path.join(DST, fn.replace('.yaml', '.conditions.yaml'))
            yaml.safe_dump(conditions_only(recs), open(cond, 'w'), sort_keys=False,
                           allow_unicode=True, width=100)
            print('exported', cond)


def open_batch(batch):
    man = _manifest()
    entry = man['batches'][batch]
    os.makedirs(OPEN, exist_ok=True)
    for fn, meta in entry['files'].items():
        enc = os.path.join(DST, fn + '.enc')
        out = os.path.join(OPEN, fn)
        subprocess.run(['openssl', 'enc', '-d', '-aes-256-cbc', '-pbkdf2', '-in', enc,
                        '-out', out, '-pass', 'file:' + KEY], check=True)
        h = hashlib.sha256(open(out, 'rb').read()).hexdigest()
        assert h == meta['sha256'], f'hash mismatch for {fn}'
        print(f'opened {fn}: hash ok')
    entry['opened'] = True
    json.dump(man, open(MAN, 'w'), indent=1)


if __name__ == '__main__':
    cmd, batch = sys.argv[1], sys.argv[2]
    if cmd == 'export':
        export_conditions()
    elif cmd == 'seal':
        seal(batch, sys.argv[3:] or None)
    else:
        open_batch(batch)
