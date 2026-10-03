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


def seal(batch):
    os.makedirs(DST, exist_ok=True)
    man = _manifest()
    sealed_files = {f for b in man['batches'].values() for f in b['files']}
    entry = man['batches'].setdefault(batch, {'files': {}, 'opened': False})
    for fn in sorted(os.listdir(SRC)):
        if not fn.endswith('.yaml') or fn in sealed_files:
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
        entry['files'][fn] = {'sha256': hashlib.sha256(raw).hexdigest(),
                              'n_records': len(recs), 'cases': cases}
        print(f'sealed {fn}: {len(recs)} records')
    json.dump(man, open(MAN, 'w'), indent=1)


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
    {'seal': seal, 'open': open_batch}[cmd](batch)
