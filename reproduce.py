"""Run offline checks, regenerate results and compare their saved SHA-256 values.

Usage: python reproduce.py
       python reproduce.py --write-expected

The second form hashes existing results without running or checking them.
See REPRODUCE.md for dependencies, expected output and separate source checks.
"""
import glob
import hashlib
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, 'src'))
EXPECTED = os.path.join(ROOT, 'tests', 'expected.json')
FILES = ['results/account.csv'] + [f'results/scores_v{i}.json' for i in range(1, 6)] + ['results/summary.json']


def sha(path):
    with open(os.path.join(ROOT, path), 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def run_tests():
    """Stop at the first failing geometry, physics or route check."""
    for t in sorted(glob.glob(os.path.join(ROOT, 'tests', 'test_*.py'))):
        r = subprocess.run([sys.executable, t], capture_output=True, text=True)
        if r.returncode:
            print(r.stdout, r.stderr)
            sys.exit(f'test failed: {os.path.relpath(t, ROOT)}')
    print('tests passed')


def headline(s):
    """Print the account scope and each version's fresh held-out result."""
    a, m, v = s['account'], s['mnal_c_textured_parent'], s['validation']
    print(f"account: {a['cells']} cells ({a['phases']} phases x {a['route_classes']} routes, including Nd2Fe14B reference): "
          + ', '.join(f'{k} {n}' for k, n in a['verdicts'].items()))
    print(f"MnAl(C) from a perfectly textured parent: easy axes at {m['easy_axis_angle_to_fibre_deg']} deg, "
          f"A_max {m['A_max_AM']} against {m['A_random']} for random grains")
    for ver in ('v1', 'v2', 'v3', 'v4', 'v5'):
        x = v[ver]
        print(f"{ver} held-out: {x['heldout_measurements']} measurements, within {x['within']}, "
              f"straddling {x['straddling']}, violated {x['violated']}")


def main():
    """Regenerate and check results, or explicitly replace the checksum baseline."""
    if '--write-expected' in sys.argv:
        with open(EXPECTED, 'w') as f:
            json.dump({p: sha(p) for p in FILES}, f, indent=1)
            f.write('\n')
        print('wrote', os.path.relpath(EXPECTED, ROOT))
        return
    run_tests()
    from align import cli
    s = cli.main()
    headline(s)
    want = json.load(open(EXPECTED))
    bad = [p for p in FILES if sha(p) != want.get(p)]
    if bad:
        print('MISMATCH:', ', '.join(bad))
        sys.exit(1)
    print('MATCH: every regenerated file is identical to the expected one')


if __name__ == '__main__':
    main()
