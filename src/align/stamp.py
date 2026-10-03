"""Provenance stamp: the commit whose inputs produced the numbers."""
import subprocess

from .data import ROOT

OUTPUTS = ('results/', 'FINDINGS.md', 'REPORT.md', 'COMPLETE.md')


def stamp():
    def git(*a):
        return subprocess.run(['git', *a], cwd=ROOT, capture_output=True, text=True).stdout
    head = git('rev-parse', '--short=12', 'HEAD').strip()
    dirty = [l[3:] for l in git('status', '--porcelain').splitlines()
             if l.strip() and not l[3:].startswith(OUTPUTS)]
    return {'commit': head, 'clean_inputs': not dirty, 'dirty_files': dirty}
