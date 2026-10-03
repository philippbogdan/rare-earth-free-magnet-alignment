#!/bin/sh
# Regenerates every number in results/, FINDINGS.md and REPORT.md from the committed inputs.
#   ./regenerate.sh            offline (uses data/verification.json as committed)
#   ./regenerate.sh --verify   also re-checks every source and quote online first
set -e
cd "$(dirname "$0")"
if [ "$1" = "--verify" ]; then
  .venv/bin/python tools/verify_sources.py data/raw/*.yaml data/curated/records*.yaml $(ls data/outcomes/*.yaml 2>/dev/null)
fi
PYTHONPATH=src .venv/bin/python -m align.cli
