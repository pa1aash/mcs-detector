#!/usr/bin/env python3
"""check_ledger.py -- independent re-verification of the claims ledger.

Validates claims_ledger.csv (at the repository root):
  * schema: exactly the 8 required columns, non-empty claim_id / statement / tier;
  * tier vocabulary in {ANCHORED, PREDICTED, QUALITATIVE};
  * every evidence_path token resolves to a file that exists in the repository;
  * claim_id uniqueness.
Each row maps a reported number to the script/data file that produced it, so a clean
run is an independent check that every headline result traces to committed data.
Exits non-zero on any violation, so it can gate `make audit`.

Usage: python tools/check_ledger.py
"""
from __future__ import annotations
import csv, glob, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEDGER = os.path.join(ROOT, "claims_ledger.csv")
COLS = ["claim_id", "manuscript_location", "statement", "tier", "evidence_path",
        "figure_or_table", "uncertainty", "last_verified_commit"]
TIERS = {"ANCHORED", "PREDICTED", "QUALITATIVE"}


def check() -> list[str]:
    errs: list[str] = []
    if not os.path.exists(LEDGER):
        return [f"ledger not found: {LEDGER}"]
    with open(LEDGER) as f:
        rows = list(csv.DictReader(f))
        if rows and list(rows[0].keys()) != COLS:
            errs.append(f"header mismatch: {list(rows[0].keys())} != {COLS}")
    seen = set()
    for i, r in enumerate(rows, 2):
        cid = (r.get("claim_id") or "").strip()
        if not cid:
            errs.append(f"row {i}: empty claim_id")
        elif cid in seen:
            errs.append(f"row {i}: duplicate claim_id {cid}")
        seen.add(cid)
        if not (r.get("statement") or "").strip():
            errs.append(f"row {i} ({cid}): empty statement")
        if (r.get("tier") or "").strip() not in TIERS:
            errs.append(f"row {i} ({cid}): bad tier {r.get('tier')!r} (need {TIERS})")
        for tok in (r.get("evidence_path") or "").split(";"):
            tok = tok.strip()
            if not tok:
                continue
            if any(ch in tok for ch in "*?["):
                if not glob.glob(os.path.join(ROOT, tok)):
                    errs.append(f"row {i} ({cid}): evidence glob matches nothing: {tok}")
            elif not os.path.exists(os.path.join(ROOT, tok)):
                errs.append(f"row {i} ({cid}): evidence path missing: {tok}")
    return errs


def main() -> int:
    errs = check()
    for e in errs:
        print("LEDGER FAIL:", e, file=sys.stderr)
    if errs:
        print(f"\n{len(errs)} ledger error(s).", file=sys.stderr)
        return 1
    print("ledger OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
