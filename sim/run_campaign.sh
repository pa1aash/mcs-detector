#!/usr/bin/env bash
# run_campaign.sh --- S3 Highland-validation + kappa_M-linearity campaign.
# Runs the locked-physics mcs_sim over (energy x thickness) on a homogeneous PLA
# slab, plus empty-frame (vacuum) baselines. Run inside the g4highland env:
#   conda activate g4highland && bash sim/run_campaign.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/sim/build/mcs_sim"
OUT="$ROOT/data/runs"
THREADS="${THREADS:-6}"
NEVT="${NEVT:-500000}"        # solid runs
NEVT_EMPTY="${NEVT_EMPTY:-200000}"
MAXSTEP="${MAXSTEP:-0.1}"     # mm (locked)
CUT="${CUT:-0.05}"           # mm (locked)
ENERGIES="${ENERGIES:-200 500 1000}"
THICKS="${THICKS:-3 4 8 16}"  # mm; {3,8,16}=validation+linearity, 4=areal control

mkdir -p "$OUT"
[ -x "$BIN" ] || { echo "ERROR: build mcs_sim first ($BIN missing)"; exit 1; }

run_one() {  # mat thick energy nevt tag
  local mat="$1" thick="$2" energy="$3" nevt="$4" tag="$5"
  local mac; mac="$(mktemp)"
  cat > "$mac" <<EOF
/mcs/det/material $mat
/mcs/det/thickness $thick mm
/mcs/det/maxStep $MAXSTEP mm
/run/setCut $CUT mm
/run/initialize
/gun/particle proton
/gun/energy $energy MeV
/analysis/setFileName $OUT/$tag
/run/printProgress 250000
/run/beamOn $nevt
EOF
  echo ">>> $tag  (mat=$mat thick=${thick}mm E=${energy}MeV N=$nevt)"
  "$BIN" "$mac" "$THREADS" > "$OUT/$tag.run.log" 2>&1
  rm -f "$mac"
}

t0=$SECONDS
for E in $ENERGIES; do
  for T in $THICKS; do
    run_one PLA "$T" "$E" "$NEVT" "solid_E${E}_t${T}"
  done
  run_one G4_Galactic 16 "$E" "$NEVT_EMPTY" "empty_E${E}"
done
echo "=== campaign done in $((SECONDS - t0))s ==="
ls -1 "$OUT"/*.root | sed "s#$OUT/##"
