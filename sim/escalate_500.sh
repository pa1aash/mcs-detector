#!/usr/bin/env bash
# S6 Stage-1 escalation: tier the 500 MeV configs to resolve the deconvolved
# Delta_kappa4 to <30% CI for the collapse fit. Resumable (driver skips configs
# already >= target). Diamond left at 1e6 (N_eff->inf corner, excluded from fit).
set -e
cd "$(dirname "$0")/.."
source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh
conda activate g4highland
DRV="python sim/campaign.py --energies 500 --threads 6"

echo "=== WAVE 0: solids+empty -> 1e7 (leverages every deconvolution) $(date +%T) ==="
$DRV --baselines-only --nevt 1e7

for tf in "rectilinear 3e6" "schwarzp 3e6" "gyroid 6e6" "voronoi 1e7"; do
  set -- $tf; topo=$1; nevt=$2
  for f in 0.2 0.3 0.4 0.5; do
    echo "=== $topo f$f -> $nevt   $(date +%T) ==="
    $DRV --only ${topo}:${f}:500 --nevt $nevt
  done
done
echo "=== ESCALATION DONE $(date +%T) ==="
