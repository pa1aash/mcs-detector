#!/usr/bin/env bash
# S6 Stage-2B: 200 MeV campaign at full statistics, tiered to resolve the deconvolved
# Delta_kappa4 to <30% CI (same tiers proven at 500 MeV in Stage 1). Resumable.
set -e
cd "$(dirname "$0")/.."
source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh
conda activate g4highland
DRV="python sim/campaign.py --energies 200 --threads 6"

echo "=== WAVE 0: solids+empty -> 1e7  $(date +%T) ==="
$DRV --baselines-only --nevt 1e7

for tf in "rectilinear 3e6" "schwarzp 3e6" "gyroid 6e6" "voronoi 3e7" "diamond 1e6"; do
  set -- $tf; topo=$1; nevt=$2
  for f in 0.2 0.3 0.4 0.5; do
    echo "=== $topo f$f -> $nevt   $(date +%T) ==="
    $DRV --only ${topo}:${f}:200 --nevt $nevt
  done
done
echo "=== 200 MeV CAMPAIGN DONE $(date +%T) ==="
