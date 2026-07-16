#!/usr/bin/env bash
# S6 Stage-1 voronoi escalation -> 3e7, to settle real-foam-excess vs deconvolution
# noise in the small-signal (high-N_eff) configs. Resumable (driver skips >=target).
set -e
cd "$(dirname "$0")/.."
source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh
conda activate g4highland
DRV="python sim/campaign.py --energies 500 --threads 6 --nevt 3e7"
for f in 0.2 0.3 0.4 0.5; do
  echo "=== voronoi f$f -> 3e7   $(date +%T) ==="
  $DRV --only voronoi:${f}:500
done
echo "=== VORONOI 3e7 DONE $(date +%T) ==="
