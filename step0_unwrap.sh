#!/bin/bash
# Whole molecules for the analysis: every production run gets traj_unwrap_whole.xtc, its trajectory centered on the
# whole system (gmx trjconv -center) and then with every molecule made whole across the periodic boundaries
# (-pbc whole). The centered trajectory is kept as traj_center.xtc.
#
# Run:  bash step0_unwrap.sh <folder of the simulations>
set -e
if [ -z "$1" ]; then echo "Run:  bash step0_unwrap.sh <folder of the simulations>"; exit 1; fi
cd "$1"

for run in all_dry_production_run/chain_length_*/coverage_*/nvtrunP all_wet_production_run/chain_length_*/coverage_*/nptrunP; do
    if [ -f $run/nvt.tpr ]; then tpr=$run/nvt.tpr; else tpr=$run/npt.tpr; fi      # dry runs NVT, wet runs NPT
    printf '0\n0\n' | gmx trjconv -s $tpr -f $run/traj_comp.xtc -o $run/traj_center.xtc -center
    printf '0\n' | gmx trjconv -s $tpr -f $run/traj_center.xtc -o $run/traj_unwrap_whole.xtc -pbc whole
done
