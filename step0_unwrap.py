"""
Step 0: whole molecules in the trajectories (GROMACS). For every system, the trajectory is centered on the whole
system (gmx trjconv -center), and then every molecule is made whole across the periodic boundaries (-pbc whole). The
centered trajectory is kept as traj_center.xtc.

Run:  python3 step0_unwrap.py
"""
import subprocess
from pathlib import Path

import settings

for system in settings.SYSTEMS:
    folder = Path(settings.SIMULATIONS) / system['folder']
    run_input = str(folder / settings.RUN_INPUT)
    centered = str(folder / 'traj_center.xtc')
    # the answers to the questions of gmx trjconv: group 0 (System) to center on and to write
    subprocess.run(['gmx', 'trjconv', '-s', run_input, '-f', str(folder / settings.RAW_TRAJECTORY), '-o', centered,
                    '-center'], input='0\n0\n', text=True, check=True)
    subprocess.run(['gmx', 'trjconv', '-s', run_input, '-f', centered, '-o', str(folder / settings.TRAJECTORY),
                    '-pbc', 'whole'], input='0\n', text=True, check=True)
