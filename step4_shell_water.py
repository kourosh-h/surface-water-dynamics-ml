"""
Step 4: the water of the first shell at each group of atoms of every system.

First shell: as in step 2. Every water molecule of the shell is given to the group of its closest atom, by the
distance divided by the contact distance of the group. For every group, the number of these water molecules per nm^2
of each analyzed face, averaged over each block of the analyzed frames.

Output: results/shell_water/<system>.csv: block, group, water_per_nm2
Run:  python3 step4_shell_water.py      (after step 1)
"""
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import settings
from first_shell import faces_and_groups, open_system


def shell_water_of_system(system):
    """First-shell water at each group, per block, of one system, written to results/shell_water/<system>.csv."""
    u, water, solid, groups = open_system(system)
    counts = []                         # first-shell water at each group, in every frame
    for g in range(len(groups)):
        counts.append([])
    area = []                           # area of a face in every frame (nm^2)
    for ts in u.trajectory[-settings.LAST_FRAMES:]:
        box = ts.dimensions.astype(np.float64)
        face, group = faces_and_groups(water.positions, groups, solid.positions[:, 2].mean(), box)
        for g in range(len(groups)):
            counts[g].append(np.count_nonzero(group == g))
        area.append(box[0] * box[1] / 100.0)

    # per nm^2 of each analyzed face, averaged over each block of frames
    rows = []
    for b, frames in enumerate(np.array_split(np.arange(len(area)), settings.N_BLOCKS)):
        for g, (name, atoms, cutoff) in enumerate(groups):
            per_nm2 = np.array(counts[g])[frames] / np.array(area)[frames] / len(settings.FACES)
            rows.append({'block': b + 1, 'group': name, 'water_per_nm2': np.mean(per_nm2)})
    pd.DataFrame(rows).to_csv(f'results/shell_water/{system["name"]}.csv', index=False)


if __name__ == '__main__':
    Path('results/shell_water').mkdir(parents=True, exist_ok=True)
    with Pool(settings.PROCESSES) as pool:
        pool.map(shell_water_of_system, settings.SYSTEMS)
