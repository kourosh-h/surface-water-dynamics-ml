"""
Water of the first shell at each group of atoms of every wet SAM system: the terminal carbons (CH3), the other chain
carbons (CH2) and the oxygens of the silica and of the head groups (oxygen).

First shell: a water molecule is in the first shell when its oxygen is within the carbon cutoff of a chain carbon or
within the oxygen cutoff of an oxygen (cutoffs from step1_RDF.py). Every water molecule of the shell is given to the
group of its closest atom, by the distance divided by the cutoff of the group. For every group, the number of these
water molecules per nm^2 of each face of the slab, averaged over each 9-ns block of the last 90 ns of the production
run (frames 20 ps apart).

Output: results/first_shell_water/L<L>_G<coverage>.csv: block, group, water_per_nm2
Run:  python3 step4_first_shell_water.py <folder of the simulations>      (after step1_RDF.py)
"""
from multiprocessing import Pool
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from systems import find_systems, run_folder, simulations_folder

PROCESSES = 8                   # systems analyzed at the same time
LAST_FRAMES = 4501              # the last 90 ns of the run, frames 20 ps apart
N_BLOCKS = 10                   # the 90 ns in ten blocks of 9 ns
GROUP_NAMES = ['CH3', 'CH2', 'oxygen']
# the cutoffs of the first shell, in Angstrom (the length unit of MDAnalysis)
cutoffs = pd.read_csv('results/RDF/cutoffs.csv', index_col='atoms').cutoff_nm
CARBON_CUTOFF = 10 * float(cutoffs['carbon'])
OXYGEN_CUTOFF = 10 * float(cutoffs['oxygen'])


def closest_distance(xyz, atoms, box, cutoff):
    """Distance from every water oxygen to the closest atom of the group (Angstrom); infinite if none is within the
    cutoff."""
    distance = np.full(len(xyz), np.inf, dtype=np.float32)
    pairs, distances = capped_distance(xyz, atoms.positions, max_cutoff=cutoff, box=box)
    np.minimum.at(distance, pairs[:, 0], distances)
    return distance


def first_shell_water_of_system(simulations, L, coverage):
    """First-shell water at each group, per block, of one system, written to results/first_shell_water/."""
    folder = run_folder(simulations, 'wet', L, coverage)
    u = mda.Universe(str(folder / 'confout.gro'), str(folder / 'traj_unwrap_whole.xtc'))
    water = u.select_atoms('resname TIP3 and name OH2')                       # the water oxygens
    groups = [(u.select_atoms(f'resname OTS and name C{L}'), CARBON_CUTOFF),                  # CH3
              (u.select_atoms(f'resname OTS and name C* and not name C{L}'), CARBON_CUTOFF),  # CH2
              (u.select_atoms('(resname SUR or resname OTS) and name O*'), OXYGEN_CUTOFF)]     # silica and head-group O

    counts = [[], [], []]               # first-shell water at each group, in every frame
    area = []                           # area of a face of the slab in every frame (nm^2)
    for ts in u.trajectory[-LAST_FRAMES:]:
        box = ts.dimensions.astype(np.float64)
        # distance to the closest atom of each group divided by its cutoff (with a tolerance for the rounding of the
        # numbers in the test of the shell)
        ratios = []
        for atoms, cutoff in groups:
            ratios.append(closest_distance(water.positions, atoms, box, cutoff) / cutoff)
        ratios = np.array(ratios)
        in_shell = (ratios <= 1.0 + 1e-7).any(axis=0)
        closest = np.argmin(ratios, axis=0)
        for g in range(3):
            counts[g].append(np.count_nonzero(in_shell & (closest == g)))
        area.append(box[0] * box[1] / 100.0)

    # per nm^2 of each of the two faces, averaged over each block of frames
    rows = []
    for b, frames in enumerate(np.array_split(np.arange(len(area)), N_BLOCKS)):
        for g, name in enumerate(GROUP_NAMES):
            per_nm2 = np.array(counts[g])[frames] / np.array(area)[frames] / 2.0
            rows.append({'block': b + 1, 'group': name, 'water_per_nm2': np.mean(per_nm2)})
    pd.DataFrame(rows).to_csv(f'results/first_shell_water/L{L}_G{coverage}.csv', index=False)


if __name__ == '__main__':
    simulations = simulations_folder()
    Path('results/first_shell_water').mkdir(parents=True, exist_ok=True)
    systems = [(simulations, L, coverage) for L, coverage in find_systems(simulations, 'wet')]
    with Pool(PROCESSES) as pool:
        pool.starmap(first_shell_water_of_system, systems)
