"""
Survival-conditioned lateral mean squared displacement (MSD) of the first-shell water of every wet SAM system, at the
lags tau = 20, 40, ..., 200 ps.

First shell: a water molecule is in the first shell when its oxygen is within the carbon cutoff of a chain carbon (CH3
or CH2) or within the oxygen cutoff of an oxygen of the silica or of the head groups (cutoffs from step1_RDF.py). Each
face of the slab has its own first shell: a water molecule in the shell is on the face of its closest atom, of the
group with the smallest distance / cutoff.

MSD(tau): over all time origins in the last 90 ns of the production run (frames 20 ps apart), the mean squared
displacement in the plane of the slab (x, y) of the water molecules that are in the first shell of a face at the origin
and stay in it, without a break, until the origin + tau; both faces together. Its uncertainty: the standard deviation
of the MSD of the ten 9-ns blocks of the run (a time origin belongs to the block it is in).

Output: results/MSD/L<L>_G<coverage>.csv: lag_ps, MSD_nm2, MSD_error_nm2
Run:  python3 step2_MSD.py <folder of the simulations>      (after step1_RDF.py)
"""
from collections import deque
from multiprocessing import Pool
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from systems import find_systems, run_folder, simulations_folder

PROCESSES = 8                   # systems analyzed at the same time
LAST_FRAMES = 4501              # the last 90 ns of the run, frames 20 ps apart
LAGS = range(1, 11)             # lags in frames: 20 ... 200 ps
N_BLOCKS = 10                   # the 90 ns in ten blocks of 9 ns
# the cutoffs of the first shell, in Angstrom (the length unit of MDAnalysis)
cutoffs = pd.read_csv('results/RDF/cutoffs.csv', index_col='atoms').cutoff_nm
CARBON_CUTOFF = 10 * float(cutoffs['carbon'])
OXYGEN_CUTOFF = 10 * float(cutoffs['oxygen'])


def closest_atom(xyz, atoms, box, cutoff):
    """For every water oxygen: the distance to the closest atom of the group (Angstrom; infinite if none is within the
    cutoff) and the number of that atom in the group."""
    pairs, distances = capped_distance(xyz, atoms.positions, max_cutoff=cutoff, box=box)
    order = np.lexsort((distances, pairs[:, 0]))           # the pairs by water and, for every water, by distance
    waters, first = np.unique(pairs[order, 0], return_index=True)
    closest = order[first]                                 # the pair of every water with its closest atom
    distance = np.full(len(xyz), np.inf, dtype=np.float32)
    distance[waters] = distances[closest]
    atom = np.zeros(len(xyz), dtype=int)
    atom[waters] = pairs[closest, 1]
    return distance, atom


def shell_faces(xyz, groups, middle, box):
    """The first shell of every water oxygen: 1 on the upper face of the slab, 2 on the lower face, 0 not in the shell.
    For every group of atoms: the closest atom, its distance divided by the cutoff (with a tolerance for the rounding
    of the numbers) and its face, above or below the middle of the silica; a water in the shell takes the face of the
    group with the smallest distance / cutoff."""
    face = np.zeros(len(xyz), dtype=int)
    smallest = np.full(len(xyz), np.inf)
    for atoms, cutoff in groups:
        distance, atom = closest_atom(xyz, atoms, box, cutoff)
        z = atoms.positions[atom, 2] - middle
        z -= box[2] * np.rint(z / box[2])
        ratio = distance / cutoff
        closer = (ratio <= 1.0 + 1e-7) & (ratio < smallest)
        face[closer] = np.where(z >= 0, 1, 2)[closer]
        smallest[closer] = ratio[closer]
    return face


def msd_of_system(simulations, L, coverage):
    """MSD of one system at every lag, written to results/MSD/L<L>_G<coverage>.csv."""
    folder = run_folder(simulations, 'wet', L, coverage)
    u = mda.Universe(str(folder / 'confout.gro'), str(folder / 'traj_unwrap_whole.xtc'))
    water = u.select_atoms('resname TIP3 and name OH2')                       # the water oxygens
    silica = u.select_atoms('resname SUR')
    groups = [(u.select_atoms(f'resname OTS and name C{L}'), CARBON_CUTOFF),                  # CH3
              (u.select_atoms(f'resname OTS and name C* and not name C{L}'), CARBON_CUTOFF),  # CH2
              (u.select_atoms('(resname SUR or resname OTS) and name O*'), OXYGEN_CUTOFF)]     # silica and head-group O
    frames = u.trajectory[-LAST_FRAMES:]

    # sums over the time origins of the squared lateral displacements of the surviving water molecules, for each face
    # of the slab, block (that of the origin) and lag; and the number of these molecules, for each block and lag
    squared_sum = np.zeros((2, N_BLOCKS, len(LAGS)))
    n_survivors = np.zeros((N_BLOCKS, len(LAGS)), dtype=int)
    times = []                                     # time of every frame (ns)
    unwrapped_xy = np.zeros((len(water), 2))       # lateral positions without the jumps over the box edges
    history = deque(maxlen=max(LAGS) + 1)          # the lateral positions of this frame and of the 10 before
    face = np.zeros(len(water), dtype=int)
    run_length = np.zeros(len(water), dtype=int)   # frames in the shell of the same face without a break, up to now

    for j, ts in enumerate(frames):
        times.append(ts.time / 1000.0)
        box = ts.dimensions.astype(np.float64)
        xyz = water.positions

        # lateral positions: add the step from the frame before, taken across the nearest box edge
        xy = xyz[:, :2].astype(np.float64)
        if j > 0:
            step = xy - xy_before
            step -= box[:2] * np.rint(step / box[:2])
            unwrapped_xy = unwrapped_xy + step
        xy_before = xy
        history.append(unwrapped_xy)

        # first shell; the frames in it without a break: one more if the water stays, 1 if it just came in
        face_before = face
        face = shell_faces(xyz, groups, silica.positions[:, 2].mean(), box)
        stayed = (face != 0) & (face == face_before)
        run_length = np.where(stayed, run_length + 1, np.where(face != 0, 1, 0))

        # this frame is the end of the lags from the origins 1 ... 10 frames before
        for i, lag in enumerate(LAGS):
            if j < lag:
                continue
            block = min(N_BLOCKS - 1, (j - lag) * N_BLOCKS // len(frames))
            for side in [1, 2]:
                survivors = (run_length >= lag + 1) & (face == side)
                squared = np.sum((unwrapped_xy[survivors] - history[-1 - lag][survivors]) ** 2, axis=1) * 0.01  # nm^2
                squared_sum[side - 1, block, i] += squared.sum()
                n_survivors[block, i] += np.count_nonzero(survivors)

    # MSD, both faces together, and the SD of the MSD of the blocks
    dt = float(np.median(np.diff(times)))          # time between the frames (ns)
    rows = []
    for i, lag in enumerate(LAGS):
        block_sums = squared_sum[:, :, i].sum(axis=0)
        rows.append({'lag_ps': lag * dt * 1000, 'MSD_nm2': block_sums.sum() / n_survivors[:, i].sum(),
                     'MSD_error_nm2': np.std(block_sums / n_survivors[:, i], ddof=1)})
    pd.DataFrame(rows).to_csv(f'results/MSD/L{L}_G{coverage}.csv', index=False)


if __name__ == '__main__':
    simulations = simulations_folder()
    Path('results/MSD').mkdir(parents=True, exist_ok=True)
    systems = [(simulations, L, coverage) for L, coverage in find_systems(simulations, 'wet')]
    with Pool(PROCESSES) as pool:
        pool.starmap(msd_of_system, systems)
