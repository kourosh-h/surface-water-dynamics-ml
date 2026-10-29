"""
Survival-conditioned lateral MSD and contact survival P of the first-shell water of every region of the first shell,
for every wet SAM system, at the lags tau = 20, 40, ..., 200 ps, per 9-ns block.

First shell: a water molecule is in the first shell when its oxygen is within the carbon cutoff of a chain carbon (CH3
or CH2) or within the oxygen cutoff of an oxygen of the silica or of the head groups (cutoffs from step 1 of
../RDF_MSD_and_P). Its group is the one of the three (CH3, CH2, oxygen) with the smallest distance / cutoff, and its
face that of the closest atom of this group.
Regions: all first-shell water; the water whose group is the oxygens; CH2; CH3.
A molecule stays in a region while it is in the first shell, on the same face, and its group belongs to the region.

For every time origin in the last 90 ns of the production run (frames 20 ps apart) and every lag: the molecules in the
region of a face at the origin (origins), those that stay in it, without a break, until the origin + tau (survivors),
and the squared displacement of the survivors in the plane of the slab (x, y); both faces together, summed over the
origins of each of the ten 9-ns blocks (a time origin belongs to the block it is in).

Output: results/per_system/L<L>_G<coverage>_blocks.csv: region, block, lag_ps, MSD_nm2 (of the block), n_survivors,
            n_origins
        results/per_system/L<L>_G<coverage>.csv: region, lag_ps, n_origins, n_survivors, P, MSD_nm2, MSD_error_nm2 (SD
            of the MSD of the blocks with 20 or more survivors, when there are 3 or more of them)
Run:  python3 step1_region_MSD_and_P.py <folder of the simulations>      (after step 1 of ../RDF_MSD_and_P)
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
# the regions: the groups (0 CH3, 1 CH2, 2 oxygen) whose water belongs to them
REGIONS = {'all': (0, 1, 2), 'oxygen': (2,), 'CH2': (1,), 'CH3': (0,)}
# the cutoffs of the first shell, in Angstrom (the length unit of MDAnalysis)
cutoffs = pd.read_csv('../RDF_MSD_and_P/results/RDF/cutoffs.csv', index_col='atoms').cutoff_nm
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


def shell_faces_and_groups(xyz, groups, middle, box):
    """For every water oxygen: its face (1 upper, 2 lower, 0 not in the first shell) and its group (0 CH3, 1 CH2,
    2 oxygen, -1 none): the group with the smallest distance / cutoff (with a tolerance for the rounding of the
    numbers), and the face of its closest atom, above or below the middle of the silica."""
    face = np.zeros(len(xyz), dtype=int)
    group = np.full(len(xyz), -1)
    smallest = np.full(len(xyz), np.inf)
    for g, (atoms, cutoff) in enumerate(groups):
        distance, atom = closest_atom(xyz, atoms, box, cutoff)
        z = atoms.positions[atom, 2] - middle
        z -= box[2] * np.rint(z / box[2])
        ratio = distance / cutoff
        closer = (ratio <= 1.0 + 1e-7) & (ratio < smallest)
        face[closer] = np.where(z >= 0, 1, 2)[closer]
        group[closer] = g
        smallest[closer] = ratio[closer]
    return face, group


def regions_of_system(simulations, L, coverage):
    """MSD and P of every region of one system, written to results/per_system/."""
    folder = run_folder(simulations, 'wet', L, coverage)
    u = mda.Universe(str(folder / 'confout.gro'), str(folder / 'traj_unwrap_whole.xtc'))
    water = u.select_atoms('resname TIP3 and name OH2')                       # the water oxygens
    silica = u.select_atoms('resname SUR')
    groups = [(u.select_atoms(f'resname OTS and name C{L}'), CARBON_CUTOFF),                  # CH3
              (u.select_atoms(f'resname OTS and name C* and not name C{L}'), CARBON_CUTOFF),  # CH2
              (u.select_atoms('(resname SUR or resname OTS) and name O*'), OXYGEN_CUTOFF)]     # silica and head-group O
    frames = u.trajectory[-LAST_FRAMES:]

    # for every region: sums over the origins of the squared displacements of the survivors, and the numbers of
    # survivors and of origins, for each block (that of the origin) and lag
    shape = (len(REGIONS), N_BLOCKS, len(LAGS))
    squared_sum, n_survivors, n_origins = np.zeros(shape), np.zeros(shape, dtype=int), np.zeros(shape, dtype=int)
    times = []                                     # time of every frame (ns)
    unwrapped_xy = np.zeros((len(water), 2))       # lateral positions without the jumps over the box edges
    history = deque(maxlen=max(LAGS) + 1)          # the lateral positions of this frame and of the 10 before
    codes = deque(maxlen=max(LAGS) + 1)            # the face of every water in every region, this frame and 10 before
    code = np.zeros((len(REGIONS), len(water)), dtype=int)
    run_length = np.zeros((len(REGIONS), len(water)), dtype=int)   # frames in the region without a break, up to now

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

        # the face of every water in every region (0: not in it); the frames in it without a break
        face, group = shell_faces_and_groups(xyz, groups, silica.positions[:, 2].mean(), box)
        code_before = code
        code = []
        for members in REGIONS.values():
            code.append(np.where(np.isin(group, members), face, 0))
        code = np.array(code)
        stayed = (code != 0) & (code == code_before)
        run_length = np.where(stayed, run_length + 1, np.where(code != 0, 1, 0))
        codes.append(code)

        # this frame is the end of the lags from the origins 1 ... 10 frames before
        for i, lag in enumerate(LAGS):
            if j < lag:
                continue
            block = min(N_BLOCKS - 1, (j - lag) * N_BLOCKS // len(frames))
            at_origin = codes[-1 - lag]
            for r in range(len(REGIONS)):
                n_origins[r, block, i] += np.count_nonzero(at_origin[r])
                survivors = (at_origin[r] != 0) & (run_length[r] >= lag + 1) & (code[r] == at_origin[r])
                squared = np.sum((unwrapped_xy[survivors] - history[-1 - lag][survivors]) ** 2, axis=1) * 0.01  # nm^2
                squared_sum[r, block, i] += squared.sum()
                n_survivors[r, block, i] += np.count_nonzero(survivors)

    # every region: per block, and the whole run with the SD of the MSD of the blocks
    dt = float(np.median(np.diff(times)))          # time between the frames (ns)
    block_rows, rows = [], []
    for r, region in enumerate(REGIONS):
        for i, lag in enumerate(LAGS):
            block_msd = np.divide(squared_sum[r, :, i], n_survivors[r, :, i], out=np.full(N_BLOCKS, np.nan),
                                  where=n_survivors[r, :, i] > 0)
            for b in range(N_BLOCKS):
                block_rows.append({'region': region, 'block': b + 1, 'lag_ps': lag * dt * 1000, 'MSD_nm2': block_msd[b],
                                   'n_survivors': n_survivors[r, b, i], 'n_origins': n_origins[r, b, i]})
            good = n_survivors[r, :, i] >= 20
            survivors, origins = n_survivors[r, :, i].sum(), n_origins[r, :, i].sum()
            if origins:
                P = survivors / origins
            else:
                P = np.nan
            if survivors:
                msd = squared_sum[r, :, i].sum() / survivors
            else:
                msd = np.nan
            if good.sum() >= 3:
                msd_error = np.std(block_msd[good], ddof=1)
            else:
                msd_error = np.nan
            rows.append({'region': region, 'lag_ps': lag * dt * 1000, 'n_origins': origins, 'n_survivors': survivors,
                         'P': P, 'MSD_nm2': msd, 'MSD_error_nm2': msd_error})
    pd.DataFrame(block_rows).to_csv(f'results/per_system/L{L}_G{coverage}_blocks.csv', index=False)
    pd.DataFrame(rows).to_csv(f'results/per_system/L{L}_G{coverage}.csv', index=False)


if __name__ == '__main__':
    simulations = simulations_folder()
    Path('results/per_system').mkdir(parents=True, exist_ok=True)
    systems = [(simulations, L, coverage) for L, coverage in find_systems(simulations, 'wet')]
    with Pool(PROCESSES) as pool:
        pool.starmap(regions_of_system, systems)
