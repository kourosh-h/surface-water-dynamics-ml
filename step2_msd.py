"""
Step 2: the lateral mean squared displacement (MSD) of the first-shell water of every system, at the lag times of
LAGS.

First shell: a water molecule is in the first shell when its oxygen is within the contact distance of an atom of one
of the groups of shell_groups() in settings.py (contact distances from step 1). It belongs to the face of its closest
atom, of the group with the smallest distance / contact distance.

MSD(tau): over all time origins of the analyzed frames, the mean squared displacement in the plane of the surface
(x, y) of the water molecules that are in the first shell of a face at the origin and stay in it, without a break,
until the origin + tau; all analyzed faces together. Its uncertainty: the standard deviation of the MSD of the blocks
(a time origin belongs to the block it is in).

Output: results/msd/<system>.csv: lag_ps, MSD_nm2, MSD_error_nm2
Run:  python3 step2_msd.py      (after step 1)
"""
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import settings
from first_shell import faces_and_groups, open_system


def msd_of_system(system):
    """MSD of one system at every lag, written to results/msd/<system>.csv."""
    u, water, solid, groups = open_system(system)
    frames = u.trajectory[-settings.LAST_FRAMES:]
    n_lags = len(settings.LAGS)

    # sums over the time origins of the squared lateral displacements of the water molecules that stay, for each face,
    # block (that of the origin) and lag; and the number of these molecules, for each block and lag
    squared_sum = np.zeros((2, settings.N_BLOCKS, n_lags))
    n_survivors = np.zeros((settings.N_BLOCKS, n_lags), dtype=int)
    times = []                                     # time of every frame (ns)
    unwrapped_xy = np.zeros((len(water), 2))       # lateral positions without the jumps over the box edges
    history = []                                   # the lateral positions of this frame and of the frames before
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
        if len(history) > max(settings.LAGS) + 1:      # keep only the frames of the longest lag
            del history[0]

        # first shell; the frames in it without a break: one more if the water stays, 1 if it just came in
        face_before = face
        face, group = faces_and_groups(xyz, groups, solid.positions[:, 2].mean(), box)
        stayed = (face != 0) & (face == face_before)
        run_length = np.where(stayed, run_length + 1, np.where(face != 0, 1, 0))

        # this frame is the end of the lags from the origins 1, 2, ... frames before
        for i, lag in enumerate(settings.LAGS):
            if j < lag:
                continue
            block = min(settings.N_BLOCKS - 1, (j - lag) * settings.N_BLOCKS // len(frames))
            for side in [1, 2]:
                survivors = (run_length >= lag + 1) & (face == side)
                squared = np.sum((unwrapped_xy[survivors] - history[-1 - lag][survivors]) ** 2, axis=1) * 0.01  # nm^2
                squared_sum[side - 1, block, i] += squared.sum()
                n_survivors[block, i] += np.count_nonzero(survivors)

    # MSD, all faces together, and the standard deviation of the MSD of the blocks
    dt = float(np.median(np.diff(times)))          # time between the frames (ns)
    rows = []
    for i, lag in enumerate(settings.LAGS):
        block_sums = squared_sum[:, :, i].sum(axis=0)
        rows.append({'lag_ps': lag * dt * 1000, 'MSD_nm2': block_sums.sum() / n_survivors[:, i].sum(),
                     'MSD_error_nm2': np.std(block_sums / n_survivors[:, i], ddof=1)})
    pd.DataFrame(rows).to_csv(f'results/msd/{system["name"]}.csv', index=False)


if __name__ == '__main__':
    Path('results/msd').mkdir(parents=True, exist_ok=True)
    with Pool(settings.PROCESSES) as pool:
        pool.map(msd_of_system, settings.SYSTEMS)
