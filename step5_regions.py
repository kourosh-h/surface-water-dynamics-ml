"""
Step 5: the lateral MSD and the contact survival P of the water of every region of the first shell (REGIONS in
settings.py), for every system, at the lag times of LAGS, per block.

First shell: as in step 2. A water molecule of the shell is in a region when its group belongs to the region, and it
stays in the region while it is in the first shell, on the same face, and its group belongs to the region.

For every time origin of the analyzed frames and every lag: the molecules in the region of a face at the origin
(origins), those that stay in it, without a break, until the origin + tau (survivors), and the squared displacement
of the survivors in the plane of the surface (x, y); all analyzed faces together, summed over the origins of each
block (a time origin belongs to the block it is in).

Output: results/regions/<system>_blocks.csv: region, block, lag_ps, MSD_nm2 (of the block), n_survivors, n_origins
        results/regions/<system>.csv: region, lag_ps, n_origins, n_survivors, P, MSD_nm2, MSD_error_nm2 (standard
            deviation of the MSD of the blocks with MIN_BLOCK_SURVIVORS or more survivors, when there are
            MIN_GOOD_BLOCKS or more of them)
Run:  python3 step5_regions.py      (after step 1)
"""
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import settings
from first_shell import faces_and_groups, open_system


def regions_of_system(system):
    """MSD and P of every region of one system, written to results/regions/."""
    u, water, solid, groups = open_system(system)
    frames = u.trajectory[-settings.LAST_FRAMES:]
    n_lags = len(settings.LAGS)

    # the groups of every region, as their numbers in the list of groups
    group_names = []
    for name, atoms, cutoff in groups:
        group_names.append(name)
    regions = []
    for region, members in settings.REGIONS.items():
        numbers = []
        for member in members:
            numbers.append(group_names.index(member))
        regions.append(numbers)

    # for every region: sums over the origins of the squared displacements of the survivors, and the numbers of
    # survivors and of origins, for each block (that of the origin) and lag
    shape = (len(regions), settings.N_BLOCKS, n_lags)
    squared_sum, n_survivors, n_origins = np.zeros(shape), np.zeros(shape, dtype=int), np.zeros(shape, dtype=int)
    times = []                                     # time of every frame (ns)
    unwrapped_xy = np.zeros((len(water), 2))       # lateral positions without the jumps over the box edges
    history = []                                   # the lateral positions of this frame and of the frames before
    codes = []                                     # the face of every water in every region, this frame and before
    code = np.zeros((len(regions), len(water)), dtype=int)
    run_length = np.zeros((len(regions), len(water)), dtype=int)   # frames in the region without a break, up to now

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
        face, group = faces_and_groups(xyz, groups, solid.positions[:, 2].mean(), box)
        code_before = code
        code = []
        for numbers in regions:
            code.append(np.where(np.isin(group, numbers), face, 0))
        code = np.array(code)
        stayed = (code != 0) & (code == code_before)
        run_length = np.where(stayed, run_length + 1, np.where(code != 0, 1, 0))
        codes.append(code)
        if len(history) > max(settings.LAGS) + 1:      # keep only the frames of the longest lag
            del history[0]
            del codes[0]

        # this frame is the end of the lags from the origins 1, 2, ... frames before
        for i, lag in enumerate(settings.LAGS):
            if j < lag:
                continue
            block = min(settings.N_BLOCKS - 1, (j - lag) * settings.N_BLOCKS // len(frames))
            at_origin = codes[-1 - lag]
            for r in range(len(regions)):
                n_origins[r, block, i] += np.count_nonzero(at_origin[r])
                survivors = (at_origin[r] != 0) & (run_length[r] >= lag + 1) & (code[r] == at_origin[r])
                squared = np.sum((unwrapped_xy[survivors] - history[-1 - lag][survivors]) ** 2, axis=1) * 0.01  # nm^2
                squared_sum[r, block, i] += squared.sum()
                n_survivors[r, block, i] += np.count_nonzero(survivors)

    # every region: per block, and the whole run with the standard deviation of the MSD of the blocks
    dt = float(np.median(np.diff(times)))          # time between the frames (ns)
    block_rows = []
    rows = []
    for r, region in enumerate(settings.REGIONS):
        for i, lag in enumerate(settings.LAGS):
            block_msd = np.full(settings.N_BLOCKS, np.nan)
            for b in range(settings.N_BLOCKS):
                if n_survivors[r, b, i] > 0:
                    block_msd[b] = squared_sum[r, b, i] / n_survivors[r, b, i]
                block_rows.append({'region': region, 'block': b + 1, 'lag_ps': lag * dt * 1000, 'MSD_nm2': block_msd[b],
                                   'n_survivors': n_survivors[r, b, i], 'n_origins': n_origins[r, b, i]})
            good = n_survivors[r, :, i] >= settings.MIN_BLOCK_SURVIVORS
            survivors, origins = n_survivors[r, :, i].sum(), n_origins[r, :, i].sum()
            if origins:
                P = survivors / origins
            else:
                P = np.nan
            if survivors:
                msd = squared_sum[r, :, i].sum() / survivors
            else:
                msd = np.nan
            if good.sum() >= settings.MIN_GOOD_BLOCKS:
                msd_error = np.std(block_msd[good], ddof=1)
            else:
                msd_error = np.nan
            rows.append({'region': region, 'lag_ps': lag * dt * 1000, 'n_origins': origins, 'n_survivors': survivors,
                         'P': P, 'MSD_nm2': msd, 'MSD_error_nm2': msd_error})
    pd.DataFrame(block_rows).to_csv(f'results/regions/{system["name"]}_blocks.csv', index=False)
    pd.DataFrame(rows).to_csv(f'results/regions/{system["name"]}.csv', index=False)


if __name__ == '__main__':
    Path('results/regions').mkdir(parents=True, exist_ok=True)
    with Pool(settings.PROCESSES) as pool:
        pool.map(regions_of_system, settings.SYSTEMS)
