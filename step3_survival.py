"""
Step 3: the contact survival probability P(tau) of the first-shell water of every system, at the lag times of LAGS.

First shell: as in step 2. P(tau): over all time origins of the analyzed frames, the fraction of the water molecules
in the first shell of a face at the origin that stay in it, without a break, until the origin + tau; all analyzed
faces together.

Output: results/survival/<system>.csv: lag_ps, P
Run:  python3 step3_survival.py      (after step 1)
"""
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import settings
from first_shell import faces_and_groups, open_system


def survival_of_system(system):
    """P of one system at every lag, written to results/survival/<system>.csv."""
    u, water, solid, groups = open_system(system)
    times = []                                     # time of every frame (ns)
    in_shell = []                                  # number of water molecules in the first shell, in every frame
    n_survivors = [0] * len(settings.LAGS)         # summed over the frames: water in the shell since lag frames before
    face = np.zeros(len(water), dtype=int)
    run_length = np.zeros(len(water), dtype=int)   # frames in the shell of the same face without a break, up to now
    for ts in u.trajectory[-settings.LAST_FRAMES:]:
        times.append(ts.time / 1000.0)
        box = ts.dimensions.astype(np.float64)
        face_before = face
        face, group = faces_and_groups(water.positions, groups, solid.positions[:, 2].mean(), box)
        stayed = (face != 0) & (face == face_before)
        run_length = np.where(stayed, run_length + 1, np.where(face != 0, 1, 0))
        in_shell.append(np.count_nonzero(face))
        for i, lag in enumerate(settings.LAGS):
            n_survivors[i] += np.count_nonzero(run_length >= lag + 1)

    # P: the survivors divided by the water molecules in the shell at the origins (all frames but the last lag ones)
    dt = float(np.median(np.diff(times)))          # time between the frames (ns)
    rows = []
    for i, lag in enumerate(settings.LAGS):
        rows.append({'lag_ps': lag * dt * 1000, 'P': n_survivors[i] / sum(in_shell[:len(in_shell) - lag])})
    pd.DataFrame(rows).to_csv(f'results/survival/{system["name"]}.csv', index=False)


if __name__ == '__main__':
    Path('results/survival').mkdir(parents=True, exist_ok=True)
    with Pool(settings.PROCESSES) as pool:
        pool.map(survival_of_system, settings.SYSTEMS)
