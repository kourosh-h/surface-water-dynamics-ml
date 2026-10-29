"""
The results of step 1 of every system in two tables, with the reliability of every value.

A value of a region, system and lag is reliable when at least 200 molecule - time origin pairs survive until the lag
(both faces together), at least 3 of the ten 9-ns blocks have 20 or more of them, and P is at least 0.005.

Output:
    results/region_curves.csv: region, L, coverage (molecules/nm^2), lag_ps, MSD_nm2, MSD_error_nm2, P, P_error
        (standard deviation of the ten block values of P), reliable
    results/region_blocks.csv: region, L, coverage, block, lag_ps, MSD_nm2, n_survivors, n_origins
Run:  python3 step2_region_tables.py      (after step 1)
"""
import numpy as np
import pandas as pd

from systems import systems_in

MIN_SURVIVORS = 200             # molecule - origin pairs that survive until the lag, both faces together
MIN_BLOCK_SURVIVORS = 20        # survivors of a good block
MIN_GOOD_BLOCKS = 3             # good blocks of a reliable value
MIN_P = 0.005                   # smallest reliable P

curves, blocks = [], []
for L, coverage in systems_in('results/per_system'):
    system = {'L': L, 'coverage': float(coverage)}
    curves.append(pd.read_csv(f'results/per_system/L{L}_G{coverage}.csv').assign(**system))
    blocks.append(pd.read_csv(f'results/per_system/L{L}_G{coverage}_blocks.csv').assign(**system))
curves = pd.concat(curves)
blocks = pd.concat(blocks)

# reliability, from the blocks; and the standard deviation of P over the ten blocks
blocks['P'] = blocks.n_survivors / blocks.n_origins
blocks['good'] = (blocks.n_survivors >= MIN_BLOCK_SURVIVORS) & np.isfinite(blocks.MSD_nm2)
per_lag = blocks.groupby(['region', 'L', 'coverage', 'lag_ps']).agg(good_blocks=('good', 'sum'), P_error=('P', 'std'))
curves = curves.merge(per_lag.reset_index(), on=['region', 'L', 'coverage', 'lag_ps'])
curves['reliable'] = ((curves.n_survivors >= MIN_SURVIVORS) & (curves.good_blocks >= MIN_GOOD_BLOCKS)
                      & (curves.P >= MIN_P))

columns = ['region', 'L', 'coverage', 'lag_ps', 'MSD_nm2', 'MSD_error_nm2', 'P', 'P_error', 'reliable']
curves[columns].to_csv('results/region_curves.csv', index=False)
blocks[['region', 'L', 'coverage', 'block', 'lag_ps', 'MSD_nm2', 'n_survivors', 'n_origins']].to_csv(
    'results/region_blocks.csv', index=False)
