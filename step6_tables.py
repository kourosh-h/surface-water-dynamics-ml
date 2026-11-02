"""
Step 6: the results of steps 1 to 5 of all systems in tables. Every part below reads the files it needs from results/.

Output (folder results):
    msd_and_P_all_lags.csv: the labels of the system, lag_ps, MSD_nm2, MSD_error_nm2, P (all first-shell water, all
        lags)
    msd_and_P.csv: the labels, lag_ps (whole ps), MSD_nm2, P at the lags up to the last of SHAPE_LAGS: the input of
        step 7
    first_shell_water_blocks.csv: the labels, block, group, water_per_nm2
    first_shell_water.csv: the labels, group, water_per_nm2 and sd_blocks (mean and standard deviation over the
        blocks)
    first_shell_water_mean_over_<MEAN_OVER>.csv: group, the other labels, water_per_nm2 and water_per_nm2_sd (mean and
        standard deviation over the systems that differ only in MEAN_OVER)
    largest_group.csv: the labels, largest_group, percent (the group with the most first-shell water, and its share of
        the first-shell water)
    share_<group>.csv: the labels, fraction (the share of the first-shell water at the group)
    region_curves.csv: region, the labels, lag_ps, MSD_nm2, MSD_error_nm2, P, P_error (standard deviation of the block
        values of P), reliable
    region_blocks.csv: region, the labels, block, lag_ps, MSD_nm2, n_survivors, n_origins
    first_shell_water_counts.csv: the labels, water_in_box, first_shell (mean number of first-shell water molecules in
        a frame, all faces together), first_shell_<group>, pairs_at_origins (molecule - time origin pairs at the first
        lag), pairs_in_contact_<lag>ps (those in contact until the lags of COUNT_LAGS)
Run:  python3 step6_tables.py      (after steps 1 to 5)
"""
import numpy as np
import pandas as pd

import settings

group_names = []                        # the groups of the first shell
for name in settings.shell_groups(settings.SYSTEMS[0]):
    group_names.append(name)

# ----------------------------------------------------------------------
# MSD and P of all first-shell water (steps 2 and 3), and the first-shell water at each group per block (step 4)
# ----------------------------------------------------------------------
msd_and_P = []
first_shell_water = []
for system in settings.SYSTEMS:
    name = system['name']
    table = pd.read_csv(f'results/msd/{name}.csv').merge(pd.read_csv(f'results/survival/{name}.csv'), on='lag_ps')
    water = pd.read_csv(f'results/shell_water/{name}.csv')
    for k, label in enumerate(settings.LABELS):
        table.insert(k, label, system[label])
        water.insert(k, label, system[label])
    msd_and_P.append(table)
    first_shell_water.append(water)
msd_and_P = pd.concat(msd_and_P)
msd_and_P.to_csv('results/msd_and_P_all_lags.csv', index=False)
pd.concat(first_shell_water).to_csv('results/first_shell_water_blocks.csv', index=False)

# the input of the clustering: the lags up to the last lag of the curve shapes, in whole picoseconds
short = msd_and_P[msd_and_P.lag_ps < max(settings.SHAPE_LAGS) + 1].copy()
short['lag_ps'] = short.lag_ps.round().astype(int)
short[settings.LABELS + ['lag_ps', 'MSD_nm2', 'P']].to_csv('results/msd_and_P.csv', index=False)

# ----------------------------------------------------------------------
# First-shell water of every system at each group: mean and standard deviation over the blocks
# ----------------------------------------------------------------------
blocks = pd.read_csv('results/first_shell_water_blocks.csv')
rows = []
for keys, part in blocks.groupby(settings.LABELS + ['group']):
    row = dict(zip(settings.LABELS + ['group'], keys))
    row['water_per_nm2'] = part.water_per_nm2.mean()
    row['sd_blocks'] = part.water_per_nm2.std()
    rows.append(row)
pd.DataFrame(rows).sort_values(['group'] + settings.LABELS).to_csv('results/first_shell_water.csv', index=False)

# mean and standard deviation over the systems that differ only in MEAN_OVER
if settings.MEAN_OVER is not None:
    water = pd.read_csv('results/first_shell_water.csv')
    others = []
    for label in settings.LABELS:
        if label != settings.MEAN_OVER:
            others.append(label)
    values = water.groupby(['group'] + others).water_per_nm2
    pd.DataFrame({'water_per_nm2': values.mean(), 'water_per_nm2_sd': values.std()}).to_csv(
        f'results/first_shell_water_mean_over_{settings.MEAN_OVER}.csv')

# ----------------------------------------------------------------------
# Composition of the first shell: the share of the first-shell water at every group, and the largest group
# ----------------------------------------------------------------------
water = pd.read_csv('results/first_shell_water.csv')
table = water.pivot(index=settings.LABELS, columns='group', values='water_per_nm2')
shares = table.div(table.sum(axis=1), axis=0)           # the fraction of the first-shell water at every group
largest = shares.idxmax(axis=1)
pd.DataFrame({'largest_group': largest, 'percent': 100 * shares.max(axis=1)}).to_csv('results/largest_group.csv')
for group in group_names:
    shares[group].rename('fraction').to_csv(f'results/share_{group}.csv')

# ----------------------------------------------------------------------
# MSD and P of every region (step 5), with the reliability of every value
# ----------------------------------------------------------------------
curves = []
blocks = []
for system in settings.SYSTEMS:
    curve = pd.read_csv(f'results/regions/{system["name"]}.csv')
    block = pd.read_csv(f'results/regions/{system["name"]}_blocks.csv')
    for label in settings.LABELS:
        curve[label] = system[label]
        block[label] = system[label]
    curves.append(curve)
    blocks.append(block)
curves = pd.concat(curves)
blocks = pd.concat(blocks)

# reliability, from the blocks; and the standard deviation of P over the blocks
keys = ['region'] + settings.LABELS + ['lag_ps']
blocks['P'] = blocks.n_survivors / blocks.n_origins
blocks['good'] = (blocks.n_survivors >= settings.MIN_BLOCK_SURVIVORS) & np.isfinite(blocks.MSD_nm2)
per_lag = blocks.groupby(keys).agg(good_blocks=('good', 'sum'), P_error=('P', 'std'))
curves = curves.merge(per_lag.reset_index(), on=keys)
curves['reliable'] = ((curves.n_survivors >= settings.MIN_SURVIVORS) & (curves.good_blocks >= settings.MIN_GOOD_BLOCKS)
                      & (curves.P >= settings.MIN_P))
columns = ['region'] + settings.LABELS + ['lag_ps', 'MSD_nm2', 'MSD_error_nm2', 'P', 'P_error', 'reliable']
curves[columns].to_csv('results/region_curves.csv', index=False)
columns = ['region'] + settings.LABELS + ['block', 'lag_ps', 'MSD_nm2', 'n_survivors', 'n_origins']
blocks[columns].to_csv('results/region_blocks.csv', index=False)

# ----------------------------------------------------------------------
# The numbers of water molecules: in the box (step 1), in the first shell, and of the molecule - time origin pairs
# ----------------------------------------------------------------------
water_in_box = pd.read_csv('results/water_molecules.csv', index_col='system').water_molecules
origins = settings.LAST_FRAMES - settings.LAGS[0]       # the time origins of the first lag
rows = []
for system in settings.SYSTEMS:
    table = pd.read_csv(f'results/regions/{system["name"]}.csv')
    table['lag'] = table.lag_ps.round().astype(int)
    first = table[table.lag == table.lag.min()].set_index('region')
    row = {}
    for label in settings.LABELS:
        row[label] = system[label]
    row['water_in_box'] = int(water_in_box[system['name']])
    # the molecule - origin pairs at the first lag divided by the number of origins: first-shell water in a frame
    row['first_shell'] = first.n_origins['all'] / origins
    for group in group_names:
        if group in settings.REGIONS:
            row[f'first_shell_{group}'] = first.n_origins[group] / origins
    row['pairs_at_origins'] = first.n_origins['all']
    for lag in settings.COUNT_LAGS:
        at_lag = table[table.lag == lag].set_index('region')
        row[f'pairs_in_contact_{lag}ps'] = at_lag.n_survivors['all']
    rows.append(row)
pd.DataFrame(rows).to_csv('results/first_shell_water_counts.csv', index=False)
