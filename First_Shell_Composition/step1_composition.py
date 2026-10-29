"""
The composition of the first hydration shell of every wet system: the share of the first-shell water at each group of
atoms (CH3, CH2 and the oxygens of the silica and the head groups), and the group with the most first-shell water.

Input: results/first_shell_water.csv of ../RDF_MSD_and_P (step 6; set by FIRST_SHELL_WATER): L, coverage
    (molecules/nm^2), group (CH3, CH2 or oxygen: the closest group of the first-shell water), water_per_nm2 (mean over
    the ten blocks), sd_blocks
Output:
    results/largest_group.csv: L, coverage, largest_group, percent (the group with the most first-shell water, and its
        share of the first-shell water)
    results/share_CH3.csv, results/share_CH2.csv, results/share_oxygen.csv: L, coverage, fraction (share of the
        first-shell water at the group)
Run:  python3 step1_composition.py      (after step 6 of ../RDF_MSD_and_P)
"""
from pathlib import Path

import pandas as pd

FIRST_SHELL_WATER = '../RDF_MSD_and_P/results/first_shell_water.csv'

Path('results').mkdir(exist_ok=True)
water = pd.read_csv(FIRST_SHELL_WATER)
table = water.pivot(index=['L', 'coverage'], columns='group', values='water_per_nm2')
shares = table.div(table.sum(axis=1), axis=0)            # fraction of the first-shell water at every group

# the largest group and its share, in percent
largest = shares.idxmax(axis=1)
pd.DataFrame({'largest_group': largest, 'percent': 100 * shares.max(axis=1)}).to_csv('results/largest_group.csv')

# share of every group
for group in ['CH3', 'CH2', 'oxygen']:
    shares[group].rename('fraction').to_csv(f'results/share_{group}.csv')
