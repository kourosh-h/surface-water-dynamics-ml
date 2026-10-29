"""
First-shell water of every system at each group of atoms: the mean over the ten blocks and its standard deviation,
from results/first_shell_water_blocks.csv (step 5).

Output: results/first_shell_water.csv: L, coverage (molecules/nm^2), group, water_per_nm2, sd_blocks
Run:  python3 step6_first_shell_means.py      (after step 5)
"""
import pandas as pd

blocks = pd.read_csv('results/first_shell_water_blocks.csv')
rows = []
for (L, coverage, group), system in blocks.groupby(['L', 'coverage', 'group']):
    rows.append({'L': L, 'coverage': coverage, 'group': group, 'water_per_nm2': system.water_per_nm2.mean(),
                 'sd_blocks': system.water_per_nm2.std()})
pd.DataFrame(rows).sort_values(['group', 'L', 'coverage']).to_csv('results/first_shell_water.csv', index=False)
