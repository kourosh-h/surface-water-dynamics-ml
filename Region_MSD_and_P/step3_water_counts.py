"""
The numbers of water molecules in the analysis of the first-shell water of every wet system.

Water molecules in the box: from the table of the simulation systems (results/table_simulation_systems.csv of
sam-system-summary, set by SYSTEM_TABLE).
First shell: the molecule - origin pairs of every region at the origins of the 20-ps lag (step 1) divided by the number
of these origins (every saved frame of the last 90 ns but the last one): the mean number of first-shell water
molecules in a frame, both faces together, all and at each group (CH3, CH2, oxygen).
Molecule - origin pairs of all first-shell water: at the origins of the 20-ps lag, and those in contact until 20, 40
and 100 ps.

Output: results/first_shell_water_counts.csv: L, coverage (molecules/nm^2), water_in_box, first_shell, first_shell_CH3,
        first_shell_CH2, first_shell_oxygen, pairs_at_origins, pairs_in_contact_20ps, pairs_in_contact_40ps,
        pairs_in_contact_100ps
Run:  python3 step3_water_counts.py      (after step 1)
"""
import pandas as pd

from systems import systems_in

SYSTEM_TABLE = '../../sam-system-summary/results/table_simulation_systems.csv'
ORIGINS = 4500                  # time origins of the 20-ps lag: the 4501 frames of the last 90 ns but the last one

systems = pd.read_csv(SYSTEM_TABLE)
systems = systems[systems.state == 'wet']
rows = []
for L, coverage in systems_in('results/per_system'):
    water = systems[(systems.L == L) & (systems.coverage_molecules_per_nm2.round(2) == round(float(coverage), 2))]
    table = pd.read_csv(f'results/per_system/L{L}_G{coverage}.csv')
    table['lag'] = table.lag_ps.round().astype(int)
    at_20 = table[table.lag == 20].set_index('region')
    at_40 = table[table.lag == 40].set_index('region')
    at_100 = table[table.lag == 100].set_index('region')
    rows.append({'L': L, 'coverage': round(float(coverage), 2), 'water_in_box': int(water.water_molecules.iloc[0]),
                 'first_shell': at_20.n_origins['all'] / ORIGINS,
                 'first_shell_CH3': at_20.n_origins['CH3'] / ORIGINS,
                 'first_shell_CH2': at_20.n_origins['CH2'] / ORIGINS,
                 'first_shell_oxygen': at_20.n_origins['oxygen'] / ORIGINS,
                 'pairs_at_origins': at_20.n_origins['all'], 'pairs_in_contact_20ps': at_20.n_survivors['all'],
                 'pairs_in_contact_40ps': at_40.n_survivors['all'],
                 'pairs_in_contact_100ps': at_100.n_survivors['all']})
pd.DataFrame(rows).to_csv('results/first_shell_water_counts.csv', index=False)
