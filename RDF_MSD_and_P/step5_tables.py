"""
MSD, P and first-shell water of all systems in one table each, from the files of steps 2, 3 and 4.

Output:
    results/msd_and_P.csv: L, coverage (molecules/nm^2), lag_ps, MSD_nm2, P at the lags 20 ... 100 ps: the input of
        the classification in ../ML_Classification
    results/msd_and_P_all_lags.csv: L, coverage, lag_ps, MSD_nm2, MSD_error_nm2, P at all lags (20 ... 200 ps)
    results/first_shell_water_blocks.csv: L, coverage, block, group, water_per_nm2

Run:  python3 step5_tables.py      (after steps 2, 3 and 4)
"""
import pandas as pd

from systems import systems_in

msd_and_P, first_shell_water = [], []
for L, coverage in systems_in('results/MSD'):
    system = pd.read_csv(f'results/MSD/L{L}_G{coverage}.csv').merge(pd.read_csv(f'results/P/L{L}_G{coverage}.csv'),
                                                                    on='lag_ps')
    system.insert(0, 'L', L)
    system.insert(1, 'coverage', float(coverage))
    msd_and_P.append(system)
    water = pd.read_csv(f'results/first_shell_water/L{L}_G{coverage}.csv')
    water.insert(0, 'L', L)
    water.insert(1, 'coverage', float(coverage))
    first_shell_water.append(water)
msd_and_P = pd.concat(msd_and_P)
msd_and_P.to_csv('results/msd_and_P_all_lags.csv', index=False)
pd.concat(first_shell_water).to_csv('results/first_shell_water_blocks.csv', index=False)

# the input of the classification: the lags 20 ... 100 ps, in whole picoseconds
ml_input = msd_and_P[msd_and_P.lag_ps < 101].copy()
ml_input['lag_ps'] = ml_input.lag_ps.round().astype(int)
ml_input[['L', 'coverage', 'lag_ps', 'MSD_nm2', 'P']].to_csv('results/msd_and_P.csv', index=False)
