"""
Two tables of the first-shell water, from the files of steps 5 and 6 and of step 1 of ../Region_MSD_and_P.

Output:
    results/first_shell_water_mean_over_L.csv: group, coverage (molecules/nm^2), water_per_nm2, water_per_nm2_sd: the
        first-shell water per nm^2 at the CH3 and CH2 groups and at the oxygens, mean and SD over the chain lengths
    results/msd_and_P_20_40ps.csv: L, coverage_molecules_per_nm2, MSD_20ps_nm2, MSD_20ps_sd_nm2, MSD_40ps_nm2,
        MSD_40ps_sd_nm2, P_20ps, P_20ps_sd, P_40ps, P_40ps_sd: MSD and P with their SDs at 20 and 40 ps of the
        systems at a few coverages (TABLE_COVERAGES); the SD of P is that of the ten block values of P, the survivors
        divided by the origins of the block (all first-shell water, results/per_system of ../Region_MSD_and_P)

Run:  python3 step7_tables.py      (after steps 5 and 6, and step 1 of ../Region_MSD_and_P)
"""
import pandas as pd

from systems import systems_in

TABLE_COVERAGES = ['0.77', '2.40', '3.90']                  # molecules/nm^2, the coverages of the MSD and P table
REGION_RESULTS = '../Region_MSD_and_P/results/per_system'   # block counts of every system (all first-shell water)

# first-shell water at each group: mean and SD over the chain lengths
water = pd.read_csv('results/first_shell_water.csv')
water = water[water.group.isin(['CH3', 'CH2', 'oxygen'])]
over_L = water.groupby(['group', 'coverage']).water_per_nm2
pd.DataFrame({'water_per_nm2': over_L.mean(), 'water_per_nm2_sd': over_L.std()}).to_csv(
    'results/first_shell_water_mean_over_L.csv')

# the coverage of every system as written in the file names
file_coverage = {}
for L, coverage in systems_in('results/MSD'):
    file_coverage[L, round(float(coverage), 2)] = coverage

# MSD and P with their SDs at 20 and 40 ps
msd_and_P = pd.read_csv('results/msd_and_P_all_lags.csv')
rows = []
for system in msd_and_P[msd_and_P.lag_ps.round() == 40].itertuples():
    if f'{system.coverage:.2f}' not in TABLE_COVERAGES:
        continue
    at_20 = msd_and_P[(msd_and_P.L == system.L) & (msd_and_P.coverage == system.coverage)
                      & (msd_and_P.lag_ps.round() == 20)].iloc[0]
    name = file_coverage[system.L, round(system.coverage, 2)]
    blocks = pd.read_csv(f'{REGION_RESULTS}/L{system.L}_G{name}_blocks.csv')
    blocks = blocks[blocks.region == 'all']
    P_sd = {}                           # SD of the block values of P at every lag (ps)
    for lag, part in blocks.groupby('lag_ps'):
        P_sd[round(lag)] = (part.n_survivors / part.n_origins).std()
    rows.append({'L': system.L, 'coverage_molecules_per_nm2': f'{system.coverage:.2f}',
                 'MSD_20ps_nm2': f'{at_20.MSD_nm2:.4f}', 'MSD_20ps_sd_nm2': f'{at_20.MSD_error_nm2:.4f}',
                 'MSD_40ps_nm2': f'{system.MSD_nm2:.3f}', 'MSD_40ps_sd_nm2': f'{system.MSD_error_nm2:.3f}',
                 'P_20ps': f'{at_20.P:.3f}', 'P_20ps_sd': f'{P_sd[20]:.3f}',
                 'P_40ps': f'{system.P:.3f}', 'P_40ps_sd': f'{P_sd[40]:.3f}'})
pd.DataFrame(rows).to_csv('results/msd_and_P_20_40ps.csv', index=False)
