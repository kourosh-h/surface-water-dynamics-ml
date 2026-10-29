"""
The two regimes of the first-shell water dynamics of the wet SAM systems: k-means grouping (scikit-learn) of the
shapes of their MSD and P curves.

The MSD and P of every system are divided by their values at 20 ps. The shape of the curves of a system: the
logarithms of these ratios at 40, 60, 80 and 100 ps, each scaled to mean 0 and SD 1 over the systems. k-means with two
groups (20 starts, random_state 123); regime 0 (faster renewal) is the group of the lower mean coverage.

Input:
    results/msd_and_P.csv of ../RDF_MSD_and_P (set by MSD_AND_P): one row per system and lag: L (chain length),
        coverage (molecules/nm^2), lag_ps (20, 40, 60, 80 and 100 ps), MSD_nm2 (survival-conditioned lateral MSD of
        the first-shell water, both faces of the slab together) and P (contact survival of the first-shell water)

Output:
    results/curves_MSD_rel.csv: L, coverage, regime (0 faster renewal, 1 slower renewal), lag_ps, MSD_rel
    results/curves_P_rel.csv: L, coverage, regime, lag_ps, P_rel
    results/regime_map.csv: the regime of every system, one row per chain length (column L), one column per coverage
    results/regime_mean_MSD_rel.csv and results/regime_mean_P_rel.csv: regime, lag_ps, mean, sd (standard deviation
        over the systems of the regime)

Run:  python3 step1_classify_water_dynamics.py      (after step 5 of ../RDF_MSD_and_P)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

MSD_AND_P = '../RDF_MSD_and_P/results/msd_and_P.csv'
SHAPE_LAGS = [40, 60, 80, 100]          # ps, the lags of the curve shapes


def curve_shapes(curves):
    """The shape of the curves of every system: the logarithms of MSD_rel and P_rel at the lags of SHAPE_LAGS."""
    table = curves.pivot(index=['L', 'coverage'], columns='lag_ps', values=['MSD_rel', 'P_rel'])
    shapes = pd.DataFrame(index=table.index)
    for lag in SHAPE_LAGS:
        shapes[f'log MSD_rel at {lag} ps'] = np.log(table['MSD_rel', lag])
        shapes[f'log P_rel at {lag} ps'] = np.log(table['P_rel', lag])
    return shapes


Path('results').mkdir(exist_ok=True)

# ----------------------------------------------------------------------
# 1. MSD and P of every system divided by their values at 20 ps
# ----------------------------------------------------------------------
curves = pd.read_csv(MSD_AND_P)
at_20_ps = curves.groupby(['L', 'coverage']).transform('first')     # the first row of every system is its 20-ps row
curves['MSD_rel'] = curves.MSD_nm2 / at_20_ps.MSD_nm2
curves['P_rel'] = curves.P / at_20_ps.P

# ----------------------------------------------------------------------
# 2. Two groups of systems (k-means) from the shapes of their curves, each shape value scaled to mean 0 and SD 1
# ----------------------------------------------------------------------
shapes = curve_shapes(curves)
X = StandardScaler().fit_transform(shapes)
groups = KMeans(n_clusters=2, random_state=123, n_init=20).fit(X).labels_
systems = shapes.index.to_frame(index=False)            # one row per system: L and coverage
# regime 0 (faster renewal) is the group of the lower mean coverage
if systems.coverage[groups == 0].mean() > systems.coverage[groups == 1].mean():
    groups = 1 - groups
systems['regime'] = groups

# the curves with the regime of their system
curves = curves.merge(systems[['L', 'coverage', 'regime']], on=['L', 'coverage'])
curves[['L', 'coverage', 'regime', 'lag_ps', 'MSD_rel']].to_csv('results/curves_MSD_rel.csv', index=False)
curves[['L', 'coverage', 'regime', 'lag_ps', 'P_rel']].to_csv('results/curves_P_rel.csv', index=False)

# map of the regimes: one row per chain length, one column per coverage
regime_map = systems.pivot(index='L', columns='coverage', values='regime')
regime_map.columns = [f'{coverage:.2f}' for coverage in regime_map.columns]
regime_map.to_csv('results/regime_map.csv')

# ----------------------------------------------------------------------
# 3. Mean and standard deviation of the curves of each regime at every lag, from the curve files
# ----------------------------------------------------------------------
for name in ['MSD_rel', 'P_rel']:
    curves = pd.read_csv(f'results/curves_{name}.csv')
    values = curves.groupby(['regime', 'lag_ps'])[name]
    table = pd.DataFrame({'mean': values.mean(), 'sd': values.std()})
    table.to_csv(f'results/regime_mean_{name}.csv')
