"""
Step 7: the two regimes of the water dynamics: two clusters of the systems from the shapes of their MSD and P curves
(k-means clustering with scikit-learn).

The MSD and P of every system are divided by their values at the first lag. The shape of the curves of a system: the
logarithms of these ratios at the lags of SHAPE_LAGS, each scaled to mean 0 and standard deviation 1 over the systems.
k-means with two clusters (20 starts). Regime 0 (faster renewal) is the cluster with the lower median of P at the last
lag of SHAPE_LAGS divided by P at the first lag.

Input:  results/msd_and_P.csv (step 6)
Output (folder results/clusters):
    regimes.csv: the labels of the system, regime (0 faster renewal, 1 slower renewal)
    curves_MSD_rel.csv: the labels, regime, lag_ps, MSD_rel
    curves_P_rel.csv: the labels, regime, lag_ps, P_rel
    regime_mean_MSD_rel.csv and regime_mean_P_rel.csv: regime, lag_ps, mean, sd (standard deviation over the systems
        of the regime)
Run:  python3 step7_clusters.py      (after step 6)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

import settings

Path('results/clusters').mkdir(parents=True, exist_ok=True)

# MSD and P of every system divided by their values at the first lag
curves = pd.read_csv('results/msd_and_P.csv')
first = curves.groupby(settings.LABELS).transform('first')     # the first row of every system is its first lag
curves['MSD_rel'] = curves.MSD_nm2 / first.MSD_nm2
curves['P_rel'] = curves.P / first.P

# the shape of the curves of every system: the logarithms of MSD_rel and P_rel at the lags of SHAPE_LAGS
table = curves.pivot(index=settings.LABELS, columns='lag_ps', values=['MSD_rel', 'P_rel'])
shapes = pd.DataFrame(index=table.index)
for lag in settings.SHAPE_LAGS:
    shapes[f'log MSD_rel at {lag} ps'] = np.log(table['MSD_rel', lag])
    shapes[f'log P_rel at {lag} ps'] = np.log(table['P_rel', lag])

# two clusters (k-means) of the shapes, each shape value scaled to mean 0 and standard deviation 1
X = StandardScaler().fit_transform(shapes)
regime = KMeans(n_clusters=2, random_state=settings.SEED, n_init=20).fit(X).labels_
# regime 0 (faster renewal) is the cluster whose P drops more until the last lag
last = shapes[f'log P_rel at {settings.SHAPE_LAGS[-1]} ps'].to_numpy()
if np.median(last[regime == 0]) > np.median(last[regime == 1]):
    regime = 1 - regime
systems = shapes.index.to_frame(index=False)            # one row per system: its labels
systems['regime'] = regime
systems.to_csv('results/clusters/regimes.csv', index=False)

# the curves with the regime of their system
curves = curves.merge(systems, on=settings.LABELS)
curves[settings.LABELS + ['regime', 'lag_ps', 'MSD_rel']].to_csv('results/clusters/curves_MSD_rel.csv', index=False)
curves[settings.LABELS + ['regime', 'lag_ps', 'P_rel']].to_csv('results/clusters/curves_P_rel.csv', index=False)

# mean and standard deviation of the curves of each regime at every lag, from the curve files
for name in ['MSD_rel', 'P_rel']:
    curves = pd.read_csv(f'results/clusters/curves_{name}.csv')
    values = curves.groupby(['regime', 'lag_ps'])[name]
    pd.DataFrame({'mean': values.mean(), 'sd': values.std()}).to_csv(f'results/clusters/regime_mean_{name}.csv')
