"""
Step 8: tests of the regimes, on the MSD and P of all first-shell water of step 5 (region 'all').

A system enters the tests when its values are reliable at the first lag and at every lag of SHAPE_LAGS. The regimes
are found as in step 7: the logarithms of MSD and P divided by their values at the first lag, at the lags of
SHAPE_LAGS, standardized; k-means with N_CLUSTERS clusters (20 starts), numbered from the fastest to the slowest
renewal (by the median of P at the last lag of SHAPE_LAGS divided by P at the first lag, from the lowest).

Tests:
    the silhouette score and the bootstrap stability (median adjusted Rand index between the clusters and those of
        k-means refitted to RESAMPLES resamples of the systems) for every number of clusters of CLUSTER_NUMBERS: the
        number of regimes N_CLUSTERS is chosen from them;
    the regimes from each half of the blocks (their molecule - time origin pairs pooled);
    the regimes from the curves up to the longest lag at which all systems are reliable;
    the renewal index P(last lag of SHAPE_LAGS) / P(first lag), and the exponent alpha of MSD ~ tau^alpha from the
        first lag to the last lag of SHAPE_LAGS (straight-line fit of log MSD against log tau).

Input:  results/region_curves.csv and results/region_blocks.csv (step 6)
Output (folder results/cluster_tests):
    k_scores.csv: k, silhouette, bootstrap_stability
    renewal_index.csv: the labels of the system, regime (0 for the fastest renewal), P_ratio
    msd_exponent.csv: the labels, regime, alpha
    cluster_tests.csv: the number of systems, and the systems that keep their regime (and the adjusted Rand index) in
        the two halves of the run and with the curves up to the longest reliable lag
Run:  python3 step8_cluster_tests.py      (after step 6)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import StandardScaler

import settings

REGION = 'all'                          # all first-shell water


def curve_shapes(table, lags, first_lag):
    """The logarithms of MSD and P divided by their values at the first lag, at the given lags: one row per system."""
    shapes = pd.DataFrame(index=table.index)
    for lag in lags:
        shapes[f'log MSD_rel at {lag} ps'] = np.log(table['MSD_nm2', lag] / table['MSD_nm2', first_lag])
        shapes[f'log P_rel at {lag} ps'] = np.log(table['P', lag] / table['P', first_lag])
    return shapes


def regimes(shapes):
    """
    k-means with N_CLUSTERS clusters on the standardized shapes, the clusters numbered from the fastest to the slowest
    renewal: by the median of log(P at the last lag of SHAPE_LAGS / P at the first lag), from the lowest.
    """
    X = StandardScaler().fit_transform(shapes)
    clusters = KMeans(n_clusters=settings.N_CLUSTERS, random_state=settings.SEED, n_init=20).fit(X).labels_
    last = shapes[f'log P_rel at {settings.SHAPE_LAGS[-1]} ps'].to_numpy()
    medians = []
    for cluster in range(settings.N_CLUSTERS):
        medians.append(np.median(last[clusters == cluster]))
    order = np.argsort(medians, kind='stable')          # the clusters from the lowest to the highest median
    regime = np.zeros(len(clusters), dtype=int)
    for number, cluster in enumerate(order):
        regime[clusters == cluster] = number
    return regime


Path('results/cluster_tests').mkdir(parents=True, exist_ok=True)
curves = pd.read_csv('results/region_curves.csv')
curves['lag_ps'] = curves.lag_ps.round().astype(int)
blocks = pd.read_csv('results/region_blocks.csv')
blocks['lag_ps'] = blocks.lag_ps.round().astype(int)
rows = curves[curves.region == REGION]
table = rows.pivot(index=settings.LABELS, columns='lag_ps', values=['MSD_nm2', 'P'])
reliable = rows.pivot(index=settings.LABELS, columns='lag_ps', values='reliable')
lags = sorted(rows.lag_ps.unique())                     # all lags (ps)
first_lag = lags[0]
last_lag = settings.SHAPE_LAGS[-1]
systems = reliable.index[reliable[[first_lag] + settings.SHAPE_LAGS].all(axis=1)]    # the systems of the tests
table = table.loc[systems]
shapes = curve_shapes(table, settings.SHAPE_LAGS, first_lag)
X = StandardScaler().fit_transform(shapes)
regime = regimes(shapes)
result = systems.to_frame(index=False)
result['regime'] = regime

# silhouette score and bootstrap stability for every number of clusters
scores = []
for k in settings.CLUSTER_NUMBERS:
    clusters = KMeans(n_clusters=k, random_state=settings.SEED, n_init=20).fit(X).labels_
    random_numbers = np.random.default_rng(settings.SEED)
    agreement = []
    for repeat in range(settings.RESAMPLES):
        sample = random_numbers.choice(len(X), len(X), replace=True)
        refitted = KMeans(n_clusters=k, random_state=settings.SEED, n_init=10).fit(X[sample])
        agreement.append(adjusted_rand_score(clusters, refitted.predict(X)))
    scores.append({'k': k, 'silhouette': silhouette_score(X, clusters), 'bootstrap_stability': np.median(agreement)})
pd.DataFrame(scores).to_csv('results/cluster_tests/k_scores.csv', index=False)

# renewal index and exponent alpha of the MSD
renewal = result.copy()
renewal['P_ratio'] = (table['P', last_lag] / table['P', first_lag]).to_numpy()
renewal.to_csv('results/cluster_tests/renewal_index.csv', index=False)
fit_lags = np.array([first_lag] + settings.SHAPE_LAGS)
alpha = []
for system in systems:
    msd = table['MSD_nm2'].loc[system, fit_lags].to_numpy(float)
    alpha.append(np.polyfit(np.log(fit_lags), np.log(msd), 1)[0])
exponent = result.copy()
exponent['alpha'] = alpha
exponent.to_csv('results/cluster_tests/msd_exponent.csv', index=False)

# tests: the two halves of the blocks, and the curves up to the longest lag at which all systems are reliable
test = {'region': REGION, 'systems': len(systems)}
half = settings.N_BLOCKS // 2
for name, half_blocks in [('first_half', range(1, half + 1)), ('second_half', range(half + 1, settings.N_BLOCKS + 1))]:
    part = blocks[(blocks.region == REGION) & blocks.block.isin(half_blocks)].copy()
    part['squared_sum'] = part.MSD_nm2 * part.n_survivors
    sums = part.groupby(settings.LABELS + ['lag_ps'])[['squared_sum', 'n_survivors', 'n_origins']].sum()
    half_table = pd.DataFrame({'MSD_nm2': sums.squared_sum / sums.n_survivors,
                               'P': sums.n_survivors / sums.n_origins}).unstack('lag_ps').loc[systems]
    half_regime = regimes(curve_shapes(half_table, settings.SHAPE_LAGS, first_lag))
    test[f'{name}_same'] = int((half_regime == regime).sum())
    test[f'{name}_ARI'] = adjusted_rand_score(regime, half_regime)
longest = last_lag
for lag in lags:
    if lag <= last_lag:
        continue
    if not reliable.loc[systems, lag].all():
        break
    longest = lag
test['longest_lag_ps'] = longest
if longest > last_lag:
    longer_lags = []
    for lag in lags:
        if first_lag < lag <= longest:
            longer_lags.append(lag)
    longer_regime = regimes(curve_shapes(table, longer_lags, first_lag))
    test['longer_same'] = int((longer_regime == regime).sum())
    test['longer_ARI'] = adjusted_rand_score(regime, longer_regime)
pd.DataFrame([test]).to_csv('results/cluster_tests/cluster_tests.csv', index=False)
