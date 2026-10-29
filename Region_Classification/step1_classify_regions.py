"""
The two regimes of all first-shell water (k-means) and their tests, from the MSD and P of every region of the first
shell (../Region_MSD_and_P).

A system enters the classification when its values are reliable at every lag from 20 to 100 ps
(results/region_curves.csv of ../Region_MSD_and_P).

Classification: the eight values of a system are the logarithms of MSD and P divided by their values at 20 ps, at 40,
60, 80 and 100 ps; standardized to mean 0 and SD 1; k-means with k = 2 (scikit-learn, 20 starts, random_state 123).
Regime 0 (faster renewal) is the group with the lower median P(100 ps)/P(20 ps).
Tests: silhouette score and bootstrap stability for k = 2, 3, 4 (200 resamples of the systems, seed 123); the
classification of each half of the run (blocks 1-5 and 6-10, their molecule - origin pairs pooled); the classification
with the curves up to the longest lag at which all systems are reliable; the renewal index P(100 ps)/P(20 ps) and the
exponent alpha of MSD ~ tau^alpha from 20 to 100 ps (straight-line fit of log MSD against log tau).

Input: results/region_curves.csv and results/region_blocks.csv of ../Region_MSD_and_P (step 2), region 'all'
Output: results/k_scores.csv: k, silhouette, bootstrap_stability
        results/renewal_index.csv: L, coverage (molecules/nm^2), regime (0 faster renewal, 1 slower renewal), P_ratio
        results/msd_exponent.csv: L, coverage, regime, alpha
        results/region_tests.csv: the number of systems, and the systems that keep their regime (and the adjusted Rand
            index) in the two halves of the run and with the curves up to the longest reliable lag
Run:  python3 step1_classify_regions.py      (after step 2 of ../Region_MSD_and_P)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import StandardScaler

REGION_RESULTS = '../Region_MSD_and_P/results'
REGION = 'all'                          # all first-shell water
LAGS = [40, 60, 80, 100]                # ps, the lags of the curve shapes
GROUP_NUMBERS = [2, 3, 4]               # the numbers of groups k of the scores
RESAMPLES = 200                         # bootstrap resamples of the systems


def curve_shapes(table, lags):
    """The logarithms of MSD and P divided by their values at 20 ps, at the given lags: one row per system."""
    shapes = pd.DataFrame(index=table.index)
    for lag in lags:
        shapes[f'log MSD_rel at {lag} ps'] = np.log(table['MSD_nm2', lag] / table['MSD_nm2', 20])
        shapes[f'log P_rel at {lag} ps'] = np.log(table['P', lag] / table['P', 20])
    return shapes


def regimes(shapes):
    """k-means with k = 2 on the standardized shapes; regime 0 is the group with the lower median P(100)/P(20)."""
    labels = KMeans(n_clusters=2, random_state=123, n_init=20).fit(StandardScaler().fit_transform(shapes)).labels_
    p_rel = shapes['log P_rel at 100 ps'].to_numpy()
    if np.median(p_rel[labels == 0]) > np.median(p_rel[labels == 1]):
        labels = 1 - labels
    return labels


Path('results').mkdir(exist_ok=True)
curves = pd.read_csv(f'{REGION_RESULTS}/region_curves.csv')
curves['lag_ps'] = curves.lag_ps.round().astype(int)
blocks = pd.read_csv(f'{REGION_RESULTS}/region_blocks.csv')
blocks['lag_ps'] = blocks.lag_ps.round().astype(int)
rows = curves[curves.region == REGION]
table = rows.pivot(index=['L', 'coverage'], columns='lag_ps', values=['MSD_nm2', 'P'])
reliable = rows.pivot(index=['L', 'coverage'], columns='lag_ps', values='reliable')
systems = reliable.index[reliable[[20] + LAGS].all(axis=1)]          # the systems of the classification
table = table.loc[systems]
shapes = curve_shapes(table, LAGS)
X = StandardScaler().fit_transform(shapes)
regime = regimes(shapes)
result = systems.to_frame(index=False).assign(regime=regime)

# silhouette score and bootstrap stability for every number of groups
scores = []
for k in GROUP_NUMBERS:
    groups = KMeans(n_clusters=k, random_state=123, n_init=20).fit(X).labels_
    random_numbers = np.random.default_rng(123)
    agreement = []
    for repeat in range(RESAMPLES):
        sample = random_numbers.choice(len(X), len(X), replace=True)
        resample_groups = KMeans(n_clusters=k, random_state=123, n_init=10).fit(X[sample])
        agreement.append(adjusted_rand_score(groups, resample_groups.predict(X)))
    scores.append({'k': k, 'silhouette': silhouette_score(X, groups), 'bootstrap_stability': np.median(agreement)})
pd.DataFrame(scores).to_csv('results/k_scores.csv', index=False)

# renewal index P(100 ps)/P(20 ps) and exponent alpha of the MSD
result.assign(P_ratio=(table['P', 100] / table['P', 20]).to_numpy()).to_csv('results/renewal_index.csv', index=False)
lags = np.array([20] + LAGS)
alpha = []
for system in systems:
    msd = table['MSD_nm2'].loc[system, lags].to_numpy(float)
    alpha.append(np.polyfit(np.log(lags), np.log(msd), 1)[0])
result.assign(alpha=alpha).to_csv('results/msd_exponent.csv', index=False)

# tests: the two halves of the run, and the curves up to the longest lag with all systems reliable
test = {'region': REGION, 'systems': len(systems)}
for name, half in [('first_half', range(1, 6)), ('second_half', range(6, 11))]:
    part = blocks[(blocks.region == REGION) & blocks.block.isin(half)].copy()
    part['squared_sum'] = part.MSD_nm2 * part.n_survivors
    sums = part.groupby(['L', 'coverage', 'lag_ps'])[['squared_sum', 'n_survivors', 'n_origins']].sum()
    half_table = pd.DataFrame({'MSD_nm2': sums.squared_sum / sums.n_survivors,
                               'P': sums.n_survivors / sums.n_origins}).unstack('lag_ps').loc[systems]
    half_regime = regimes(curve_shapes(half_table, LAGS))
    test[f'{name}_same'] = int((half_regime == regime).sum())
    test[f'{name}_ARI'] = adjusted_rand_score(regime, half_regime)
longest = 100
while longest < 200 and reliable.loc[systems, longest + 20].all():
    longest += 20
test['longest_lag_ps'] = longest
if longest > 100:
    longer_regime = regimes(curve_shapes(table, range(40, longest + 1, 20)))
    test['longer_same'] = int((longer_regime == regime).sum())
    test['longer_ARI'] = adjusted_rand_score(regime, longer_regime)
pd.DataFrame([test]).to_csv('results/region_tests.csv', index=False)
