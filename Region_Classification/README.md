# Tests of the two regimes

Tests of the two regimes of all first-shell water, from the MSD and *P* of every region of the first shell
(`../Region_MSD_and_P`). A system enters the classification when its values are reliable at every lag from 20 to
100 ps. The classification is k-means with two groups on the logarithms of MSD and *P* divided by their values at
20 ps, at 40, 60, 80 and 100 ps.

## Steps

| Step | What it does | Writes |
|---|---|---|
| `step1_classify_regions.py` | the regimes; silhouette score and bootstrap stability for 2, 3 and 4 groups; the regimes of the two halves of the run and with longer lags; the renewal index *P*(100 ps)/*P*(20 ps) and the exponent α of MSD ~ τ<sup>α</sup> | `results/k_scores.csv`, `renewal_index.csv`, `msd_exponent.csv`, `region_tests.csv` |
| `plot_figures.py` | the figure | `figures/water_regime_tests.pdf` |

The columns of every file are explained at the top of the scripts.

## Figure

<p align="center"><img src="figures/water_regime_tests.png" width="100%" alt="Tests of the two regimes"></p>

(a) Silhouette score and bootstrap stability for *k* = 2, 3 and 4 groups (shaded: *k* = 2). (b) Renewal index
*P*(100 ps)/*P*(20 ps) and (c) exponent α of MSD ~ τ<sup>α</sup> (20 to 100 ps) of every system over the chain lengths
and coverages; white lines: the border between the two regimes.

## Run

```bash
python3 step1_classify_regions.py      # after step 2 of ../Region_MSD_and_P
python3 plot_figures.py
```
