# Two regimes of the water dynamics

The two regimes of the first-shell water dynamics of the wet systems, from the shapes of their MSD and *P* curves. The
MSD and *P* of every system are divided by their values at 20 ps; the logarithms of these ratios at 40, 60, 80 and
100 ps, scaled to mean 0 and SD 1, are grouped by k-means with two groups (scikit-learn). Regime 0 (faster renewal) is
the group of the lower mean coverage.

## Steps

| Step | What it does | Writes |
|---|---|---|
| `step1_classify_water_dynamics.py` | the regime of every system; the curves with their regimes; the mean and SD of the curves of each regime (input: `results/msd_and_P.csv` of `../RDF_MSD_and_P`) | `results/curves_MSD_rel.csv`, `curves_P_rel.csv`, `regime_map.csv`, `regime_mean_MSD_rel.csv`, `regime_mean_P_rel.csv` |
| `plot_figures.py` | the figure | `figures/water_dynamics_regimes.pdf` |

The columns of every file are explained at the top of the scripts.

## Figure

<p align="center"><img src="figures/water_dynamics_regimes.png" width="100%" alt="Two regimes of the water dynamics"></p>

(a) The regime of the water dynamics of every system over the chain lengths and coverages. (b), (c) Survival-conditioned
lateral MSD and contact survival *P*(τ) of the first-shell water, each divided by its value at 20 ps: mean over the
systems of each regime (line) and SD (band).

## Run

```bash
python3 step1_classify_water_dynamics.py      # after step 5 of ../RDF_MSD_and_P
python3 plot_figures.py
```
