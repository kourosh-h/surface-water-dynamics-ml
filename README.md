<p align="center">
  <img src="images/first_shell_steps.png" width="100%" alt="A snapshot of a wet system, the contact distance from the RDF, and the first hydration shell">
</p>

<h1 align="center">sam-water-dynamics</h1>

<p align="center">
  The first hydration shell of alkylsilane SAMs on silica: where the water sits, how it moves and how long it stays
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.9-0C54B4" alt="Python 3.9">
  <img src="https://img.shields.io/badge/MDAnalysis-2-E4840C" alt="MDAnalysis 2">
  <img src="https://img.shields.io/badge/scikit--learn-0.22-E4840C" alt="scikit-learn 0.22">
  <img src="https://img.shields.io/badge/license-MIT-555555" alt="MIT license">
</p>

<p align="center">
  <a href="#the-systems">Systems</a> ·
  <a href="#parts">Parts</a> ·
  <a href="#results">Results</a> ·
  <a href="#how-to-run">How to run</a> ·
  <a href="#inputs">Inputs</a> ·
  <a href="#requirements">Requirements</a>
</p>

---

This repository holds the analysis code of the water next to self-assembled monolayers (SAMs) of alkylsilanes on
hydroxylated silica. From the molecular dynamics trajectories of the wet systems, it finds the first hydration shell of
every system (picture above: the contact distance of a group of atoms is the first minimum of its radial distribution
function, and the shell is the union of the contact regions of all atoms), counts the water at each group of atoms,
measures how the water of the shell moves along the surface (survival-conditioned lateral mean squared displacement,
MSD) and how long it stays in contact with the surface (contact survival *P*), and sorts the systems into two regimes of
the water dynamics with k-means.

## The systems

Alkylsilane chains of length *L* = 4, 5, 6, 7 and 8 carbons, grafted on both faces of an α-quartz slab at ten coverages
Γ = 0.08 to 3.90 molecules/nm², in water (50 wet systems). The systems were built with
[apadana](https://github.com/kourosh-h/apadana-v4) and simulated with GROMACS. The code finds the systems from the
names of their folders, so it runs on any set of chain lengths and coverages.

## Parts

| Part | What it does |
|---|---|
| [RDF_MSD_and_P](RDF_MSD_and_P) | RDFs of the water oxygens, the cutoffs of the first shell, MSD and *P* of the first-shell water, and the water at each group of atoms |
| [Region_MSD_and_P](Region_MSD_and_P) | MSD and *P* of every region of the first shell (all water, and the water at CH<sub>3</sub>, CH<sub>2</sub> and the oxygens), per 9-ns block, with the reliability of every value |
| [ML_Classification](ML_Classification) | two regimes of the water dynamics, from the shapes of the MSD and *P* curves (k-means) |
| [Region_Classification](Region_Classification) | tests of the two regimes: the number of groups, bootstrap stability, the two halves of the run, longer lags, the renewal index and the MSD exponent |
| [First_Shell_Composition](First_Shell_Composition) | the share of the first-shell water at each group of atoms |

Every part has its own README with its steps, the files they write and its figures.

## Results

<table>
<tr>
<td width="38%"><img src="RDF_MSD_and_P/figures/first_shell_water_per_group.png" alt="First-shell water at each group of atoms"></td>
<td width="62%"><img src="ML_Classification/figures/water_dynamics_regimes.png" alt="Two regimes of the water dynamics"></td>
</tr>
<tr>
<td valign="top">First-shell water per nm² at the CH<sub>3</sub> and CH<sub>2</sub> groups and at the oxygens of the
silica and the head groups, mean over the chain lengths (band: SD over the chain lengths).</td>
<td valign="top">(a) The regime of the water dynamics of every system. (b), (c) MSD and <i>P</i> divided by their
values at 20 ps, mean over the systems of each regime (band: SD).</td>
</tr>
</table>

<p align="center"><img src="Region_MSD_and_P/figures/first_shell_water_msd.png" width="100%" alt="MSD of the first-shell water"></p>

Survival-conditioned lateral MSD of the first-shell water against the lag time, one panel per coverage and one line per
chain length (band: SD of the ten 9-ns blocks).

## How to run

Every part is a chain of numbered scripts that read the files of the steps before them; run each script from its own
folder. The steps that read the trajectories take the folder of the simulations as the argument:

```bash
# whole molecules: traj_unwrap_whole.xtc of every run
bash step0_unwrap.sh <folder of the simulations>

cd RDF_MSD_and_P
python3 step1_RDF.py <folder of the simulations>
python3 step2_MSD.py <folder of the simulations>
python3 step3_P.py <folder of the simulations>
python3 step4_first_shell_water.py <folder of the simulations>
python3 step5_tables.py
python3 step6_first_shell_means.py

cd ../Region_MSD_and_P
python3 step1_region_MSD_and_P.py <folder of the simulations>
python3 step2_region_tables.py
python3 step3_water_counts.py

cd ../RDF_MSD_and_P
python3 step7_tables.py

cd ../ML_Classification && python3 step1_classify_water_dynamics.py
cd ../Region_Classification && python3 step1_classify_regions.py
cd ../First_Shell_Composition && python3 step1_composition.py
```

Then `python3 plot_figures.py` in every part. The trajectory steps analyze 8 systems at a time (`PROCESSES` at the top
of these scripts).

## Inputs

The data files are not part of the repository. The folder of the simulations holds every production run as

```
all_wet_production_run/chain_length_<L>/coverage_<coverage>/nptrunP/
```

with `confout.gro` and `traj_unwrap_whole.xtc` (made by `step0_unwrap.sh` from `traj_comp.xtc` and `npt.tpr`); the
trajectory steps use the last 90 ns of the run, with frames 20 ps apart (`LAST_FRAMES` at the top of the scripts). For
another layout of the folders, change `run_folder()` and `find_systems()` in `systems.py`. The water counts of
`Region_MSD_and_P` also read the table of the simulation systems, `results/table_simulation_systems.csv` of
sam-system-summary (`SYSTEM_TABLE` at the top of `step3_water_counts.py`).

## Requirements

Python 3.9 with NumPy, SciPy, pandas, scikit-learn, Matplotlib (Arial font) and, for the steps that read the
trajectories, MDAnalysis; GROMACS for step 0.

```bash
pip install -r requirements.txt
```

The versions in `requirements.txt` are those the results were made with.

---

<p align="center">
  <sub>Written by Kourosh Hasheminejad and released under the <a href="LICENSE">MIT license</a>.</sub>
</p>
