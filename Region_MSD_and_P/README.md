# MSD and P of every region of the first shell

From the trajectories of the wet systems to the survival-conditioned lateral MSD and the contact survival *P* of the
first-shell water of every region of the first shell (all first-shell water; the water closest to the oxygens of the
silica and the head groups, to CH<sub>2</sub> and to CH<sub>3</sub>), per 9-ns block. The tables of step 2 are the
input of `../Region_Classification`.

## Steps

| Step | What it does | Writes |
|---|---|---|
| `step1_region_MSD_and_P.py` | MSD and *P* of every region, system, lag and block (MDAnalysis, 8 systems at a time; the cutoffs from step 1 of `../RDF_MSD_and_P`) | `results/per_system/` |
| `step2_region_tables.py` | all systems in two tables, with the reliability of every value | `results/region_curves.csv`, `results/region_blocks.csv` |
| `step3_water_counts.py` | the numbers of water molecules in the box, in the first shell and at each group, and of the molecule - origin pairs | `results/first_shell_water_counts.csv` |
| `plot_figures.py` | the figures of all first-shell water | `figures/` |

`systems.py` holds the layout of the folder of the simulations and finds the systems in it. Step 3 also reads the table
of the simulation systems, `results/table_simulation_systems.csv` of sam-system-summary (`SYSTEM_TABLE`). The columns of
every file are explained at the top of the scripts.

## Run

```bash
# after step 1 of ../RDF_MSD_and_P (the cutoffs of the first shell)
python3 step1_region_MSD_and_P.py <folder of the simulations>
python3 step2_region_tables.py
python3 step3_water_counts.py
python3 plot_figures.py
```
