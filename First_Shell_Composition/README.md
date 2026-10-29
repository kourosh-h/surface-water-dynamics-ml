# Composition of the first hydration shell

The share of the first-shell water of every wet system at each group of atoms: the terminal carbons (CH<sub>3</sub>),
the other chain carbons (CH<sub>2</sub>) and the oxygens of the silica and the head groups. Every water molecule of the
shell belongs to the group of its closest atom (step 4 of `../RDF_MSD_and_P`).

## Steps

| Step | What it does | Writes |
|---|---|---|
| `step1_composition.py` | the share of every group, and the group with the most first-shell water (input: `results/first_shell_water.csv` of `../RDF_MSD_and_P`) | `results/largest_group.csv`, `share_CH3.csv`, `share_CH2.csv`, `share_oxygen.csv` |
| `plot_figures.py` | the figure | `figures/first_shell_composition.pdf` |

The columns of every file are explained at the top of the scripts.

## Figure

<p align="center"><img src="figures/first_shell_composition.png" width="100%" alt="Composition of the first hydration shell"></p>

(a) The group closest to most of the first-shell water (color), shaded by its share of the first-shell water (below
1/2, 1/2 to 2/3, above 2/3), interpolated between the systems; thick lines where the largest group changes. (b)-(d) The
fraction of the first-shell water closest to CH<sub>3</sub>, CH<sub>2</sub> and the oxygens, one line per chain length.

## Run

```bash
python3 step1_composition.py      # after step 6 of ../RDF_MSD_and_P
python3 plot_figures.py
```
