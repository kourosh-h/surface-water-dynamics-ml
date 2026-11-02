"""
Step 1: radial distribution functions (RDFs) of the water oxygens around the groups of atoms of rdf_groups() in
settings.py, their first peaks and minima, and the contact distances of the first shell. The number of water molecules
of every system is written too.

The RDFs are taken from every RDF_EVERY-th analyzed frame, with bins of RDF_BIN up to RDF_RANGE, and normalized as in
gmx rdf (by the number of pairs per box volume). Of every RDF, smoothed by a running mean over five bins: the first
peak (below PEAK_BEFORE) and the first minimum after it (below MINIMUM_BEFORE).

Contact distances: for each entry of CONTACT_DISTANCES, the median over the systems of the first minima of its RDFs,
rounded to 0.01 nm.

Output:
    results/rdf/<system>.csv: r_nm and one column per group of atoms
    results/rdf/first_minima.csv: the labels of the system, pair, first_peak_nm, first_minimum_nm
    results/cutoffs.csv: atoms (the name of the contact distance), cutoff_nm
    results/water_molecules.csv: system, the labels of the system, water_molecules (in the box)
Run:  python3 step1_rdf.py
"""
from multiprocessing import Pool
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

import settings

EDGES = np.arange(0.0, settings.RDF_RANGE + 1e-9, settings.RDF_BIN)       # edges of the RDF bins (nm)


def first_peak_and_minimum(r, g):
    """The first peak of g(r) and the first minimum after it, from g(r) smoothed over five bins."""
    smooth = np.convolve(g, np.ones(5) / 5, mode='same')
    peak = np.argmax(np.where(r < settings.PEAK_BEFORE, smooth, -np.inf))
    after = np.flatnonzero((r > r[peak]) & (r < settings.MINIMUM_BEFORE))
    minimum = after[np.argmin(smooth[after])]
    return r[peak], r[minimum]


def rdf_of_system(system):
    """
    The RDFs of one system (written to results/rdf/<system>.csv), their first peaks and minima, and the number of
    water molecules.
    """
    folder = Path(settings.SIMULATIONS) / system['folder']
    u = mda.Universe(str(folder / settings.STRUCTURE), str(folder / settings.TRAJECTORY))
    water = u.select_atoms(settings.WATER_OXYGEN)
    groups = {}
    for name, selection in settings.rdf_groups(system).items():
        groups[name] = u.select_atoms(selection)
    counts = {}             # pairs per distance bin
    density = {}            # number of pairs per box volume, summed over the frames
    for name in groups:
        counts[name] = np.zeros(len(EDGES) - 1)
        density[name] = 0.0
    for ts in u.trajectory[-settings.LAST_FRAMES::settings.RDF_EVERY]:
        box = ts.dimensions.astype(np.float64)
        volume = box[0] * box[1] * box[2] / 1000.0          # nm^3
        for name, atoms in groups.items():
            pairs, distances = capped_distance(atoms.positions, water.positions, max_cutoff=10 * settings.RDF_RANGE,
                                               box=box)
            counts[name] += np.histogram(distances / 10.0, bins=EDGES)[0]
            density[name] += len(atoms) * len(water) / volume

    # pairs per bin divided by the pairs per volume and by the volume of the spherical shell of the bin
    r = 0.5 * (EDGES[:-1] + EDGES[1:])
    shell_volume = 4.0 / 3.0 * np.pi * (EDGES[1:] ** 3 - EDGES[:-1] ** 3)
    rdf = pd.DataFrame({'r_nm': r})
    minima = []
    for name in groups:
        rdf[name] = counts[name] / (density[name] * shell_volume)
        peak, minimum = first_peak_and_minimum(r, rdf[name].to_numpy())
        row = {}
        for label in settings.LABELS:
            row[label] = system[label]
        row['pair'] = f'{name}-OW'
        row['first_peak_nm'] = peak
        row['first_minimum_nm'] = minimum
        minima.append(row)
    rdf.to_csv(f'results/rdf/{system["name"]}.csv', index=False)
    return pd.DataFrame(minima), len(water)


if __name__ == '__main__':
    Path('results/rdf').mkdir(parents=True, exist_ok=True)
    with Pool(settings.PROCESSES) as pool:
        results = pool.map(rdf_of_system, settings.SYSTEMS)

    minima = []
    water_molecules = []
    for system, (system_minima, n_water) in zip(settings.SYSTEMS, results):
        minima.append(system_minima)
        row = {'system': system['name']}
        for label in settings.LABELS:
            row[label] = system[label]
        row['water_molecules'] = n_water
        water_molecules.append(row)
    minima = pd.concat(minima)
    minima.to_csv('results/rdf/first_minima.csv', index=False)
    pd.DataFrame(water_molecules).to_csv('results/water_molecules.csv', index=False)

    # contact distances of the first shell: the median first minimum over the systems, rounded to 0.01 nm
    names = []
    cutoffs = []
    for name, rdf_names in settings.CONTACT_DISTANCES.items():
        pairs = []
        for rdf_name in rdf_names:
            pairs.append(f'{rdf_name}-OW')
        names.append(name)
        cutoffs.append(round(minima[minima.pair.isin(pairs)].first_minimum_nm.median(), 2))
    pd.DataFrame({'atoms': names, 'cutoff_nm': cutoffs}).to_csv('results/cutoffs.csv', index=False)
