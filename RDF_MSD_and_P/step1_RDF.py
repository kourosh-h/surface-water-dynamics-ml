"""
Radial distribution functions (RDFs) of the water oxygens around four groups of atoms of every wet SAM system, their
first minima, and the cutoffs of the first hydration shell.

Groups: the terminal carbons of the chains (CH3), the other chain carbons (CH2), the silanol oxygens of the silica (OD)
and the oxygens of the head groups (head_O). The RDFs are taken from every 10th frame (every 200 ps) of the last 90 ns
of the production run, with 0.005-nm bins up to 1.0 nm, and normalized as in gmx rdf (by the number of pairs per box
volume). Of every RDF, smoothed by a running mean over five bins: the first peak (r < 0.6 nm) and the first minimum
after it (r < 0.8 nm).

Cutoffs of the first shell: the median first minimum over the systems, rounded to 0.01 nm. One cutoff for all chain
carbons, from the CH3 RDFs (the first minimum of CH2 varies more between the systems), and one for all oxygens, from the
OD and head_O RDFs.

Output (folder results/RDF):
    L<L>_G<coverage>.csv: the RDFs of one system: r_nm, CH3, CH2, OD, head_O
    first_minima.csv: L, coverage, pair, first_peak_nm, first_minimum_nm
    cutoffs.csv: atoms (carbon or oxygen), cutoff_nm

Run:  python3 step1_RDF.py <folder of the simulations>
"""
from multiprocessing import Pool
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from systems import find_systems, run_folder, simulations_folder

PROCESSES = 8                                   # systems analyzed at the same time
LAST_FRAMES = 4501                              # the last 90 ns of the run, frames 20 ps apart
RDF_EVERY = 10                                  # the RDFs use every 10th of these frames
EDGES = np.arange(0.0, 1.0 + 1e-9, 0.005)       # edges of the RDF bins (nm)


def first_peak_and_minimum(r, g):
    """First peak (r < 0.6 nm) and the first minimum after it (r < 0.8 nm) of g(r), smoothed over five bins."""
    smooth = np.convolve(g, np.ones(5) / 5, mode='same')
    peak = np.argmax(np.where(r < 0.6, smooth, -np.inf))
    after = np.flatnonzero((r > r[peak]) & (r < 0.8))
    minimum = after[np.argmin(smooth[after])]
    return r[peak], r[minimum]


def rdf_of_system(simulations, L, coverage):
    """The four RDFs of one system (written to results/RDF/L<L>_G<coverage>.csv) and their first peaks and minima."""
    folder = run_folder(simulations, 'wet', L, coverage)
    u = mda.Universe(str(folder / 'confout.gro'), str(folder / 'traj_unwrap_whole.xtc'))
    water = u.select_atoms('resname TIP3 and name OH2')                       # the water oxygens
    groups = {'CH3': u.select_atoms(f'resname OTS and name C{L}'),            # terminal carbon of every chain
              'CH2': u.select_atoms(f'resname OTS and name C* and not name C{L}'),
              'OD': u.select_atoms('resname SUR and name OD'),
              'head_O': u.select_atoms('resname OTS and name O*')}
    counts = {name: np.zeros(len(EDGES) - 1) for name in groups}     # pairs per distance bin
    density = {name: 0.0 for name in groups}       # number of pairs per box volume, summed over the frames
    for ts in u.trajectory[-LAST_FRAMES::RDF_EVERY]:
        box = ts.dimensions.astype(np.float64)
        volume = box[0] * box[1] * box[2] / 1000.0          # nm^3
        for name, atoms in groups.items():
            pairs, distances = capped_distance(atoms.positions, water.positions, max_cutoff=10.0, box=box)
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
        minima.append({'L': L, 'coverage': float(coverage), 'pair': f'{name}-OW', 'first_peak_nm': peak,
                       'first_minimum_nm': minimum})
    rdf.to_csv(f'results/RDF/L{L}_G{coverage}.csv', index=False)
    return pd.DataFrame(minima)


if __name__ == '__main__':
    simulations = simulations_folder()
    Path('results/RDF').mkdir(parents=True, exist_ok=True)
    systems = [(simulations, L, coverage) for L, coverage in find_systems(simulations, 'wet')]
    with Pool(PROCESSES) as pool:
        minima = pd.concat(pool.starmap(rdf_of_system, systems))
    minima.to_csv('results/RDF/first_minima.csv', index=False)

    # cutoffs of the first shell: the median first minimum over the systems, rounded to 0.01 nm
    carbon = minima[minima.pair == 'CH3-OW'].first_minimum_nm.median()
    oxygen = minima[minima.pair.isin(['OD-OW', 'head_O-OW'])].first_minimum_nm.median()
    cutoffs = pd.DataFrame({'atoms': ['carbon', 'oxygen'], 'cutoff_nm': [round(carbon, 2), round(oxygen, 2)]})
    cutoffs.to_csv('results/RDF/cutoffs.csv', index=False)
