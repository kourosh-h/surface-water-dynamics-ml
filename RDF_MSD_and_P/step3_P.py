"""
Contact survival P(tau) of the first-shell water of every wet SAM system, at the lags tau = 20, 40, ..., 200 ps.

First shell: a water molecule is in the first shell when its oxygen is within the carbon cutoff of a chain carbon (CH3
or CH2) or within the oxygen cutoff of an oxygen of the silica or of the head groups (cutoffs from step1_RDF.py). Each
face of the slab has its own first shell: a water molecule in the shell is on the face of its closest atom, of the
group with the smallest distance / cutoff.

P(tau): over all time origins in the last 90 ns of the production run (frames 20 ps apart), the fraction of the water
molecules in the first shell of a face at the origin that stay in it, without a break, until the origin + tau; both
faces together.

Output: results/P/L<L>_G<coverage>.csv: lag_ps, P
Run:  python3 step3_P.py <folder of the simulations>      (after step1_RDF.py)
"""
from multiprocessing import Pool
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from systems import find_systems, run_folder, simulations_folder

PROCESSES = 8                   # systems analyzed at the same time
LAST_FRAMES = 4501              # the last 90 ns of the run, frames 20 ps apart
LAGS = range(1, 11)             # lags in frames: 20 ... 200 ps
# the cutoffs of the first shell, in Angstrom (the length unit of MDAnalysis)
cutoffs = pd.read_csv('results/RDF/cutoffs.csv', index_col='atoms').cutoff_nm
CARBON_CUTOFF = 10 * float(cutoffs['carbon'])
OXYGEN_CUTOFF = 10 * float(cutoffs['oxygen'])


def closest_atom(xyz, atoms, box, cutoff):
    """For every water oxygen: the distance to the closest atom of the group (Angstrom; infinite if none is within the
    cutoff) and the number of that atom in the group."""
    pairs, distances = capped_distance(xyz, atoms.positions, max_cutoff=cutoff, box=box)
    order = np.lexsort((distances, pairs[:, 0]))           # the pairs by water and, for every water, by distance
    waters, first = np.unique(pairs[order, 0], return_index=True)
    closest = order[first]                                 # the pair of every water with its closest atom
    distance = np.full(len(xyz), np.inf, dtype=np.float32)
    distance[waters] = distances[closest]
    atom = np.zeros(len(xyz), dtype=int)
    atom[waters] = pairs[closest, 1]
    return distance, atom


def shell_faces(xyz, groups, middle, box):
    """The first shell of every water oxygen: 1 on the upper face of the slab, 2 on the lower face, 0 not in the shell.
    For every group of atoms: the closest atom, its distance divided by the cutoff (with a tolerance for the rounding
    of the numbers) and its face, above or below the middle of the silica; a water in the shell takes the face of the
    group with the smallest distance / cutoff."""
    face = np.zeros(len(xyz), dtype=int)
    smallest = np.full(len(xyz), np.inf)
    for atoms, cutoff in groups:
        distance, atom = closest_atom(xyz, atoms, box, cutoff)
        z = atoms.positions[atom, 2] - middle
        z -= box[2] * np.rint(z / box[2])
        ratio = distance / cutoff
        closer = (ratio <= 1.0 + 1e-7) & (ratio < smallest)
        face[closer] = np.where(z >= 0, 1, 2)[closer]
        smallest[closer] = ratio[closer]
    return face


def P_of_system(simulations, L, coverage):
    """P of one system at every lag, written to results/P/L<L>_G<coverage>.csv."""
    folder = run_folder(simulations, 'wet', L, coverage)
    u = mda.Universe(str(folder / 'confout.gro'), str(folder / 'traj_unwrap_whole.xtc'))
    water = u.select_atoms('resname TIP3 and name OH2')                       # the water oxygens
    silica = u.select_atoms('resname SUR')
    groups = [(u.select_atoms(f'resname OTS and name C{L}'), CARBON_CUTOFF),                  # CH3
              (u.select_atoms(f'resname OTS and name C* and not name C{L}'), CARBON_CUTOFF),  # CH2
              (u.select_atoms('(resname SUR or resname OTS) and name O*'), OXYGEN_CUTOFF)]     # silica and head-group O

    times = []                          # time of every frame (ns)
    in_shell = []                       # number of water molecules in the first shell, in every frame
    n_survivors = [0] * len(LAGS)       # summed over the frames: water molecules in the shell since lag frames before
    face = np.zeros(len(water), dtype=int)
    run_length = np.zeros(len(water), dtype=int)   # frames in the shell of the same face without a break, up to now
    for ts in u.trajectory[-LAST_FRAMES:]:
        times.append(ts.time / 1000.0)
        box = ts.dimensions.astype(np.float64)
        face_before = face
        face = shell_faces(water.positions, groups, silica.positions[:, 2].mean(), box)
        stayed = (face != 0) & (face == face_before)
        run_length = np.where(stayed, run_length + 1, np.where(face != 0, 1, 0))
        in_shell.append(np.count_nonzero(face))
        for i, lag in enumerate(LAGS):
            n_survivors[i] += np.count_nonzero(run_length >= lag + 1)

    # P: the survivors divided by the water molecules in the shell at the origins (all frames but the last lag ones)
    dt = float(np.median(np.diff(times)))          # time between the frames (ns)
    rows = []
    for i, lag in enumerate(LAGS):
        rows.append({'lag_ps': lag * dt * 1000, 'P': n_survivors[i] / sum(in_shell[:len(in_shell) - lag])})
    pd.DataFrame(rows).to_csv(f'results/P/L{L}_G{coverage}.csv', index=False)


if __name__ == '__main__':
    simulations = simulations_folder()
    Path('results/P').mkdir(parents=True, exist_ok=True)
    systems = [(simulations, L, coverage) for L, coverage in find_systems(simulations, 'wet')]
    with Pool(PROCESSES) as pool:
        pool.starmap(P_of_system, systems)
