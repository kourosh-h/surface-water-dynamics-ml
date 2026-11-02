"""
The first shell of the water, the same in steps 2 to 5: open the trajectory of a system, and find in one frame the
water molecules of the first shell, their face and their group.
"""
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

import settings


def open_system(system):
    """
    The trajectory of one system: the universe, the water oxygens, the atoms of the solid and the groups of the first
    shell as a list of (name, atoms, contact distance in Angstrom, the length unit of MDAnalysis).
    """
    folder = Path(settings.SIMULATIONS) / system['folder']
    u = mda.Universe(str(folder / settings.STRUCTURE), str(folder / settings.TRAJECTORY))
    water = u.select_atoms(settings.WATER_OXYGEN)
    solid = u.select_atoms(settings.SOLID)
    cutoffs = pd.read_csv('results/cutoffs.csv', index_col='atoms').cutoff_nm      # from step 1
    groups = []
    for name, (selection, distance) in settings.shell_groups(system).items():
        groups.append((name, u.select_atoms(selection), 10 * float(cutoffs[distance])))
    return u, water, solid, groups


def closest_atom(xyz, atoms, box, cutoff):
    """
    For every water oxygen: the distance to the closest atom of the group (Angstrom; infinite if no atom is within the
    cutoff) and the number of that atom in the group.
    """
    pairs, distances = capped_distance(xyz, atoms.positions, max_cutoff=cutoff, box=box)
    order = np.lexsort((distances, pairs[:, 0]))           # the pairs by water and, for every water, by distance
    waters, first = np.unique(pairs[order, 0], return_index=True)
    closest = order[first]                                 # the pair of every water with its closest atom
    distance = np.full(len(xyz), np.inf, dtype=np.float32)
    distance[waters] = distances[closest]
    atom = np.zeros(len(xyz), dtype=int)
    atom[waters] = pairs[closest, 1]
    return distance, atom


def faces_and_groups(xyz, groups, middle, box):
    """
    For every water oxygen: its face (1 upper, 2 lower, 0 not in the first shell) and its group (the number of the
    group in groups, -1 for none). The group is the one with the smallest distance / contact distance (with a
    tolerance for the rounding of the numbers), and the face that of its closest atom, above or below the middle of
    the solid. Water at a face that is not in settings.FACES is left out of the first shell.
    """
    face = np.zeros(len(xyz), dtype=int)
    group = np.full(len(xyz), -1)
    smallest = np.full(len(xyz), np.inf)
    for g, (name, atoms, cutoff) in enumerate(groups):
        distance, atom = closest_atom(xyz, atoms, box, cutoff)
        z = atoms.positions[atom, 2] - middle
        z -= box[2] * np.rint(z / box[2])                   # the closest periodic image
        ratio = distance / cutoff
        closer = (ratio <= 1.0 + 1e-7) & (ratio < smallest)
        face[closer] = np.where(z >= 0, 1, 2)[closer]
        group[closer] = g
        smallest[closer] = ratio[closer]
    if 'upper' not in settings.FACES:
        group[face == 1] = -1
        face[face == 1] = 0
    if 'lower' not in settings.FACES:
        group[face == 2] = -1
        face[face == 2] = 0
    return face, group
