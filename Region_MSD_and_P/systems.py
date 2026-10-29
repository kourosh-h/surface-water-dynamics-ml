"""
Where the production runs are, and which systems there are.

The folder of the simulations is given on the command line. Inside it every production run has its own folder:

    all_dry_production_run/chain_length_<L>/coverage_<coverage>/nvtrunP/     dry systems (NVT)
    all_wet_production_run/chain_length_<L>/coverage_<coverage>/nptrunP/     wet systems (NPT)

L is the number of carbon atoms of the alkyl chain and coverage the number of chains per nm^2 on each face of the
silica slab, as written in the folder name (for example 0.45 or 3.9). all_wet_production_run/chain_length_0/coverage_0
is the bare silica in water. For another layout of the folders, change run_folder() and find_systems().
"""
import sys
from pathlib import Path


def simulations_folder():
    """The folder of the simulations, from the command line."""
    if len(sys.argv) < 2:
        sys.exit(f'Run:  python3 {Path(sys.argv[0]).name} <folder of the simulations>')
    return Path(sys.argv[1])


def run_folder(simulations, state, L, coverage):
    """The folder of one production run: dry runs NVT (nvtrunP), wet runs NPT (nptrunP)."""
    if state == 'dry':
        return simulations / 'all_dry_production_run' / f'chain_length_{L}' / f'coverage_{coverage}' / 'nvtrunP'
    return simulations / 'all_wet_production_run' / f'chain_length_{L}' / f'coverage_{coverage}' / 'nptrunP'


def find_systems(simulations, state, bare_silica=False):
    """
    The systems of one state (dry or wet) found in the folder of the simulations, as (L, coverage) with the coverage
    as written in the folder name, sorted by chain length and coverage. The bare silica (L = 0) is left out unless
    bare_silica is True.
    """
    found = []
    for folder in (simulations / f'all_{state}_production_run').glob('chain_length_*/coverage_*'):
        L = int(folder.parent.name.replace('chain_length_', ''))
        coverage = folder.name.replace('coverage_', '')
        if L > 0 or bare_silica:
            found.append((L, float(coverage), coverage))
    found.sort()
    return [(L, coverage) for L, value, coverage in found]


def systems_in(folder):
    """
    The systems of the result files L<L>_G<coverage>.csv in a folder, as (L, coverage) with the coverage as written in
    the file name, sorted by chain length and coverage. Files with more after the coverage (L4_G0.08_blocks.csv) are
    left out.
    """
    found = []
    for path in Path(folder).glob('L*_G*.csv'):
        L, coverage = path.stem[1:].split('_G')         # L4_G0.08 -> 4, 0.08
        if '_' in coverage:
            continue
        found.append((int(L), float(coverage), coverage))
    found.sort()
    return [(L, coverage) for L, value, coverage in found]
