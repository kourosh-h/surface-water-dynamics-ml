"""
Settings of the analysis. This is the only file to change for your own systems: where the simulations are, which
systems to analyze, which atoms make the first shell, and the numbers of the method.

The values below are those of the example: self-assembled monolayers (SAMs) of alkylsilanes on both faces of a silica
slab in water, with chain lengths L of 4 to 8 carbon atoms and coverages of 0.08 to 3.9 molecules/nm^2.
"""

# ----------------------------------------------------------------------
# Systems
# ----------------------------------------------------------------------
SIMULATIONS = '../simulations'          # the folder that holds the folders of all systems
STRUCTURE = 'confout.gro'               # the structure file in the folder of every system
TRAJECTORY = 'traj_unwrap_whole.xtc'    # the trajectory with whole molecules in the folder of every system (step 0)
RUN_INPUT = 'npt.tpr'                   # step 0: the run input file of GROMACS
RAW_TRAJECTORY = 'traj_comp.xtc'        # step 0: the trajectory written by GROMACS

# Every system has a name (for its result files), its folder inside SIMULATIONS and the values that describe it. The
# names of these values are listed in LABELS; they are the first columns of the tables.
# Change this part to your own systems: the lists and the loop below make the systems of the example from the chain
# lengths and the coverages of its SAMs and from the names of their folders.
CHAIN_LENGTHS = [4, 5, 6, 7, 8]
COVERAGES = ['0.08', '0.45', '0.77', '1.20', '1.65', '2.0', '2.40', '2.92', '3.45', '3.9']    # as in the folder names

SYSTEMS = []
for L in CHAIN_LENGTHS:
    for coverage in COVERAGES:
        SYSTEMS.append({'name': f'L{L}_G{coverage}',
                        'folder': f'all_wet_production_run/chain_length_{L}/coverage_{coverage}/nptrunP',
                        'L': L,
                        'coverage': float(coverage)})
LABELS = ['L', 'coverage']

# ----------------------------------------------------------------------
# Atoms (MDAnalysis selections)
# ----------------------------------------------------------------------
WATER_OXYGEN = 'resname TIP3 and name OH2'      # the oxygen of every water molecule
SOLID = 'resname SUR'                            # the atoms of the solid, here the silica slab

# The faces of the solid whose water is analyzed: 'upper' (above the middle of the solid), 'lower' (below it) or both.
# The surface lies in the x-y plane.
FACES = ['upper', 'lower']


def rdf_groups(system):
    """The groups of atoms of the RDFs with the water oxygens."""
    L = system['L']
    groups = {}
    groups['CH3'] = f'resname OTS and name C{L}'                    # the last carbon of every chain
    groups['CH2'] = f'resname OTS and name C* and not name C{L}'    # the other carbons of the chains
    groups['OD'] = 'resname SUR and name OD'                        # the oxygens of the silanol groups of the silica
    groups['head_O'] = 'resname OTS and name O*'                    # the oxygens of the head groups
    return groups


# The contact distances of the first shell. Each one is the median over the systems of the first minima of the RDFs
# listed with it, rounded to 0.01 nm.
CONTACT_DISTANCES = {'carbon': ['CH3'], 'oxygen': ['OD', 'head_O']}


def shell_groups(system):
    """
    The groups of atoms that make the first shell, each with the name of its contact distance. A water molecule is in
    the first shell when its oxygen is within the contact distance of an atom of one of these groups. Its group is the
    one with the smallest distance divided by the contact distance.
    """
    L = system['L']
    groups = {}
    groups['CH3'] = (f'resname OTS and name C{L}', 'carbon')
    groups['CH2'] = (f'resname OTS and name C* and not name C{L}', 'carbon')
    groups['oxygen'] = ('(resname SUR or resname OTS) and name O*', 'oxygen')    # silica and head-group oxygens
    return groups


# The regions of the first shell, each with the groups whose water belongs to it. The region 'all' must stay.
REGIONS = {'all': ['CH3', 'CH2', 'oxygen'], 'oxygen': ['oxygen'], 'CH2': ['CH2'], 'CH3': ['CH3']}

# ----------------------------------------------------------------------
# Frames, lag times and blocks
# ----------------------------------------------------------------------
LAST_FRAMES = 4501      # the analyzed frames: the last 4501 frames of every trajectory (90 ns, frames 20 ps apart)
LAGS = range(1, 11)     # the lag times in frames: 20, 40, ..., 200 ps
N_BLOCKS = 10           # the analyzed frames in ten blocks of 9 ns, for the standard deviations
PROCESSES = 8           # the number of systems analyzed at the same time (CPU cores)

# ----------------------------------------------------------------------
# Step 1: RDFs
# ----------------------------------------------------------------------
RDF_EVERY = 10          # the RDFs use every 10th analyzed frame
RDF_BIN = 0.005         # nm, the width of the bins
RDF_RANGE = 1.0         # nm, the largest distance
PEAK_BEFORE = 0.6       # nm: the first peak of an RDF is looked for at shorter distances
MINIMUM_BEFORE = 0.8    # nm: the first minimum after the peak is looked for at shorter distances

# ----------------------------------------------------------------------
# Steps 5 and 6: regions and tables
# ----------------------------------------------------------------------
# A value of MSD and P of a region is reliable when at least MIN_SURVIVORS molecule - time origin pairs stay in the
# region until the lag, at least MIN_GOOD_BLOCKS blocks have MIN_BLOCK_SURVIVORS or more of them, and P is at least
# MIN_P. The standard deviation of the MSD of a region uses the blocks with MIN_BLOCK_SURVIVORS or more.
MIN_SURVIVORS = 200
MIN_BLOCK_SURVIVORS = 20
MIN_GOOD_BLOCKS = 3
MIN_P = 0.005
MEAN_OVER = 'L'                 # the first-shell water at each group is also averaged over this label (None: not)
COUNT_LAGS = [20, 40, 100]      # ps, the lags of the counts of the molecule - time origin pairs in contact

# ----------------------------------------------------------------------
# Steps 7 and 8: k-means clustering
# ----------------------------------------------------------------------
SHAPE_LAGS = [40, 60, 80, 100]  # ps, the lags of the shapes of the MSD and P curves
N_CLUSTERS = 2                  # the number of regimes; choose it from the scores of step 8
CLUSTER_NUMBERS = [2, 3, 4]     # step 8: the numbers of clusters compared
RESAMPLES = 200                 # step 8: the bootstrap resamples of the systems
SEED = 123                      # random seed of k-means and of the resamples

# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
X_LABEL = 'coverage'            # the label on the horizontal axis, and the colors of the curves of the systems
LINE_LABEL = 'L'                # one line (or row of panels) for every value of this label (None: no such lines)
LABEL_TITLES = {'L': '$L$', 'coverage': 'coverage (molecules/nm$^2$)'}     # how the labels are written
