"""
Figures of the first-shell water and of the RDFs of the water oxygens, from the result files of steps 1 and 7.

    figures/first_shell_water_per_group.pdf   first-shell water per nm^2 at the CH3 and CH2 groups and at the oxygens of
                                              the silica and the head groups against the coverage, mean over the chain
                                              lengths (band: SD over the chain lengths)
    figures/water_rdfs.pdf                    RDFs g(r) between the water oxygen and the terminal CH3 carbons, the CH2
                                              carbons, the silanol oxygens of the silica and the head-group oxygens,
                                              one row per chain length and one color per coverage

Every figure is also written as PNG.
Run:  python3 plot_figures.py      (after steps 1 and 7)
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')                   # write files only, no window
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

from systems import systems_in

COVERAGE_TITLE = r'$\Gamma$ (molecules/nm$^2$)'
# the groups of the first-shell water: name in the result files, color, marker and label
GROUPS = [('CH3', '#7B4FA0', '^', r'at CH$_3$'), ('CH2', '#5F6B73', 'o', r'at CH$_2$'),
          ('oxygen', '#D2473A', 's', 'at silica and\nhead-group O')]
# the atom pairs of the RDFs: column in the RDF files and title
PAIRS = [('CH3', r'C$_{\rm CH_3}$–O$_{\rm w}$'), ('CH2', r'C$_{\rm CH_2}$–O$_{\rm w}$'),
         ('OD', r'O$_{\rm silanol}$–O$_{\rm w}$'), ('head_O', r'O$_{\rm head}$–O$_{\rm w}$')]
R_RANGE = (0.2, 1.0)                    # nm, r range of the RDFs

# Fonts, lines and ticks
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'pdf.fonttype': 42, 'mathtext.default': 'it',
    'axes.unicode_minus': False, 'axes.labelsize': 8.0, 'xtick.labelsize': 8.0, 'ytick.labelsize': 8.0,
    'legend.fontsize': 8.0, 'legend.title_fontsize': 8.0, 'axes.linewidth': 1.2, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.labelpad': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 4.2, 'ytick.major.size': 4.2, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
    'xtick.minor.size': 2.0, 'ytick.minor.size': 2.0, 'xtick.minor.width': 0.7, 'ytick.minor.width': 0.7,
    'lines.linewidth': 1.45, 'lines.markeredgewidth': 0.7, 'lines.dashed_pattern': [4.0, 2.2],
    'lines.scale_dashes': False, 'legend.frameon': False, 'legend.handlelength': 1.6, 'legend.handletextpad': 0.4,
    'legend.columnspacing': 0.9, 'legend.labelspacing': 0.28})


def coverage_colors(coverages):
    """One color per coverage, in the order of the coverages: the matplotlib tab colors (a color map for more than
    ten coverages)."""
    if len(coverages) <= 10:
        colors = plt.get_cmap('tab10').colors
    else:
        colors = plt.get_cmap('viridis')([k / (len(coverages) - 1) for k in range(len(coverages))])
    return {coverage: colors[k] for k, coverage in enumerate(coverages)}


def save(fig, name):
    """The figure as PDF and PNG in figures/."""
    fig.savefig(f'figures/{name}.pdf', bbox_inches='tight')
    fig.savefig(f'figures/{name}.png', bbox_inches='tight', dpi=200)
    plt.close(fig)


Path('figures').mkdir(exist_ok=True)

# first-shell water at each group, mean over the chain lengths
data = pd.read_csv('results/first_shell_water_mean_over_L.csv')
fig, ax = plt.subplots(figsize=(2.6, 2.0))
for group, color, marker, label in GROUPS:
    rows = data[data.group == group]
    ax.fill_between(rows.coverage, rows.water_per_nm2 - rows.water_per_nm2_sd,
                    rows.water_per_nm2 + rows.water_per_nm2_sd, color=color, alpha=0.18, linewidth=0, zorder=2)
    ax.plot(rows.coverage, rows.water_per_nm2, color=color, ls='--', marker=marker, ms=4.25, zorder=3, label=label)
ax.set_xlabel(COVERAGE_TITLE)
ax.set_ylabel(r'water (molecules/nm$^2$)')
ax.set_ylim(bottom=0)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_title('first-shell water, mean over the chain lengths', fontsize=8.0, pad=4)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.1), ncol=3)
save(fig, 'first_shell_water_per_group')

# RDFs: one row per chain length, one column per atom pair; the panels of a column share the g(r) range
systems = systems_in('results/RDF')
chain_lengths = []
coverages = []
for L, coverage in systems:
    if L not in chain_lengths:
        chain_lengths.append(L)
    if float(coverage) not in coverages:
        coverages.append(float(coverage))
colors = coverage_colors(sorted(coverages))
fig, axes = plt.subplots(len(chain_lengths), len(PAIRS), figsize=(7.0, 1.25 * len(chain_lengths) + 0.6),
                         sharex=True, sharey='col', squeeze=False)
fig.subplots_adjust(wspace=0.32, hspace=0.12)
largest = [0.0] * len(PAIRS)            # the largest g(r) of every column beyond the start of the r range
for L, coverage in systems:
    rdf = pd.read_csv(f'results/RDF/L{L}_G{coverage}.csv')
    row = chain_lengths.index(L)
    for column, (pair, title) in enumerate(PAIRS):
        axes[row, column].plot(rdf.r_nm, rdf[pair], color=colors[float(coverage)])
        largest[column] = max(largest[column], rdf[rdf.r_nm > R_RANGE[0]][pair].max())
for column, (pair, title) in enumerate(PAIRS):
    axes[0, column].set_ylim(0, 1.05 * largest[column])
    axes[0, column].set_title(title, fontsize=8.0, pad=3)
    axes[-1, column].set_xlabel('$r$ (nm)')
for row, L in enumerate(chain_lengths):
    axes[row, 0].set_ylabel('$g(r)$')
    axes[row, -1].text(1.05, 0.5, f'$L$ = {L}', transform=axes[row, -1].transAxes, ha='left', va='center')
for ax in axes.flat:
    ax.set_xlim(R_RANGE)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
entries = []
for coverage, color in colors.items():
    entries.append(Line2D([], [], color=color, label=f'{coverage:.2f}'))
fig.legend(handles=entries, loc='lower center', bbox_to_anchor=(0.5, 0.9), ncol=len(entries), title=COVERAGE_TITLE)
save(fig, 'water_rdfs')
