"""
Figures of the results, written to figures/ as PDF and PNG:

    rdfs                 the RDFs of the water oxygens around every group of atoms of step 1, with the contact
                         distances (dashed); one row of panels for every value of LINE_LABEL, one color for every value
                         of X_LABEL
    first_shell_water    the first-shell water per nm^2 at each group (top) and its share of the first shell (bottom)
                         against X_LABEL, one line for every value of LINE_LABEL (bars: standard deviation of the
                         blocks)
    msd_and_survival     the lateral MSD (top) and the contact survival P (bottom) of all first-shell water against the
                         lag time, one panel for every value of LINE_LABEL, one color for every value of X_LABEL (bars:
                         standard deviation of the blocks; reliable values only)
    regimes              (a) the regime of every system; (b), (c) MSD and P divided by their values at the first lag,
                         mean over the systems of each regime (band: standard deviation)
    cluster_tests        (a) silhouette score and bootstrap stability against the number of clusters; (b) the renewal
                         index against the MSD exponent of every system, colored by its regime

Run:  python3 plot_figures.py      (after steps 7 and 8)
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')                   # write files only, no window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FormatStrFormatter, LogLocator, NullFormatter

import settings

REGIME_COLORS = ['#2BA58A', '#1F5F99']
REGIME_NAMES = ['faster renewal', 'slower renewal']
LINE_COLORS = ['#4DBBD5', '#E68613', '#756BB1', '#B7371B', '#0D4F8B']

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'pdf.fonttype': 42, 'mathtext.default': 'it',
    'axes.unicode_minus': False, 'axes.labelsize': 8.0, 'axes.titlesize': 8.0, 'xtick.labelsize': 8.0,
    'ytick.labelsize': 8.0, 'legend.fontsize': 8.0, 'legend.title_fontsize': 8.0, 'axes.linewidth': 1.0,
    'axes.spines.top': False, 'axes.spines.right': False, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'legend.frameon': False, 'lines.linewidth': 1.2, 'lines.markersize': 3.5})


def save(fig, name):
    """The figure as PDF and PNG in figures/."""
    fig.savefig(f'figures/{name}.pdf', bbox_inches='tight')
    fig.savefig(f'figures/{name}.png', bbox_inches='tight', dpi=200)
    plt.close(fig)


def values_of(label):
    """The values of a label over all systems, sorted, each once ([None] for no label)."""
    if label is None:
        return [None]
    values = []
    for system in settings.SYSTEMS:
        if system[label] not in values:
            values.append(system[label])
    return sorted(values)


def colors_of(values, palette=None):
    """One color for every value, in order: the palette (or the matplotlib tab colors), a color map for more values."""
    if palette is None:
        palette = plt.get_cmap('tab10').colors
    if len(values) > len(palette):
        palette = plt.get_cmap('viridis')(np.linspace(0, 0.9, len(values)))
    colors = {}
    for k, value in enumerate(values):
        colors[value] = palette[k]
    return colors


def same_value(column, value):
    """True for the rows of a table column that hold the value (numbers with a decimal point within a tolerance)."""
    if isinstance(value, float):
        return np.isclose(column.to_numpy(float), value)
    return (column == value).to_numpy()


def rows_of(table, system):
    """The rows of a table that belong to one system: those with all its labels."""
    keep = np.ones(len(table), dtype=bool)
    for label in settings.LABELS:
        keep &= same_value(table[label], system[label])
    return table[keep]


def value_text(value):
    """A value of a label as text: numbers without trailing zeros (0.45, 2), other values as they are."""
    if isinstance(value, (int, float, np.integer, np.floating)):
        return f'{value:g}'
    return str(value)


def line_name(value):
    """The name of a value of LINE_LABEL in the legends, for example L = 4."""
    return f'{settings.LABEL_TITLES[settings.LINE_LABEL]} = {value_text(value)}'


def color_legend(fig, colors, title, y=1.0):
    """A legend above the figure: one colored line for every value."""
    entries = []
    for value, color in colors.items():
        entries.append(Line2D([], [], color=color, label=value_text(value)))
    fig.legend(handles=entries, loc='lower center', bbox_to_anchor=(0.5, y), ncol=min(len(entries), 10), title=title,
               handlelength=1.4, columnspacing=1.0)


Path('figures').mkdir(exist_ok=True)
x_values = values_of(settings.X_LABEL)
line_values = values_of(settings.LINE_LABEL)
x_colors = colors_of(x_values)
line_colors = colors_of(line_values, LINE_COLORS)
x_title = settings.LABEL_TITLES[settings.X_LABEL]
group_names = []
for name in settings.shell_groups(settings.SYSTEMS[0]):
    group_names.append(name)

# ----------------------------------------------------------------------
# RDFs, with the contact distances
# ----------------------------------------------------------------------
rdf_names = []
for name in settings.rdf_groups(settings.SYSTEMS[0]):
    rdf_names.append(name)
cutoffs = pd.read_csv('results/cutoffs.csv', index_col='atoms').cutoff_nm
rdf_cutoff = {}                         # the contact distance of every RDF that sets one
for distance, names in settings.CONTACT_DISTANCES.items():
    for name in names:
        rdf_cutoff[name] = cutoffs[distance]
size = (1.75 * len(rdf_names), 1.2 * len(line_values) + 0.6)
fig, axes = plt.subplots(len(line_values), len(rdf_names), figsize=size, sharex=True, sharey='col', squeeze=False)
fig.subplots_adjust(wspace=0.3, hspace=0.12)
for system in settings.SYSTEMS:
    rdf = pd.read_csv(f'results/rdf/{system["name"]}.csv')
    row = 0
    if settings.LINE_LABEL is not None:
        row = line_values.index(system[settings.LINE_LABEL])
    for column, name in enumerate(rdf_names):
        axes[row, column].plot(rdf.r_nm, rdf[name], color=x_colors[system[settings.X_LABEL]], lw=1.0)
for column, name in enumerate(rdf_names):
    axes[0, column].set_title(f'{name} – water O')
    axes[-1, column].set_xlabel('$r$ (nm)')
    if name in rdf_cutoff:
        for row in range(len(line_values)):
            axes[row, column].axvline(rdf_cutoff[name], color='0.3', ls='--', lw=0.8)
for row, value in enumerate(line_values):
    axes[row, 0].set_ylabel('$g(r)$')
    if value is not None:
        axes[row, -1].text(1.05, 0.5, line_name(value), transform=axes[row, -1].transAxes, ha='left', va='center')
for ax in axes.flat:
    ax.set_xlim(0.2, settings.RDF_RANGE)
    ax.set_ylim(bottom=0)
color_legend(fig, x_colors, x_title, y=0.9)
save(fig, 'rdfs')

# ----------------------------------------------------------------------
# First-shell water at each group, and its share of the first shell
# ----------------------------------------------------------------------
water = pd.read_csv('results/first_shell_water.csv')
fig, axes = plt.subplots(2, len(group_names), figsize=(2.3 * len(group_names), 3.8), sharex=True, squeeze=False)
fig.subplots_adjust(wspace=0.35, hspace=0.12)
for column, group in enumerate(group_names):
    rows = water[water.group == group]
    share = pd.read_csv(f'results/share_{group}.csv')
    for value in line_values:
        if value is None:
            part, part_share, style, label = rows, share, 'none', None
        else:
            part = rows[same_value(rows[settings.LINE_LABEL], value)]
            part_share = share[same_value(share[settings.LINE_LABEL], value)]
            style, label = '-', line_name(value)
        part = part.sort_values(settings.X_LABEL)
        part_share = part_share.sort_values(settings.X_LABEL)
        axes[0, column].errorbar(part[settings.X_LABEL], part.water_per_nm2, yerr=part.sd_blocks,
                                 color=line_colors[value], marker='o', ls=style, capsize=1.5, elinewidth=0.8,
                                 label=label)
        axes[1, column].plot(part_share[settings.X_LABEL], part_share.fraction, color=line_colors[value], marker='o',
                             ls=style)
    axes[0, column].set_title(f'first-shell water at {group}')
    axes[1, column].set_xlabel(x_title)
    axes[0, column].set_ylim(bottom=0)
    axes[1, column].set_ylim(0, 1)
axes[0, 0].set_ylabel('water (molecules/nm$^2$)')
axes[1, 0].set_ylabel('share of the first shell')
if settings.LINE_LABEL is not None:
    axes[0, 0].legend(loc='lower left', bbox_to_anchor=(0.0, 1.12), ncol=len(line_values), handlelength=1.4)
save(fig, 'first_shell_water')

# ----------------------------------------------------------------------
# MSD and P of all first-shell water against the lag time
# ----------------------------------------------------------------------
curves = pd.read_csv('results/region_curves.csv')
curves = curves[(curves.region == 'all') & curves.reliable]
fig, axes = plt.subplots(2, len(line_values), figsize=(1.6 * len(line_values) + 0.8, 3.6), sharex=True, sharey='row',
                         squeeze=False)
fig.subplots_adjust(wspace=0.12, hspace=0.12)
for system in settings.SYSTEMS:
    part = rows_of(curves, system)
    column = 0
    if settings.LINE_LABEL is not None:
        column = line_values.index(system[settings.LINE_LABEL])
    color = x_colors[system[settings.X_LABEL]]
    axes[0, column].errorbar(part.lag_ps, part.MSD_nm2, yerr=part.MSD_error_nm2, color=color, marker='o', ms=2.5,
                             capsize=1.2, elinewidth=0.7, lw=1.0)
    axes[1, column].errorbar(part.lag_ps, part.P, yerr=part.P_error, color=color, marker='o', ms=2.5, capsize=1.2,
                             elinewidth=0.7, lw=1.0)
for column, value in enumerate(line_values):
    if value is not None:
        axes[0, column].set_title(line_name(value))
    axes[1, column].set_xlabel(r'$\tau$ (ps)')
    axes[1, column].set_yscale('log')
axes[0, 0].set_ylabel('MSD (nm$^2$)')
axes[1, 0].set_ylabel(r'$P(\tau)$')
color_legend(fig, x_colors, x_title, y=0.93)
save(fig, 'msd_and_survival')

# ----------------------------------------------------------------------
# The two regimes
# ----------------------------------------------------------------------
regimes = pd.read_csv('results/clusters/regimes.csv')
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.3), gridspec_kw={'width_ratios': [2.0, 1, 1]})
fig.subplots_adjust(wspace=0.5)
# (a) the regime of every system: one cell for every system, one row for every value of LINE_LABEL
grid = np.full((len(line_values), len(x_values)), np.nan)
for system in settings.SYSTEMS:
    rows = rows_of(regimes, system)
    if len(rows) == 0:
        continue
    row = 0
    if settings.LINE_LABEL is not None:
        row = line_values.index(system[settings.LINE_LABEL])
    grid[row, x_values.index(system[settings.X_LABEL])] = rows.regime.iloc[0]
axes[0].pcolormesh(np.arange(len(x_values) + 1) - 0.5, np.arange(len(line_values) + 1) - 0.5,
                   np.ma.masked_invalid(grid), cmap=ListedColormap(REGIME_COLORS),
                   norm=BoundaryNorm([-0.5, 0.5, 1.5], 2), edgecolors='white', linewidth=0.8)
tick_names = []
for value in x_values:
    tick_names.append(value_text(value))
axes[0].set_xticks(range(len(x_values)))
axes[0].set_xticklabels(tick_names, rotation=90)
axes[0].set_xlabel(x_title)
tick_names = []
for value in line_values:
    if value is None:
        tick_names.append('')
    else:
        tick_names.append(line_name(value))
axes[0].set_yticks(range(len(line_values)))
axes[0].set_yticklabels(tick_names)
entries = []
for color, name in zip(REGIME_COLORS, REGIME_NAMES):
    entries.append(Patch(facecolor=color, label=name))
axes[0].legend(handles=entries, loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2, handlelength=1.1)
# (b), (c) the curves of the two regimes: mean and standard deviation
panels = [(axes[1], 'MSD_rel', 'MSD relative to the first lag'), (axes[2], 'P_rel', '$P$ relative to the first lag')]
for ax, name, title in panels:
    data = pd.read_csv(f'results/clusters/regime_mean_{name}.csv')
    for regime, color in enumerate(REGIME_COLORS):
        rows = data[data.regime == regime]
        ax.fill_between(rows.lag_ps, rows['mean'] - rows.sd, rows['mean'] + rows.sd, color=color, alpha=0.2, lw=0)
        ax.plot(rows.lag_ps, rows['mean'], color=color, lw=1.6)
    ax.set_xlabel(r'$\tau$ (ps)')
    ax.set_title(title)
axes[2].set_yscale('log')
axes[2].yaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))
axes[2].yaxis.set_major_formatter(FormatStrFormatter('%g'))
axes[2].yaxis.set_minor_formatter(NullFormatter())
for ax, letter in zip(axes, 'abc'):
    ax.text(-0.1, 1.2, f'({letter})', transform=ax.transAxes, fontsize=9.0, va='bottom', ha='right')
save(fig, 'regimes')

# ----------------------------------------------------------------------
# Tests of the clusters
# ----------------------------------------------------------------------
scores = pd.read_csv('results/cluster_tests/k_scores.csv')
renewal = pd.read_csv('results/cluster_tests/renewal_index.csv')
exponent = pd.read_csv('results/cluster_tests/msd_exponent.csv')
fig, axes = plt.subplots(1, 2, figsize=(5.4, 2.3))
fig.subplots_adjust(wspace=0.45)
axes[0].plot(scores.k, scores.silhouette, color='k', marker='o', label='silhouette score')
axes[0].plot(scores.k, scores.bootstrap_stability, color='k', marker='s', ls='--', label='bootstrap stability')
axes[0].set_xticks(scores.k)
axes[0].set_xlabel('number of clusters $k$')
axes[0].set_ylim(0, 1.05)
axes[0].legend(loc='lower left')
for regime, (color, name) in enumerate(zip(REGIME_COLORS, REGIME_NAMES)):
    keep = renewal.regime == regime
    axes[1].scatter(exponent.alpha[keep], renewal.P_ratio[keep], s=14, color=color, edgecolor='0.2', linewidth=0.3,
                    label=name)
axes[1].set_xlabel(r'MSD exponent $\alpha$')
axes[1].set_ylabel('renewal index')
axes[1].legend(loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2, handletextpad=0.2)
for ax, letter in zip(axes, 'ab'):
    ax.text(-0.25, 1.12, f'({letter})', transform=ax.transAxes, fontsize=9.0, va='bottom', ha='right')
save(fig, 'cluster_tests')
