"""
The figure of the tests of the two regimes of all first-shell water, from the result files of step 1.

    figures/water_regime_tests.pdf: (a) silhouette score and bootstrap stability for k = 2, 3 and 4 groups (shaded:
        k = 2); (b) renewal index P(100 ps)/P(20 ps) and (c) exponent alpha of MSD ~ tau^alpha (20 to 100 ps) of every
        system over the chain lengths and coverages; white lines: border between the two regimes

The figure is also written as PNG.
Run:  python3 plot_figures.py      (after step 1)
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')                   # write files only, no window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# Fonts, lines and ticks
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'pdf.fonttype': 42, 'mathtext.default': 'it',
    'axes.unicode_minus': False, 'axes.labelsize': 8.0, 'xtick.labelsize': 8.0, 'ytick.labelsize': 8.0,
    'legend.fontsize': 7.5, 'axes.linewidth': 1.2, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.labelpad': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.major.size': 4.2,
    'ytick.major.size': 4.2, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 'xtick.minor.size': 2.0,
    'ytick.minor.size': 2.0, 'xtick.minor.width': 0.7, 'ytick.minor.width': 0.7, 'lines.linewidth': 1.45,
    'lines.markersize': 4.5, 'lines.markeredgewidth': 0.7, 'legend.frameon': False, 'legend.handlelength': 1.6,
    'legend.handletextpad': 0.4, 'legend.labelspacing': 0.28})


def draw_map(ax, data, column):
    """One value per system over the chain lengths and coverages, with the border between the two regimes drawn in
    white along the cell edges."""
    values = data.pivot(index='L', columns='coverage', values=column)
    regimes = data.pivot(index='L', columns='coverage', values='regime').to_numpy(float)
    n_rows, n_columns = values.shape
    image = ax.pcolormesh(np.arange(n_columns + 1) - 0.5, np.arange(n_rows + 1) - 0.5, values.to_numpy(float),
                          cmap='viridis', shading='flat', edgecolors='none', antialiased=False)
    for i in range(n_rows):
        for j in range(n_columns - 1):
            if regimes[i, j] != regimes[i, j + 1]:
                ax.plot([j + 0.5, j + 0.5], [i - 0.5, i + 0.5], color='white', solid_capstyle='butt', zorder=5)
    for i in range(n_rows - 1):
        for j in range(n_columns):
            if regimes[i, j] != regimes[i + 1, j]:
                ax.plot([j - 0.5, j + 0.5], [i + 0.5, i + 0.5], color='white', solid_capstyle='butt', zorder=5)
    ax.set_xlim(-0.5, n_columns - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.set_xticks(range(n_columns))
    ax.set_xticklabels([f'{coverage:.2f}' for coverage in values.columns], rotation=90)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([f'$L$ = {L}' for L in values.index])
    ax.set_xlabel(r'$\Gamma$ (molecules/nm$^2$)')
    ax.tick_params(pad=1.5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.9)
    return image


Path('figures').mkdir(exist_ok=True)
fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.1), gridspec_kw={'width_ratios': [1, 1.45, 1.45]})
fig.subplots_adjust(wspace=0.85)

# (a) silhouette score and bootstrap stability against the number of groups; k = 2 shaded
ax = axes[0]
choice = pd.read_csv('results/k_scores.csv')
ax.axvspan(1.75, 2.25, color='#E6E6E6', lw=0, zorder=0)
ax.plot(choice.k, choice.silhouette, color='k', marker='o', zorder=3, label='silhouette score')
ax.plot(choice.k, choice.bootstrap_stability, color='#7A7A7A', marker='s', mfc='white', mec='#7A7A7A', zorder=3,
        label='bootstrap stability')
ax.set_xlim(choice.k.min() - 0.4, choice.k.max() + 0.4)
ax.set_xticks(choice.k)
ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('number of groups $k$')
ax.set_ylabel('score')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0))

# (b) renewal index and (c) exponent of every system, with the border between the regimes
for ax, file, column, label in [(axes[1], 'renewal_index', 'P_ratio',
                                 r'$\mathrm{P}(100\,\mathrm{ps})/\mathrm{P}(20\,\mathrm{ps})$'),
                                (axes[2], 'msd_exponent', 'alpha', r'exponent $\alpha$')]:
    image = draw_map(ax, pd.read_csv(f'results/{file}.csv'), column)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.06, pad=0.03)
    colorbar.set_label(label)
for ax, letter in zip(axes, 'abc'):
    ax.text(-0.3, 1.2, f'({letter})', transform=ax.transAxes, fontsize=9.0, va='bottom')

fig.savefig('figures/water_regime_tests.pdf', bbox_inches='tight')
fig.savefig('figures/water_regime_tests.png', bbox_inches='tight', dpi=200)
plt.close(fig)
