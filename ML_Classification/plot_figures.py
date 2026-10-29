"""
The figure of the two regimes of the first-shell water dynamics, from the result files of step 1.

    figures/water_dynamics_regimes.pdf: (a) the regime of every system over the chain lengths and coverages; (b), (c)
        survival-conditioned lateral MSD and contact survival P(tau) of the first-shell water, each divided by its
        value at 20 ps: mean over the systems of each regime (line) and standard deviation (band)

The figure is also written as PNG.
Run:  python3 plot_figures.py      (after step 1)
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')                   # write files only, no window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter, LogLocator, NullFormatter

REGIME_COLORS = ['#2BA58A', '#1F5F99']
REGIME_NAMES = ['faster renewal', 'slower renewal']

# Fonts, lines and ticks
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'pdf.fonttype': 42, 'mathtext.default': 'it',
    'axes.unicode_minus': False, 'axes.labelsize': 8.0, 'xtick.labelsize': 8.0, 'ytick.labelsize': 8.0,
    'legend.fontsize': 8.0, 'legend.title_fontsize': 8.0, 'axes.linewidth': 1.2, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.labelpad': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 4.2, 'ytick.major.size': 4.2, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
    'xtick.minor.size': 2.0, 'ytick.minor.size': 2.0, 'xtick.minor.width': 0.7, 'ytick.minor.width': 0.7,
    'legend.frameon': False, 'legend.handletextpad': 0.4, 'legend.columnspacing': 1.2})


def draw_regime_map(ax, regimes, chain_lengths, coverages):
    """The regime of every system as colored cells, with white lines along the cell edges between neighboring systems
    of different regimes."""
    n_rows, n_columns = regimes.shape
    ax.pcolormesh(np.arange(n_columns + 1) - 0.5, np.arange(n_rows + 1) - 0.5, regimes,
                  cmap=ListedColormap(REGIME_COLORS), norm=BoundaryNorm(np.arange(-0.5, 2), 2), shading='flat',
                  edgecolors='none', antialiased=False)
    for i in range(n_rows):
        for j in range(n_columns - 1):
            if regimes[i, j] != regimes[i, j + 1]:
                ax.plot([j + 0.5, j + 0.5], [i - 0.5, i + 0.5], color='white', lw=1.45, solid_capstyle='butt')
    for i in range(n_rows - 1):
        for j in range(n_columns):
            if regimes[i, j] != regimes[i + 1, j]:
                ax.plot([j - 0.5, j + 0.5], [i + 0.5, i + 0.5], color='white', lw=1.45, solid_capstyle='butt')
    ax.set_xlim(-0.5, n_columns - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.set_xticks(range(n_columns))
    ax.set_xticklabels(coverages, rotation=90)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([f'$L$ = {L}' for L in chain_lengths])
    ax.set_xlabel(r'$\Gamma$ (molecules/nm$^2$)')
    ax.tick_params(pad=1.5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.9)


Path('figures').mkdir(exist_ok=True)
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.3), gridspec_kw={'width_ratios': [2.3, 1, 1]})
fig.subplots_adjust(wspace=0.55)

# (a) map of the regimes
regime_map = pd.read_csv('results/regime_map.csv', index_col='L')
draw_regime_map(axes[0], regime_map.to_numpy(), regime_map.index, regime_map.columns)
entries = []
for color, name in zip(REGIME_COLORS, REGIME_NAMES):
    entries.append(Rectangle((0, 0), 1, 1, fc=color, ec='k', lw=0.4, label=name))
axes[0].legend(handles=entries, loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2, handlelength=1.1,
               borderaxespad=0.2, title='water-dynamics regimes')

# (b), (c) the curves of the two regimes: mean and standard deviation
for ax, name in ((axes[1], 'MSD_rel'), (axes[2], 'P_rel')):
    data = pd.read_csv(f'results/regime_mean_{name}.csv')
    for regime, color in enumerate(REGIME_COLORS):
        rows = data[data.regime == regime]
        ax.fill_between(rows.lag_ps, rows['mean'] - rows.sd, rows['mean'] + rows.sd, color=color, alpha=0.18, lw=0)
        ax.plot(rows.lag_ps, rows['mean'], color=color, lw=1.6)
    ax.set_xlim(15, 105)
    ax.set_xticks([20, 60, 100])
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.set_xlabel(r'$\tau$ (ps)')
    ax.set_title('first-shell water', loc='right', fontsize=8.0, pad=2)
axes[1].set_ylabel(r'MSD$^{2D}_{\rm surv}(\tau)$/MSD$^{2D}_{\rm surv}(20\,{\rm ps})$')
axes[1].yaxis.set_minor_locator(AutoMinorLocator())
axes[2].set_ylabel(r'$\mathrm{P}(\tau)/\mathrm{P}(20\,\mathrm{ps})$')
axes[2].set_yscale('log')                # P decays about exponentially
axes[2].set_ylim(0.1, 1.15)
axes[2].set_yticks([0.1, 0.2, 0.5, 1.0])
axes[2].yaxis.set_major_formatter(FormatStrFormatter('%g'))
axes[2].yaxis.set_minor_formatter(NullFormatter())
axes[2].yaxis.set_minor_locator(LogLocator(subs=np.arange(2, 10)))
for ax, letter in zip(axes, 'abc'):
    ax.text(-0.12 if letter == 'a' else -0.42, 1.25, f'({letter})', transform=ax.transAxes, fontsize=9.0, va='bottom')

fig.savefig('figures/water_dynamics_regimes.pdf', bbox_inches='tight')
fig.savefig('figures/water_dynamics_regimes.png', bbox_inches='tight', dpi=200)
plt.close(fig)
