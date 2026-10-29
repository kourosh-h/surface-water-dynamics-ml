"""
Figures of the MSD and the contact survival P of all first-shell water, from results/region_curves.csv (step 2).

    figures/first_shell_water_msd.pdf        survival-conditioned lateral MSD against the lag time, one panel per
                                             coverage and one line per chain length
    figures/first_shell_water_survival.pdf   contact-survival probability P(tau) against the lag time, on a
                                             logarithmic axis

Bands: standard deviation of the ten 9-ns block values. Every figure is also written as PNG.
Run:  python3 plot_figures.py      (after step 2)
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')                   # write files only, no window
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator

REGION = 'all'                          # all first-shell water
PANELS_PER_ROW = 5
LENGTH_COLORS = {4: '#4DBBD5', 5: '#E68613', 6: '#756BB1', 7: '#B7371B', 8: '#0D4F8B'}

# Fonts, lines and ticks
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'pdf.fonttype': 42, 'mathtext.default': 'it',
    'axes.unicode_minus': False, 'axes.labelsize': 8.0, 'xtick.labelsize': 8.0, 'ytick.labelsize': 8.0,
    'legend.fontsize': 8.0, 'legend.title_fontsize': 8.0, 'axes.linewidth': 1.2, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.labelpad': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 4.2, 'ytick.major.size': 4.2, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
    'xtick.minor.size': 2.0, 'ytick.minor.size': 2.0, 'xtick.minor.width': 0.7, 'ytick.minor.width': 0.7,
    'lines.linewidth': 1.45, 'lines.markeredgewidth': 0.7, 'legend.frameon': False, 'legend.handlelength': 1.8,
    'legend.handletextpad': 0.4, 'legend.columnspacing': 1.0, 'legend.labelspacing': 0.28})


def length_style(L, k):
    """Color and marker of a chain length: its own color (the k-th tab color for another L); circles for even L,
    squares for odd L."""
    color = LENGTH_COLORS.get(L, plt.get_cmap('tab10').colors[k % 10])
    if L % 2 == 0:
        return color, 'o'
    return color, 's'


def save(fig, name):
    """The figure as PDF and PNG in figures/."""
    fig.savefig(f'figures/{name}.pdf', bbox_inches='tight')
    fig.savefig(f'figures/{name}.png', bbox_inches='tight', dpi=200)
    plt.close(fig)


def curves_figure(curves, column, error, label, name, log_scale):
    """One panel per coverage, one line per chain length with its band; the panels share the axes."""
    coverages = sorted(curves.coverage.unique())
    chain_lengths = sorted(curves.L.unique())
    n_rows = (len(coverages) + PANELS_PER_ROW - 1) // PANELS_PER_ROW
    n_columns = min(len(coverages), PANELS_PER_ROW)
    fig, axes = plt.subplots(n_rows, n_columns, figsize=(1.4 * n_columns + 0.3, 1.45 * n_rows + 0.3), sharex=True,
                             sharey=True, squeeze=False)
    fig.subplots_adjust(wspace=0.12, hspace=0.35)
    for n, coverage in enumerate(coverages):
        ax = axes[n // PANELS_PER_ROW, n % PANELS_PER_ROW]
        for k, L in enumerate(chain_lengths):
            color, marker = length_style(L, k)
            rows = curves[(curves.coverage == coverage) & (curves.L == L)]
            ax.fill_between(rows.lag_ps, rows[column] - rows[error], rows[column] + rows[error], color=color,
                            alpha=0.18, linewidth=0, zorder=2)
            ax.plot(rows.lag_ps, rows[column], color=color, marker=marker, ms=2.4, zorder=3)
        ax.set_title(rf'$\Gamma$ = {coverage:.2f}', fontsize=8.0, pad=3)
    for ax in axes.flat[len(coverages):]:           # no data for the last places of the last row
        ax.axis('off')
    for ax in axes.flat:
        ax.set_xlim(0, 210)
        ax.set_xticks([0, 100, 200])
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        if log_scale:
            ax.set_yscale('log')
        else:
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    for row in range(n_rows):
        axes[row, 0].set_ylabel(label)
    for column_index in range(n_columns):
        axes[-1, column_index].set_xlabel(r'$\tau$ (ps)')
    entries = []
    for k, L in enumerate(chain_lengths):
        color, marker = length_style(L, k)
        entries.append(Line2D([], [], color=color, marker=marker, ms=4, lw=1.2, label=f'$L$ = {L}'))
    entries.append(Line2D([], [], color='none', label=r'$\Gamma$ in molecules/nm$^2$'))
    fig.legend(handles=entries, loc='lower center', bbox_to_anchor=(0.5, 0.92), ncol=len(entries))
    save(fig, name)


Path('figures').mkdir(exist_ok=True)
curves = pd.read_csv('results/region_curves.csv')
curves = curves[curves.region == REGION]
curves_figure(curves, 'MSD_nm2', 'MSD_error_nm2', r'MSD$^{2D}_{\rm surv}(\tau)$ (nm$^2$)', 'first_shell_water_msd',
              False)
curves_figure(curves, 'P', 'P_error', r'$\mathrm{P}(\tau)$', 'first_shell_water_survival', True)
