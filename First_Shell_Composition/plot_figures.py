"""
The figure of the composition of the first hydration shell, from the result files of step 1.

    figures/first_shell_composition.pdf: (a) the group closest to most of the first-shell water molecules (color),
        shaded by the share of all first-shell water molecules of the system closest to it (below 1/2, 1/2 to 2/3,
        above 2/3), interpolated between the systems with a shape-preserving cubic interpolation (PCHIP), first along
        the coverages and then along the chain lengths; thick lines where the largest group changes. (b)-(d) Fraction
        of all first-shell water molecules of a system that are closest to CH3, CH2 and the oxygens, one line per
        chain length.

The figure is also written as PNG.
Run:  python3 plot_figures.py      (after step 1)
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')                   # write files only, no window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from scipy.interpolate import PchipInterpolator

GROUPS = ['CH3', 'CH2', 'oxygen']
GROUP_NAMES = {'CH3': r'CH$_3$', 'CH2': r'CH$_2$', 'oxygen': 'silica and headgroup O'}
GROUP_COLORS = {'CH3': '#7B4FA0', 'CH2': '#5F6B73', 'oxygen': '#D2473A'}
LENGTH_COLORS = {4: '#4DBBD5', 5: '#E68613', 6: '#756BB1', 7: '#B7371B', 8: '#0D4F8B'}
LEVELS = [0.5, 2 / 3]                   # the largest share below 1/2, from 1/2 to 2/3, above 2/3
SHADES = [0.22, 0.4, 0.6]               # the share of the group color in the color of each level (the rest white)
COVERAGE_TITLE = r'$\Gamma$ (molecules/nm$^2$)'

# Fonts, lines and ticks
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'pdf.fonttype': 42, 'mathtext.default': 'it',
    'axes.unicode_minus': False, 'axes.labelsize': 8.0, 'xtick.labelsize': 8.0, 'ytick.labelsize': 8.0,
    'legend.fontsize': 7.5, 'axes.linewidth': 1.2, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.labelpad': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.major.size': 4.2,
    'ytick.major.size': 4.2, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 'xtick.minor.size': 2.0,
    'ytick.minor.size': 2.0, 'xtick.minor.width': 0.7, 'ytick.minor.width': 0.7, 'lines.markeredgewidth': 0.7,
    'lines.dashed_pattern': [4.0, 2.2], 'lines.scale_dashes': False, 'legend.frameon': False,
    'legend.handlelength': 1.8, 'legend.handletextpad': 0.4, 'legend.labelspacing': 0.28, 'legend.columnspacing': 1.0})


def shade_of(color, amount):
    """The color mixed with white: amount 0 is white, 1 the color."""
    rgb = to_rgb(color)
    return (1 - amount * (1 - rgb[0]), 1 - amount * (1 - rgb[1]), 1 - amount * (1 - rgb[2]))


def length_style(L, k):
    """Color and marker of a chain length: its own color (the k-th tab color for another L); circles for even L,
    squares for odd L."""
    color = LENGTH_COLORS.get(L, plt.get_cmap('tab10').colors[k % 10])
    if L % 2 == 0:
        return color, 'o'
    return color, 's'


Path('figures').mkdir(exist_ok=True)
shares = {}
for group in GROUPS:
    shares[group] = pd.read_csv(f'results/share_{group}.csv').pivot(index='L', columns='coverage', values='fraction')
chain_lengths = list(shares['CH3'].index)
coverages = list(shares['CH3'].columns)

fig = plt.figure(figsize=(7.0, 4.3))
grid = GridSpec(2, 6, figure=fig, height_ratios=[1, 1], hspace=0.95, wspace=1.2)

# (a) the group closest to most first-shell water molecules: the shares between the systems by PCHIP in the coverage
# and then in the chain length; every point in the color of its largest share, in the shade of its level, with thin
# lines between the levels and thick lines where the largest group changes
ax = fig.add_subplot(grid[0, 1:5])
x = np.linspace(-0.5, len(coverages) - 0.5, 400)       # map position; column j is the coverage coverages[j]
y = np.linspace(-0.5, len(chain_lengths) - 0.5, 200)
coverage_of_x = np.interp(x, np.arange(len(coverages)), coverages)
surfaces = []
for group in GROUPS:
    along_coverages = PchipInterpolator(coverages, shares[group].to_numpy(), axis=1)(coverage_of_x)
    surfaces.append(PchipInterpolator(np.arange(len(chain_lengths)), along_coverages, axis=0)(
        np.clip(y, 0, len(chain_lengths) - 1)))
surfaces = np.array(surfaces)           # group, y, x
first = surfaces.argmax(axis=0)         # the group of the largest share
largest_share = surfaces.max(axis=0)
level = np.digitize(largest_share, LEVELS)                          # 0: below 1/2, 1: 1/2 to 2/3, 2: above 2/3
for g, group in enumerate(GROUPS):          # every group and level as a filled area
    for n, shade in enumerate(SHADES):
        mask = ((first == g) & (level == n)).astype(float)
        if mask.any():
            ax.contourf(x, y, mask, levels=[0.5, 1.5], colors=[shade_of(GROUP_COLORS[group], shade)])
ax.contour(x, y, largest_share, levels=LEVELS, colors='#34404c', linewidths=0.8)
for a, b in [(0, 1), (1, 2), (0, 2)]:
    other = 3 - a - b
    difference = np.where(surfaces[other] < np.maximum(surfaces[a], surfaces[b]), surfaces[a] - surfaces[b], np.nan)
    ax.contour(x, y, difference, levels=[0], colors='#34404c', linewidths=2.0)
ax.set_xlim(-0.5, len(coverages) - 0.5)
ax.set_ylim(-0.5, len(chain_lengths) - 0.5)
ax.set_xticks(range(len(coverages)))
ax.set_xticklabels([f'{coverage:.2f}' for coverage in coverages], rotation=45, ha='right', rotation_mode='anchor')
ax.set_yticks(range(len(chain_lengths)))
ax.set_yticklabels([f'$L$ = {L}' for L in chain_lengths])
ax.set_xlabel(COVERAGE_TITLE)
ax.tick_params(pad=2)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
for side in ('top', 'right'):
    ax.spines[side].set_visible(True)
    ax.spines[side].set_linewidth(0.9)
# the legend above the map: the groups (in their middle shade), then the three shades
entries = []
for group in GROUPS:
    entries.append(Rectangle((0, 0), 1, 1, fc=shade_of(GROUP_COLORS[group], SHADES[1]), ec='0.5', lw=0.4,
                             label=GROUP_NAMES[group]))
legend = ax.legend(handles=entries, loc='lower left', bbox_to_anchor=(-0.02, 1.13), ncol=3, handlelength=1.3,
                   borderaxespad=0.15)
ax.add_artist(legend)
entries = []
for shade, label in zip(SHADES, ['share < 1/2', '1/2 to 2/3', '> 2/3']):
    entries.append(Rectangle((0, 0), 1, 1, fc=shade_of('#606060', shade), ec='0.5', lw=0.4, label=label))
ax.legend(handles=entries, loc='lower left', bbox_to_anchor=(-0.02, 0.99), ncol=3, handlelength=1.3,
          columnspacing=0.8, borderaxespad=0.15)
ax.text(-0.22, 1.3, '(a)', transform=ax.transAxes, fontsize=9.0, va='bottom')

# (b)-(d) share of the first-shell water at every group, one line per chain length
for n, (letter, group) in enumerate(zip('bcd', GROUPS)):
    ax = fig.add_subplot(grid[1, 2 * n:2 * n + 2])
    table = shares[group]
    for k, L in enumerate(chain_lengths):
        color, marker = length_style(L, k)
        ax.plot(table.columns, table.loc[L], color=color, ls='--', lw=1.2, marker=marker, ms=3.8, mec=color)
    ax.set_xlabel(COVERAGE_TITLE)
    ax.set_ylim(0, 0.92)
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    if n == 0:
        ax.set_ylabel('fraction of all\nfirst-shell water\nclosest to the group')
    if group == 'oxygen':
        ax.set_title('closest to the oxygens', fontsize=8.0, pad=3)
    else:
        ax.set_title(f'closest to {GROUP_NAMES[group]}', fontsize=8.0, pad=3)
    ax.text(-0.3 if n == 0 else -0.15, 1.08, f'({letter})', transform=ax.transAxes, fontsize=9.0, ha='right',
            va='bottom')
entries = []
for k, L in enumerate(chain_lengths):
    color, marker = length_style(L, k)
    entries.append(Line2D([], [], color=color, ls='--', lw=1.2, marker=marker, ms=3.8, label=f'$L$ = {L}'))
fig.legend(handles=entries, loc='center', bbox_to_anchor=(0.5, 0.485), ncol=len(entries))

fig.savefig('figures/first_shell_composition.pdf', bbox_inches='tight')
fig.savefig('figures/first_shell_composition.png', bbox_inches='tight', dpi=200)
plt.close(fig)
