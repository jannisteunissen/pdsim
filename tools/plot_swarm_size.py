#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import geom
import argparse
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Plot avalanche statistics')
parser.add_argument('avalanche_file', type=str,
                    help='Input file with avalanche sizes')
parser.add_argument('-binwidth', type=int,
                    help='Width of each bin')
parser.add_argument('-range', type=float, nargs=2,
                    help='Range for avalanche size')
parser.add_argument('-quantile', type=float, default=0.99,
                    help='Take this quantile as upper limit')
parser.add_argument('-savefig', type=str,
                    help='Save figure to this file')
parser.add_argument('-text', type=str,
                    help='Add this text to the figure')
args = parser.parse_args()

p_m1, K_star, x, y, z = np.loadtxt(args.avalanche_file, max_rows=1)
avalanche_sizes = np.loadtxt(args.avalanche_file, skiprows=1, dtype=int)
n_avalanches = len(avalanche_sizes)
mean_size = avalanche_sizes.mean()
mean_size_std = avalanche_sizes.std() / np.sqrt(n_avalanches-1)

# Remove size 1 avalanches
mask = (avalanche_sizes == 1)
p_size_one = mask.sum() / len(avalanche_sizes)
p_size_one_std = np.sqrt(p_size_one * (1 - p_size_one) / (n_avalanches-1))
avalanche_sizes = avalanche_sizes[~mask]

if len(avalanche_sizes) == 0:
    raise ValueError('All avalanches have size one')

min_size = 2
max_size = int(np.ceil(np.quantile(avalanche_sizes, 0.99)))

if args.binwidth == None:
    # Automatically determine bin width
    args.binwidth = max(1, max_size//30)

bins = np.arange(min_size, max_size+1, args.binwidth)
bins = bins - 0.5

cm = 1/2.54
fig, ax = plt.subplots(figsize=(13*cm, 8*cm), layout='constrained')

counts, bins, _ = ax.hist(avalanche_sizes, bins=bins,
                          density=True)

ax.set_xlabel('M (number of ionizations)')
ax.set_ylabel('P(X = M)')

# Enforce scientific notation with scale factor and fixed precision
class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self):  # Override function that finds format to use.
        self.format = "%2.2f"  # Give format here

yfmt = ScalarFormatterForceFormat()
yfmt.set_powerlimits((0,0))
ax.yaxis.set_major_formatter(yfmt)

p_geom = (1 - p_m1)/(np.exp(K_star) - 1)
expected_size = np.exp(K_star)

x = np.arange(2, bins.max())
y = np.zeros(len(x))

ax.text(0.5, 0.5, r'$\bar{M}_\mathrm{sim} =$' +
        f'{mean_size:.2e} ({mean_size_std:.1e})\n' +
        r'$\bar{M} =$' + f'{expected_size:.2e}\n' +
        r"$P'_{1, \mathrm{sim}}$ = " +
        f'{p_size_one:.1e} ({p_size_one_std:.1e})' + '\n' +
        f"$P'_1$ = " + f'{p_m1:.1e}\n',
        transform=ax.transAxes, verticalalignment='center',
        horizontalalignment='left')

if args.text is not None:
    ax.text(0.5, 1.02, args.text, transform=ax.transAxes,
            verticalalignment='bottom', horizontalalignment='center')

# The avalanche size is given by m+1, with a geometric distribution for m
y = geom.pmf(x-1, p_geom)

ax.plot(x, y, label='geom. pmf')
ax.legend()

if args.savefig is not None:
    rc('font', size=10)
    rc('font', family='serif')
    rc('legend', fontsize=10)
    plt.savefig(args.savefig, dpi=200, bbox_inches='tight', pad_inches=0.01)
    print(f'Saved {args.savefig}')
else:
    plt.show()
