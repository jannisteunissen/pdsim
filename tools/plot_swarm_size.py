#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import geom
import argparse
from matplotlib import rc

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Plot avalanche statistics')
parser.add_argument('avalanche_file', type=str,
                    help='Input file with avalanche sizes')
parser.add_argument('-binwidth', type=int,
                    help='Width of each bin')
parser.add_argument('-range', type=float, nargs=2,
                    help='Range for avalanche size')
parser.add_argument('-savefig', type=str,
                    help='Save figure to')
args = parser.parse_args()

p_m1, K_star, x, y, z = np.loadtxt(args.avalanche_file, max_rows=1)
avalanche_sizes = np.loadtxt(args.avalanche_file, skiprows=1, dtype=int)
mean_size = avalanche_sizes.mean()

# Remove size 1 avalanches
mask = (avalanche_sizes == 1)
p_size_one = mask.sum() / len(avalanche_sizes)
avalanche_sizes = avalanche_sizes[~mask]

if len(avalanche_sizes) == 0:
    raise ValueError('All avalanches have size one')

min_size = 2
max_size = avalanche_sizes.max()

if args.binwidth == None:
    # Automatically determine bin width
    args.binwidth = max(1, max_size//30)

bins = np.arange(min_size, max_size+1, args.binwidth)
bins = bins - 0.5

cm = 1/2.54
fig, ax = plt.subplots(figsize=(13*cm, 9*cm), layout='constrained')

counts, bins, _ = ax.hist(avalanche_sizes, bins=bins,
                          density=True)

ax.set_xlabel('M (number of ionizations)')
ax.set_ylabel('P(X = M)')

expected_size = np.exp(K_star)
ax.text(0.7, 0.5, r'$\bar{M}_\mathrm{sim} =$' + f'{mean_size:.02e}\n' +
        r'$\bar{M} =$' + f'{expected_size:.02e}\n' +
        r"$P'_{1, \mathrm{sim}}$ = " + f'{p_size_one:.3f}' + '\n' +
        f"$P'_1$ = " + f'{p_m1:.3f}\n',
        transform=ax.transAxes, verticalalignment='center',
        horizontalalignment='center')

x = np.arange(1, bins.max())
y = np.zeros(len(x))

p_geom = (1 - p_m1)/(np.exp(K_star) - 1)

# Probability of no further ionization (only the initial one)
y[0] = p_m1

# The avalanche size is given by m+1, with a geometric distribution for m
y[1:] = (1 - p_m1) * geom.pmf(x[1:]-1, p_geom)

ax.plot(x[1:], y[1:], label='geometric')
ax.legend()

if args.savefig is not None:
    rc('font', size=10)
    rc('font', family='serif')
    rc('legend', fontsize=10)
    plt.savefig(args.savefig, dpi=200, bbox_inches='tight', pad_inches=0.01)
    print(f'Saved {args.savefig}')
else:
    plt.show()
