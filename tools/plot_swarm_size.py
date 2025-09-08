#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import geom
import argparse


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Plot avalanche statistics')
parser.add_argument('avalanche_file', type=str,
                    help='Input file with avalanche sizes')
parser.add_argument('-binwidth', type=int, default=1,
                    help='Width of each bin')
parser.add_argument('-range', type=float, nargs=2,
                    help='Range for avalanche size')
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
bins = np.arange(min_size, max_size+1, args.binwidth)
bins = bins - 0.5

fig, ax = plt.subplots()

counts, bins, _ = ax.hist(avalanche_sizes, bins=bins,
                          density=True)

ax.set_xlabel('Number of ionizations')
ax.set_ylabel('Count')
ax.vlines(mean_size, 0., counts.max(), color='black',
          label=f'mean = {mean_size:.03e}')

expected_size = np.exp(K_star)

ax.vlines(expected_size, 0., counts.max(), color='gray',
          label=f'exp(K*) = {expected_size:.03e}')

x = np.arange(1, bins.max())
y = np.zeros(len(x))

p_geom = (1 - p_m1)/(np.exp(K_star) - 1)

# Probability of no further ionization (only the initial one)
y[0] = p_m1

# The avalanche size is given by m+1, with a geometric distribution for m
y[1:] = (1 - p_m1) * geom.pmf(x[1:]-1, p_geom)

ax.plot(x[1:], y[1:], label='Approximation')

print(f'Expected Pr(X=1): {p_m1:.2e}, simulation: {p_size_one:.2e}')

ax.legend()
plt.show()
