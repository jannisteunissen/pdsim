#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Plot avalanche statistics')
parser.add_argument('avalanche_file', type=str,
                    help='Input file with avalanche sizes')
parser.add_argument('-bins', type=int, default=20, help='Number of bins')
parser.add_argument('-range', type=float, nargs=2,
                    help='Range for avalanche size')
parser.add_argument('-alpha', type=float,
                    help='Townsend ionization coefficient (1/m)')
parser.add_argument('-eta', type=float,
                    help='Townsend attachment coefficient (1/m)')
parser.add_argument('-L', type=float,
                    help='Length of domain')
args = parser.parse_args()

stats = np.loadtxt(args.avalanche_file, skiprows=1)

r0 = stats[:, 0:3]
avalanche_sizes = stats[:, 3]
mean_size = avalanche_sizes.mean()
n_avalanches = len(avalanche_sizes)

fig, ax = plt.subplots()

counts, bins, patches = ax.hist(avalanche_sizes, bins=args.bins,
                                range=args.range, density=True)
ax.set_xlabel('Avalanche size')
ax.set_ylabel('Count')
ax.vlines(mean_size, 0., counts.max(), color='black',
          label=f'mean = {mean_size:.03e}')

if args.alpha is not None:
    K = (args.alpha - args.eta) * args.L
    expected_size = np.exp(K)

    ax.vlines(expected_size, 0., counts.max(), color='gray',
              label=f'exp(K) = {expected_size:.03e}')

    x = np.arange(0.0, bins.max(), 1.0)
    dx = x[1] - x[0]
    y = 1/expected_size * (1 - 1/expected_size)**(x-1)

    ax.plot(x, y, label='eq (7)')

    n_e = expected_size / (1 - args.eta/args.alpha)
    y = (1 - args.eta/args.alpha) * np.exp(-x/n_e)/n_e
    y[0] = args.eta/args.alpha
    ax.plot(x, y, label='eq (26)')

ax.legend()
plt.show()
