#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import subprocess

import argparse

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Utility to find the approximate inception voltage')
p.add_argument('cfg', type=str, help='Config file')
p.add_argument('-p', type=float, default=1e-3,
               help='Threshold inception probability')
p.add_argument('-fbound', type=float, nargs=2, required=True,
               help='Bounds for the field scale factor')
args = p.parse_args()


def get_sample(factor):
    proc = subprocess.run(['./pdsim', args.cfg,
                           r'-output%verbosity=0',
                           r'-input%field_scale_factor=' + f'{factor}'],
                          capture_output=True, check=True)
    prob, _, _ = map(float, proc.stdout.split())
    return prob - args.p


# Assume a < b, fa = f(a) < 0, fb = f(b) > 0.
# Taken from: https://www.johndcook.com/blog/2012/06/14/root-finding-with-noisy-functions/
def noisy_bisect(f, a, b, fa, fb, tolerance):
    if b-a < tolerance:
        return (a, b)

    mid = 0.5 * (a+b)
    fmid = f(mid)

    if fmid < fa or fmid > fb:
        # Monotonicity violated.
        # Reached resolution of noise.
        return (a, b)
    if fmid < 0:
        a, fa = mid, fmid
    else:
        b, fb = mid, fmid
    return noisy_bisect(f, a, b, fa, fb, tolerance)


fa = get_sample(args.fbound[0])
fb = get_sample(args.fbound[1])
bracket = noisy_bisect(get_sample, args.fbound[0], args.fbound[1],
                       fa, fb, 1e-3)

print(f'field_scale_factor is about: {0.5 * sum(bracket):.3e}')
