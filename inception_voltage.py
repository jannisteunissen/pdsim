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
p.add_argument('-refine_steps', type=int, default=0,
               help='Refine for this number of steps')
args = p.parse_args()


# Function to find a root of
def target_function(factor):
    proc = subprocess.run(['./pdsim', args.cfg,
                           r'-output%verbosity=0',
                           r'-input%field_scale_factor=' + f'{factor}'],
                          capture_output=True, check=True)
    p_inception, _, _ = map(float, proc.stdout.split())
    return p_inception - args.p


# Assume a < b, fa = f(a) < 0, fb = f(b) > 0.
# Taken from: https://www.johndcook.com/blog/2012/06/14/root-finding-with-noisy-functions/
def noisy_bisect(f, a, b, fa, fb, tolerance):
    if b-a < tolerance:
        return (a, b)

    mid = 0.5 * (a+b)
    fmid = f(mid)

    if fmid < fa or fmid > fb:
        # Monotonicity violated, reached resolution of noise
        return (a, b)
    if fmid < 0:
        a, fa = mid, fmid
    else:
        b, fb = mid, fmid
    return noisy_bisect(f, a, b, fa, fb, tolerance)


# Improve approximation of root of f with Robbins Monro method
def Robbins_Monro(f, x0, n_steps, c_0, c_1=0.):

    x = x0
    for n in range(1, n_steps+1):
        a = c_0 / (c_1 + n)
        fx = f(x)
        x = x - a * fx

    return x


fa = target_function(args.fbound[0])
fb = target_function(args.fbound[1])
bracket = noisy_bisect(target_function, args.fbound[0], args.fbound[1],
                       fa, fb, 1e-3)
factor_estimate = 0.5 * sum(bracket)

print(f'field_scale_factor estimate: {factor_estimate:.4e}')

if args.refine_steps > 0:
    inv_deriv_estimate = (bracket[1] - bracket[0])/args.p
    factor_refined = Robbins_Monro(target_function, factor_estimate,
                                   args.refine_steps, inv_deriv_estimate)
    print(f'field_scale_factor refined: {factor_refined:.4e}')
