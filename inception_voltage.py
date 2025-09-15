#!/usr/bin/env python3

import numpy as np
import subprocess
import argparse
import pathlib
import os

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Utility to find the approximate inception voltage')
p.add_argument('cfg', type=str, help='Config file')
p.add_argument('-p', type=float, default=1e-3,
               help='Threshold inception probability')
p.add_argument('-fbound', type=float, nargs=2, required=True,
               help='Bounds for the field scale factor')
p.add_argument('-initial_steps', type=int, default=10,
               help='Number of steps for initial search')
p.add_argument('-refine_steps', type=int, default=0,
               help='Refine for this number of steps')
p.add_argument('-verbosity', type=int, default=0,
               help='How much information to print')
args = p.parse_args()


# Function to find a root of
def target_function(factor):
    this_folder = pathlib.Path(__file__).parent.resolve()
    pdsim_executable = os.path.join(this_folder, 'pdsim')
    proc = subprocess.run([pdsim_executable, args.cfg,
                           r'-output%verbosity=0',
                           r'-output%level=0',
                           r'-input%field_scale_factor=' + f'{factor}'],
                          capture_output=True, check=True)
    p_inception, _, _ = map(float, proc.stdout.split())
    return p_inception - args.p


# Assume a < b, fa = f(a) < 0, fb = f(b) > 0.
# Taken from: https://www.johndcook.com/blog/2012/06/14/root-finding-with-noisy-functions/
def noisy_bisect(f, a, b, fa, fb, tolerance, verbosity=0):
    if b-a < tolerance:
        return (a, b), (fa, fb)

    mid = 0.5 * (a+b)
    fmid = f(mid)

    if verbosity > 0:
        print(f'factor: {mid:.5e}, p_inception - threshold: {fmid:.5e}')

    if fmid < fa or fmid > fb:
        # Monotonicity violated, reached resolution of noise
        return (a, b), (fa, fb)
    if fmid < 0:
        a, fa = mid, fmid
    else:
        b, fb = mid, fmid
    return noisy_bisect(f, a, b, fa, fb, tolerance, verbosity)


# Improve approximation of root of f with Robbins Monro method
def Robbins_Monro(f, x0, n_steps, c_0, c_1=0., verbosity=0):
    x = np.zeros(n_steps+1)
    x[0] = x0

    for n in range(1, n_steps+1):
        a = c_0 / (c_1 + n)
        fx = f(x[n-1])

        if verbosity > 0:
            print(f'x: {x:.5e}, f(x): {fx:.5e}')

        x[n] = x[n-1] - a * fx

    return x[-1], x[-10:].std()


# Scan for lowest field_scale_factor for which the inception probability
# exceeds the threshold
factors = np.linspace(args.fbound[0], args.fbound[1], args.initial_steps)
val_lo, val_hi = None, None

for factor in factors:
    val = target_function(factor)

    if val < 0:
        factor_lo = factor
        val_lo = val
    else:
        factor_hi = factor
        val_hi = val
        break

if val_lo is None or val_hi is None:
    raise ValueError('The range fbound does not include a sign change')

bracket, values = noisy_bisect(target_function, factor_lo, factor_hi,
                               val_lo, val_hi, 1e-3, args.verbosity)

factor_estimate = 0.5 * sum(bracket)
err = 0.5 * (bracket[1] - bracket[0])

print(f'field_scale_factor estimate: {factor_estimate:.4e} +- {err:.1e}')

if args.refine_steps > 0:
    denom = values[1] - values[0]

    # Avoid division by a very small value
    if abs(denom) < args.p:
        denom = args.p * np.sign(denom)

    inv_deriv_estimate = (bracket[1] - bracket[0])/denom

    factor_refined, std = Robbins_Monro(target_function, factor_estimate,
                                        args.refine_steps, inv_deriv_estimate,
                                        verbosity=args.verbosity)
    print(f'field_scale_factor refined:  {factor_refined:.4e} +- {std:.1e}')
