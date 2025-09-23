#!/usr/bin/env python3

import numpy as np
import subprocess
import argparse
import pathlib
import os

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Utility to find the approximate inception voltage')
p.add_argument('cfgs', type=str, nargs='+', help='Config file(s)')
p.add_argument('-p', type=float, default=1e-4,
               help='Threshold inception probability')
p.add_argument('-fbound', type=float, nargs=2, required=True,
               help='Bounds for the field scale factor')
p.add_argument('-initial_steps', type=int, default=10,
               help='Number of steps for initial search')
p.add_argument('-refine_steps', type=int, default=0,
               help='Refine for this number of steps')
p.add_argument('-n_runs', type=int, default=10,
               help='Use this number of runs per start location')
p.add_argument('-verbosity', type=int, default=0,
               help='How much information to print')
args = p.parse_args()


# Function to find a root of
def target_function(factor):
    this_folder = pathlib.Path(__file__).parent.resolve()
    pdsim_executable = os.path.join(this_folder, 'pdsim')

    other_args = [
        r'-output%verbosity=0',
        r'-output%level=0',
        r'-avalanche%n_runs=' + f'{args.n_runs}',
        r'-input%field_scale_factor=' + f'{factor}'
    ]

    proc = subprocess.run([pdsim_executable] + args.cfgs + other_args,
                          capture_output=True, check=True)
    p_inception, _, _ = map(float, proc.stdout.split())
    return p_inception - args.p


# Assume a < b, fa = f(a) < 0, fb = f(b) > 0.
# Taken from: https://www.johndcook.com/blog/2012/06/14/root-finding-with-noisy-functions/
def noisy_bisect(f, a, b, fa, fb, tolerance, verbosity=0):
    if abs(b-a) < tolerance:
        return (a, b), (fa, fb)

    mid = 0.5 * (a+b)
    fmid = f(mid)

    if verbosity > 0:
        print(f'Bisect x: {mid:12.5e}, f(x): {fmid:.5e}')

    if fmid < fa or fmid > fb:
        # Monotonicity violated, reached resolution of noise
        return (a, b), (fa, fb)
    if fmid < 0:
        a, fa = mid, fmid
    else:
        b, fb = mid, fmid
    return noisy_bisect(f, a, b, fa, fb, tolerance, verbosity)


# Improve approximation of root of f with Robbins Monro method
def Robbins_Monro(f, x0, n_steps, inv_deriv, c_1, verbosity=0,
                  fx_abs_limit=None):
    x = np.zeros(n_steps+1)
    fx = np.zeros(n_steps)
    x[0] = x0
    n_sign_changes = 1

    for n in range(1, n_steps+1):
        fx[n-1] = f(x[n-1])

        if n > 2 and fx[n-1] * fx[n-2] < 0:
            n_sign_changes += 1

        a = inv_deriv * (1 + c_1) / (n_sign_changes + c_1)

        if verbosity > 0:
            print(f'{n:4d} x: {x[n-1]:12.5e}, f(x): {fx[n-1]:12.5e}, '
                  f'a: {a:.3e}')

        if fx_abs_limit is not None:
            # Limit magnitude. Can be helpful for asymmetric functions
            fx[n-1] = np.clip(fx[n-1], -fx_abs_limit, fx_abs_limit)

        x[n] = x[n-1] - a * fx[n-1]

    return x[-1], x[-n_steps//2:].std()


if args.fbound[0] * args.fbound[1] < 0:
    raise ValueError('fbound[0] and fbound[1] should have the same sign')

# Scan for lowest field_scale_factor for which the inception probability
# exceeds the threshold
sign = np.sign(args.fbound[0] + args.fbound[1])
a, b = sign * np.min(np.abs(args.fbound)), sign * np.max(np.abs(args.fbound))
factors = np.linspace(a, b, args.initial_steps)
val_lo, val_hi = None, None

for factor in factors:
    val = target_function(factor)

    if args.verbosity > 0:
        print(f'Scanning x: {factor:.4e}, f(x): {val:.4e}')

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

if args.verbosity > 0:
    print(f'factor estimate: {factor_estimate:.4e} {err:.1e}')

if args.refine_steps > 0:
    inv_deriv_est = (factor_hi - factor_lo) / \
        max(abs(val_hi - val_lo), args.p)

    factor_estimate, std = Robbins_Monro(target_function, factor_estimate,
                                        args.refine_steps, inv_deriv_est,
                                        0, verbosity=args.verbosity,
                                        fx_abs_limit=2*args.p)

if args.verbosity > 0:
    print('p_threshold factor     error')

print(f'{args.p:.4e}  {factor_estimate:.4e} {err:.1e}')
