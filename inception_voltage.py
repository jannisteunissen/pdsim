#!/usr/bin/env python3

import numpy as np
import subprocess
import argparse
import pathlib
import os
import sys

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Utility to find the approximate inception voltage')
p.add_argument('cfgs', type=str, nargs='+', help='Config file(s)')

g = p.add_mutually_exclusive_group(required=True)
g.add_argument('-p_avg', type=float,
               help='Threshold for average inception probability')
g.add_argument('-p_vol', type=float,
               help='Threshold for volume integral of inception probability')
g.add_argument('-p_nonzero', action='store_true',
               help='Search lowest voltage with non-zero inception probabily')

p.add_argument('-fbound', type=float, nargs=2, required=True,
               help='Bounds for the field scale factor')
p.add_argument('-initial_steps', type=int, default=10,
               help='Number of steps for initial search')
p.add_argument('-bisect_tolerance', type=float, default=1e-2,
               help='Stop bisection after interval is reduced by this factor')
p.add_argument('-n_runs', type=int, default=10,
               help='Use this number of runs per start location')
p.add_argument('-show_pmax', action='store_true',
               help='Show p_max for final estimated voltage')
p.add_argument('-settings', type=str, nargs='+',
               help=r'Other settings, e.g. "gas%%pressure=1.0" ' +
               r'"photoi%%enabled=T"')
p.add_argument('-verbosity', type=int, default=0,
               help='How much information to print')
args = p.parse_args()


# Function to find a root of
def target_function(factor):
    this_folder = pathlib.Path(__file__).parent.resolve()
    pdsim_executable = os.path.join(this_folder, 'pdsim')

    env = {}

    other_args = [
        r'-output%verbosity=0',
        r'-output%level=0',
        r'-avalanche%n_runs=' + f'{args.n_runs}',
        r'-input%field_scale_factor=' + f'{factor}'
    ]

    if args.settings is not None:
        other_args += ['-' + arg for arg in args.settings]

    if args.p_nonzero:
        other_args += ['-avalanche%use_early_exit=t']
        env['OMP_CANCELLATION'] = 'true'

    try:
        proc = subprocess.run(
            [pdsim_executable] + args.cfgs + other_args,
            capture_output=True,
            check=True,
            text=True,
            env=env
        )
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"\nError running pdsim:\n"
                         f"  Command: {' '.join(e.cmd)}\n"
                         f"  Return code: {e.returncode}\n"
                         f"  Stdout:\n{e.stdout}\n"
                         f"  Stderr:\n{e.stderr}\n")
        raise

    p_avg, stddev, volume, p_max = map(float, proc.stdout.split())

    if args.p_nonzero:
        p_max = 1 if p_avg > 0 else 0
        return (p_max, p_max)
    elif args.p_vol:
        return p_avg * volume - args.p_vol, p_max
    else:
        return p_avg - args.p_avg, p_max


# Assume a < b, fa = f(a) < 0, fb = f(b) > 0.
# Taken from: https://www.johndcook.com/blog/2012/06/14/root-finding-with-noisy-functions/
def noisy_bisect(f, a, b, fa, fb, tolerance, verbosity=0):
    if abs(b-a) < tolerance:
        return (a, b), (fa, fb)

    mid = 0.5 * (a+b)
    fmid, _ = f(mid)

    if verbosity > 0:
        print(f'Bisect x: {mid:12.5e}, f(x): {fmid:.5e}')

    if fmid < fa or fmid > fb:
        # Monotonicity violated, reached resolution of noise
        return (a, b), (fa, fb)

    if fmid <= 0:
        a, fa = mid, fmid
    else:
        b, fb = mid, fmid

    return noisy_bisect(f, a, b, fa, fb, tolerance, verbosity)


if args.fbound[0] * args.fbound[1] < 0:
    raise ValueError('fbound[0] and fbound[1] should have the same sign')

# Scan for lowest field_scale_factor for which the inception probability
# exceeds the threshold
sign = np.sign(args.fbound[0] + args.fbound[1])
a, b = sign * np.min(np.abs(args.fbound)), sign * np.max(np.abs(args.fbound))
factors = np.linspace(a, b, args.initial_steps)
val_lo, val_hi = None, None

for factor in factors:
    val, _ = target_function(factor)

    if args.verbosity > 0:
        print(f'Scanning x: {factor:.4e}, f(x): {val:.4e}')

    if val <= 0:
        factor_lo = factor
        val_lo = val
    else:
        factor_hi = factor
        val_hi = val
        break

if val_lo is None or val_hi is None:
    raise ValueError(f'The range fbound {args.fbound} does not lead to '
                     'a sign change')

bracket, _ = noisy_bisect(target_function, factor_lo, factor_hi,
                          val_lo, val_hi,
                          args.bisect_tolerance * (factor_hi - factor_lo),
                          args.verbosity)

factor_estimate = 0.5 * sum(bracket)
err = 0.5 * (bracket[1] - bracket[0])

if args.verbosity > 0:
    print('factor     error')

print(f'{factor_estimate:.4e} {err:.1e}')

if args.show_pmax:
    p, p_max = target_function(factor_estimatein)
    print(f'p_max: {p_max:.4f}')
