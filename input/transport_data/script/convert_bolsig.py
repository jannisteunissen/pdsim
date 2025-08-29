#!/usr/bin/env python3

import numpy as np
import argparse


def get_data(lines, key, n_skip=0):
    try:
        i0 = lines.index(key)
    except ValueError:
        return None

    # Skip current line and optionally additional lines
    i0 = i0 + 1 + n_skip
    i1 = i0 + lines[i0:].index('')
    text = ' '.join(lines[i0:i1])
    return np.fromstring(text, sep=' ')


def get_attachment_rates(lines):
    keys = []
    data = {}

    for line in lines:
        if line.startswith('C') and 'Attachment' in line:
            keys.append(line)

    for key in keys:
        data[key] = get_data(lines, key, 1)

    return data


def write_entry(entry_name, x_data, y_data, f):
    if y_data is None:
        y_data = np.zeros_like(x_data)

    out_data = np.column_stack([x_data, y_data])

    hdr = entry_name
    hdr += '\n-----------------------\n'
    ftr = '-----------------------\n\n'
    f.write(hdr.encode('ascii'))
    np.savetxt(f, out_data)
    f.write(ftr.encode('ascii'))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description="Converter for Bolsig+ data")
    p.add_argument("infile", type=str, help="input file")
    p.add_argument("outfile", type=str, help="output file")
    args = p.parse_args()

    keys = [
        'Electric field / N (Td)',
        'Mobility *N (1/m/V/s)',
        'Townsend ioniz. coef. alpha/N (m2)',
        'Townsend attach. coef. eta/N (m2)'
    ]

    with open(args.infile, 'r') as f:
        lines = f.read().splitlines()
        lines = [line.strip() for line in lines]

    data = {}

    for key in keys:
        data[key] = get_data(lines, key)

    # Check for three-body attachment in filename
    tmp = args.infile.split('_3ba_')

    if len(tmp) > 1:
        three_body_factor = float(tmp[1].split('.')[0])
        attachment_rates = get_attachment_rates(lines)

        # Find attachment process with lowest rate
        max_rates = {key: rate.max() for key, rate in
                     attachment_rates.items()}
        lowest = min(max_rates, key=max_rates.get)

        data['Three-body attachment rate (m6/s)'] = \
            three_body_factor * attachment_rates[lowest]
    else:
        data['Three-body attachment rate (m6/s)'] = None

    with open(args.outfile, 'wb') as f:
        x_key = 'Electric field / N (Td)'

        for key in data:
            if key != x_key:
                write_entry(key, data[x_key], data[key], f)

    print(f'Wrote {args.outfile}')
