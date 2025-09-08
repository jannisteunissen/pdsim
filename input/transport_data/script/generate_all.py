#!/usr/bin/env python3

import numpy as np
import json

# Compontens, fractions, authors, three-body attachment factor
all_gases = [
    (['N2'], [1.0], ['Phelps'], None),
    (['N2'], [1.0], ['Trinity'], None),
    (['N2'], [1.0], ['IST_Lisbon'], None),
    (['N2'], [1.0], ['Biagi'], None),
    (['Ar'], [1.0], ['Biagi'], None),
    (['N2', 'O2'], [0.8, 0.2], ['Phelps', 'Phelps'], 1e-6),
    (['N2', 'O2'], [0.8, 0.2], ['Phelps', 'MuroranIT'], 1.0),
    (['N2', 'O2'], [0.8, 0.2], ['Trinity', 'MuroranIT'], 1.0),
    (['N2', 'O2', 'Ar'], [0.78, 0.21, 0.01],
     ['Phelps', 'Phelps', 'Phelps'], 1e-6),
    (['N2', 'O2', 'H2O'], [0.792, 0.198, 0.01],
     ['Phelps', 'Phelps', 'MuroranIT'], 1e-6),
    (['N2', 'O2', 'H2O'], [0.784, 0.196, 0.02],
     ['Phelps', 'Phelps', 'MuroranIT'], 1e-6),
    (['N2', 'O2', 'H2O'], [0.776, 0.194, 0.03],
     ['Phelps', 'Phelps', 'MuroranIT'], 1e-6),
    (['N2', 'O2', 'Ar', 'H2O'], [0.7722, 0.2079, 0.0099, 0.01],
     ['Phelps', 'Phelps', 'Phelps', 'MuroranIT'], 1e-6),
    (['N2', 'O2', 'Ar', 'H2O'], [0.7644, 0.2058, 0.0098, 0.02],
     ['Phelps', 'Phelps', 'Phelps', 'MuroranIT'], 1e-6),
    (['N2', 'O2', 'Ar', 'H2O'], [0.7566, 0.2037, 0.0097, 0.03],
     ['Phelps', 'Phelps', 'Phelps', 'MuroranIT'], 1e-6),
]


def generate_script(files, species, fractions, script_name, bolsig_output,
                    min_Td, max_Td, num_Td, scale_Td, extrapolate):

    if len(files) != len(species) or \
       len(files) != len(fractions):
        raise ValueError('Unequal number of species, files and fractions')

    if not np.isclose(sum(fractions), 1.0):
        raise ValueError('Fractions are not normalized')

    template = r'''
/NOSCREEN
/NOLOGFILE

_readcollisions_

CONDITIONS
10.       / Electric field / N (Td)
0.        / Angular field frequency / N (m3/s)
0.        / Cosine of E-B field angle
300.      / Gas temperature (K)
300.      / Excitation temperature (K)
0.        / Transition energy (eV)
0         / Ionization degree
0         / Gas particle density (1/m3)
1.        / Ion charge parameter
1.        / Ion/neutral mass ratio
0         / e-e momentum effects & modified Coulomb logarithm: 0=No&No; 1=Yes&No; 2=No&Yes; 3=Yes&Yes*
1         / Energy sharing: 1=Equal*; 2=One takes all
1         / Growth: 1=Temporal*; 2=Spatial; 3=Not included; 4=Grad-n expansion
0.        / Maxwellian mean energy (eV)
400       / # of grid points
0         / Manual grid: 0=No; 1=Linear; 2=Parabolic
200.      / Manual maximum energy (eV)
1e-10     / Precision
1e-4      / Convergence
1000      / Maximum # of iterations
_fractions_  / Gas composition fractions
1         / Normalize composition to unity: 0=No; 1=Yes

RUNSERIES
1          / Variable: 1=E/N; 2=Mean energy; 3=Maxwellian energy
_min_Td_ _max_Td_  / Min Max
_num_Td_         / Number
_scale_Td_          / Type: 1=Linear; 2=Quadratic; 3=Exponential

SAVERESULTS
_output_file_        / File
3        / Format: 1=Run by run; 2=Combined; 3=E/N; 4=Energy; 5=SIGLO; 6=PLASIMO
1        / Conditions: 0=No; 1=Yes
1        / Transport coefficients: 0=No; 1=Yes
1        / Rate coefficients: 0=No; 1=Yes
0        / Reverse rate coefficients: 0=No; 1=Yes
0        / Energy loss coefficients: 0=No; 1=Yes
0        / Distribution function: 0=No; 1=Yes
0        / Skip failed runs: 0=No; 1=Yes
0        / Include cross sections: 0=No; 1=Yes

END
'''

    template_readcollision = r'''
    READCOLLISIONS
    _cross_section_file_
    _species_                  / Species
    _extrapolate_              / Extrapolate: 0= No 1= Yes

    '''

    readcollisions = ''
    for s, f in zip(species, files):
        tmp = template_readcollision
        tmp = tmp.replace('_cross_section_file_', "'" + f + "'")
        tmp = tmp.replace('_species_', s)
        tmp = tmp.replace('_extrapolate_', str(extrapolate))
        readcollisions += tmp

    script = template
    script = script.replace('_readcollisions_', readcollisions)
    fractions = ' '.join([str(f) for f in fractions])
    script = script.replace('_fractions_', fractions)
    script = script.replace('_output_file_', bolsig_output)
    script = script.replace('_min_Td_', str(min_Td))
    script = script.replace('_max_Td_', str(max_Td))
    script = script.replace('_num_Td_', str(num_Td))
    script = script.replace('_scale_Td_', str(scale_Td))

    with open(script_name, 'w') as f:
        f.write(script)

if __name__ == '__main__':

    for species, fractions, authors, f_3body in all_gases:

        basename = '_'.join([f'{s}_{f}_{a}' for s, f, a in
                            zip(species, fractions, authors)])

        files = [f'../../cross_sections/{s}_{a}.txt' for s, a in
                 zip(species, authors)]

        bolsig_output = f'bolsig_result_{basename}.txt'
        script_name = f'bolsig_script_{basename}.txt'

        min_Td = 1.0
        max_Td = 1.5e3
        num_Td = 100
        scale_Td = 2
        extrapolate = 1

        generate_script(files, species, fractions, script_name, bolsig_output,
                        min_Td, max_Td, num_Td, scale_Td, extrapolate)

        metadata = {}

        metadata['species'] = species
        metadata['fractions'] = fractions
        metadata['authors'] = authors
        metadata['f_3body'] = f_3body

        metadata_file = bolsig_output + '.json'
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, ensure_ascii=False, indent=4)
