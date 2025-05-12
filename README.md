# pdsim: Partial Discharge Simulation

The goal of this tool is to simulate the inception of partial discharges in a given electric field configuration. Most of the work still needs to be done, but these are the intended features:

1. Support loading of different types of unstructured grids that contain the electric field distribution
2. Support for classic particle-in-cell type model for electrons as well as a probabilistic avalanche model
3. Support for ions and photoionization, which both can cause secondary avalanches
4. Output: the probability of discharge inception per initial electron (assuming a uniform distribution), where discharge inception is defined as the number of electrons exceeding a threshold

## Obtaining an electric field distribution

1. The user computes an electric field profile using finite element software
2. Using the `interpolate_unstructured/convert_to_binary` script, the unstructured grid is converted into a format the code can parse
