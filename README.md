# pdsim: Monte Carlo simulation of partial discharge inception

This tool can be used to simulate the inception of partial discharges in a given electric field configuration. Given an electric field on an unstructured grid, the code can compute the probability of discharge inception for any initial electron position. The main features are:

* A simplified electron avalanche model, that allows for very fast sampling of avalanches (on the order of millions per second).
* Secondary electron emission due to ions, photoionization or photoemission can be included, with the option of having different emission coefficients for different types of materials.
* A conventional particle model for electrons can be used to test the assumptions in the simplified electron avalanche model.
* The probability of discharge inception and the average inception time are estimated for any initial electron position, and stored on the unstructured grid.

## Installation

Compilation requires a recent Fortran compiler, such as `gfortran`. If another compiler is used, the `Makefile` might have to be edited.

Step 1: Download the source code:

    git clone --recurse-submodules https://github.com/jannisteunissen/pdsim.git

Step 2: compile the two libraries that are required:

    make -C particle_core
    make -C interpolate_unstructured

Step 3: compile the main program `pdsim`:

    make

## General usage

### Computing an electric field

Compute an electric field profile using finite element software. The field should be defined at the vertices/points of the unstructured grid. Software such as COMSOL Multiphysics can for example be used.

### Converting the unstructured grid to binary data

Using the `interpolate_unstructured/convert_to_binary.py` script, convert the unstructured grid to a simple binary format that the code can parse:

    interpolate_unstructured/convert_to_binary.py <path_to_grid>

### Create a configuration file

To see which settings can be set in the configuration, run the following:

    mkdir output
    ./pdsim

The code should stop with an error message, but it produces a file `output/pdsim_config.cfg` that shows all the settings and their default values. This file can be used as a template for a new configuration file, and only the parameters that need to be changed have to be specified.
For a description of the syntax of the configuration files, see [config_fortran](https://github.com/jannisteunissen/config_fortran).

### Input data

As input, the model requires the following input data:

* Reduced ionization coefficient `alpha/N` (unit: `m^2`)
* Reduced attachment coefficient `eta/N` (unit: `m^2`)
* Reduced electron mobility `mu*N` (unit: `1/(m V s)`), only used for estimating the inception time.
* Three-body attachment rate (unit: `m^6/s`), for example relevant for air-like mixtures

These coefficients can be computed from input cross sections using a Boltzmann solver. Several examples can be found in the `input/transport_data` folder.

### Run a simulation

There are three main simulation modes that can be specified in the configuration file:

1. `simulate = integral`: computes ionization integrals and related quantities (this should always be fast)
2. `simulate = avalanches`: simulates avalanches for any initial electron position. The cost depends on `avalanche%n_runs`, the number of points in the mesh, and how likely inception is for these points.
3. `simulate = particles`: this simulates electron avalanches starting from a given initial position, and can be used to test the simplified avalanche model used in the code.

After the desired mode has been specified, the simulation can be started with:

    ./pdsim <config_file>

A `.vtu` file with the following output will then be saved:

* `inception_prob`: the estimated inception probability per initial location
* `inception_time`: the average inception time per initial location
* `alpha`, `eta`, `alpha_eff`: ionization, attachment and effective ionization coefficient
* `avalanche_p_m1`: probability that the initial electron produces no additional ionization
* `K_star`: logarithm of the expected avalanche size (total number of ionizations produced), usually very similar to `K_integral`
* `K_integral`: integral over `alpha_eff` along a field line
* `avalanche_time`: estimated travel time for an avalanche
* `avalanche_x1`, `avalanche_x2`, `avalanche_x3`: location where the avalanche reaches half its final size
* `ion_time`: estimated travel time for positive ions to hit a boundary
* `ion_x1`, `ion_x2`, `ion_x3`: location where positive ions reach a boundary
* `ion_gamma`: secondary electron emission coefficient for positive ions starting from this location, depends on the type of boundary that they will hit

This file can be visualized with tools such as Visit, pyvista or Paraview.

### Finding the partial discharge inception voltage (PDIV)

The `inception_voltage.py` script can be used to automatically find the partial discharge inception voltage (PDIV). As input, the script requires an upper and lower bound for the "field scale factor", with which the electric field on the input mesh is multiplied, and a configuration file.

# Examples

## Inception around a sphere (2D axisymmetric)

This example can be found under `examples/sphere_2d`. A `.vtu` file is provided as input, which first has to be converted:

    cd examples/sphere_2d

    ../../interpolate_unstructured/convert_to_binary.py sphere_axisymmetric_2d.vtu

Then run the example with:

    ../../pdsim sphere_2d_air.cfg

An example of the output is shown below

![inception probability](https://github.com/jannisteunissen/pdsim/blob/main/documentation/sphere_2d_example.png?raw=true)

## Other examples

Other examples can be run in the same way as the above one:

* `examples/plate_plate_2d`: Inception in plate-plate 2D geometry
