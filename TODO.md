TODO items
=======

Variable secondary emission coefficient for positive ions
==

Idea: Store the coefficient as a variable on the mesh, in front of the surfaces

Database with transport data
==

Include data for different gases, taking care to include 3-body attachment in air. For mixtures, we would need to couple with a Boltzmann solver. For some mixtures, the inclusion of additional reactions might be important (e.g. 3-body processes in humid air).

Database with secondary emission coefficients
==

Include estimates for different materials and ions.

Detachment from negative ions
==

Electrons lost to attachment (e.g., those produced by secondary emission near a surface) create a negative ion. In some gases, there can be efficient detachment from these negative ions. Maybe we can define a threshold field for such detachment and then include the resulting avalanches.

Photoionization generalization
==

Allow to set the follow parameters that control photoionization in any gas mixture:

* proportionality factor to number of ionizations
* quenching pressure
* absorption lengths (two assuming uniform distribution between these?)

Photoemission
==

* Detect photoionization photons hitting surfaces (but then what to do?)
* Can also include a separate photoemission process
