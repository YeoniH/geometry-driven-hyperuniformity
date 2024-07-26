# geometry-driven-hyperuniformity

This repository contains Python codes for computational analysis on evolving point configurations, where a randomly disordered point pattern is taken as input to the Lloyd's algorithm to transition into an effectively hyperuniform configuration.
Geometry-driven hyperuniformity refers to the resulting configuration obtained through such an ordering transition along the Lloyd's algorithm. 

The Jupyter Notebooks are provided for users to run the Lloyd's iterations and conduct computational analyses using functions defined in Analysis module.
For now, the simulation is only for 2D, but the 3D counterpart will be included in the near future — Stay tuned!
The provided code form the basis of the investiagation on "geometry-driven hyperuniformity" in the [doctoral thesis](https://openresearch-repository.anu.edu.au/items/71d9e451-34a6-42bf-a389-cb56188bbc18) of [Sungyeon Hong](https://cybernetics.anu.edu.au/people/sungyeon-hong/) along with her publications below.

To smoothly run the codes, the following packages are required:
* [numpy](https://numpy.org/)
* [scipy](https://docs.scipy.org/doc/scipy/)
* [freud-anlysis](https://freud.readthedocs.io/en/stable/#)
* [collections](https://docs.python.org/3/library/collections.html)

References:
* S. Hong, M. A. Klatt, G. Schröder-Turk, N. François, M. Saadatfar. (2021). [Dynamical arrest of topological defects in 2D hyperuniform disk packings](https://www.epj-conferences.org/articles/epjconf/abs/2021/03/epjconf_pg2021_15002/epjconf_pg2021_15002.html). In Proceedings of the Powders and Grains 2021. EPJ Web of Conferences 249, 15002.
* M. A. Klatt, J. Lovrić, D. Chen et al. (2019) [Universal hidden order in amorphous cellular geometries](https://doi.org/10.1038/s41467-019-08360-5). Nat Commun 10, 811.
