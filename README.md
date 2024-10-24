# geometry-driven-hyperuniformity

## Background

Geometry-driven hyperuniformity refers to the state of the resulting configuration obtained through an ordering transition along the Lloyd's algorithm.
Given a point pattern as input, the algorithm computes the centroid of the Voronoi cell of each point to replace every point by the associated centroid at each iteration step. 
As it repeates the process, the point pattern exhibits effective hyperuniformity (Klatt *et al*. 2019) that is characterised by an anomalous suppression of density fluctuations at large length scales.

The animation below illustrates the ordering transition of a cellular pattern up to Lloyd iteration step 300, where vectors represent cell orientations $\Theta$ and (Voronoi) cells at each iteration step are coloured according to $\Theta$. 

![](https://github.com/YeoniH/geometry-driven-hyperuniformity/blob/main/N100_Poi-1_t0-300_loop.gif)

Upon convergence towards hyperuniformity, the Voronoi landscape consists of almost regular hexagonal cells with a small fraction of pentagonal and heptagonal defects.

## About this repository

This repository contains Python codes for computational analysis on evolving point configurations, where a randomly disordered point pattern is taken as input to the Lloyd's algorithm to transition into an effectively hyperuniform configuration.
The provided codes form the basis of the investiagation on "geometry-driven hyperuniformity" in the [doctoral thesis](https://openresearch-repository.anu.edu.au/items/71d9e451-34a6-42bf-a389-cb56188bbc18) of [Sungyeon Hong](https://cybernetics.anu.edu.au/people/sungyeon-hong/) along with her publications below.

The Jupyter Notebooks are provided for users to run the Lloyd's iterations and conduct computational analyses using functions defined in Analysis module.
For now, the simulation is only for 2D, but the 3D counterpart will be included in the near future — Stay tuned!

To smoothly run the codes, the following packages are required:
* [numpy](https://numpy.org/)
* [scipy](https://docs.scipy.org/doc/scipy/)
* [freud-anlysis](https://freud.readthedocs.io/en/stable/#)
* [collections](https://docs.python.org/3/library/collections.html)

Also, some external software will be in need for the computation of structural entropy in publication [3]:
* [nauty](https://pallini.di.uniroma1.it/)

## References
[1] S. Hong, M. A. Klatt, G. Schröder-Turk, N. François, M. Saadatfar. (2021). [Dynamical arrest of topological defects in 2D hyperuniform disk packings](https://www.epj-conferences.org/articles/epjconf/abs/2021/03/epjconf_pg2021_15002/epjconf_pg2021_15002.html). In Proceedings of the Powders and Grains 2021. EPJ Web of Conferences 249, 15002.

[2] M. A. Klatt, J. Lovrić, D. Chen et al. (2019) [Universal hidden order in amorphous cellular geometries](https://doi.org/10.1038/s41467-019-08360-5). Nat Commun 10, 811.

[3] S. Hong, C. Nerse, S. Oberst, M. Saadatfar. (in press). Topological mechanical states in geometry-driven hyperuniform materials. PNAS Nexus
