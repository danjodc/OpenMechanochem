# Mechanochemistry

## Introduction
Mechanochemistry module includes functions that can be used to simulate mechanochemical phenomena via the addition of external forces.Cancel Changes
This is designed to be used with the *Atomic Simulation Environment* (ASE). Together with ASE, this module can be combined with various quantum mechanical calculators for geometry optimizations and nudged elastic band calculations to sample the effect of mechanical forces to the potential energy hypersurface. Currently, this module can execute:

* **Force Modified Potential Energy Surface** (FMPES) :<br/>
 	J. Am. Chem. Soc., 131, 18, 6377–6379 (2009) 
* **External Force Explicitly Included** (EFEI) :<br/>
	Angew. Chem. Int. Ed., 48, 4190 (2009)<br/>
 	J. Am. Chem. Soc.132, 10609-10614 (2010)
* **Enforced Geometry Optimization** (EGO) :<br/>
 	Molecular Physics,Vol 107, 22, (2009)<br/>
 	Molecular Physics, 1098, 14 (2010)  
<sub><sup>*Disclaimer: The author of this repository is not affiliated to the proponents of formalism above. Code development and testing was done independently*<sub><sup>	

Additionally, this module can also be used to add a **finite wall potential in single direction**. Current implementations includes four different types of potential, (1) linear, (2) inverse, (3) lennard-jones, and (4) $1/r^6$ potentials.

## Usage and Citation

Citation can be done as
```
De Chavez, D., Mechanochem, (2021), GitHub repository, 
	https://github.com/danjodc/Mechnochem
```

or similary depending on your citation style. Alternatively, the BibTex library can be appended with

```
@misc{DeChavez2021,
  author = {De Chavez, Danjo.},
  title = {Mechanochem},
  year = {2021},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/danjodc/mechanochem}},
}
```

## Tutorial

## Function Requests

For function request related to mechanochemistry simulations and force analysis tools please contact author.
