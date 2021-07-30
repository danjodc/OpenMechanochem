# OpenMechanochem

## Introduction
OpenMechanochemistry includes functions that can be used to simulate mechanochemical phenomena via the addition of external forces.
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

Additionally, this module can also be used to add a **finite wall potential in single direction**. Current implementations includes four different types of potential, (1) linear, (2) inverse, and (3) lennard-jones potential.

##  Citations

Citation can be done as
```
De Chavez, D., Mechanochem, (2021), GitHub repository, 
	https://github.com/danjodc/OpenMechanochem
```

or similarly depending on your citation style. For LaTex users, the BibTex library can be appended with

```
@misc{DeChavez2021,
  author = {De Chavez, Danjo.},
  title = {OpenMechanochem},
  year = {2021},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/danjodc/mechanochem}},
}
```

## Usage and Tutorial

### LinearPull Class
In this tutorial, we will use FMPES and EFEI formalism to pull hydrogen molecule along the bond coordinate.

The mechanochem classes LinearPull and WallPot inherit from the ASE atoms object and hence an Atoms instance is required.  
For an instance named mol, this can easily be done as

```
import mechanochem as mc

pull = mc.LinearPull(mol)
```

At this point, user should provide specific keyword arguments for the LinearPull object which would depend on the method key.
As stated earlier the Mechanochem classes inherits from Atoms object and accepts similar parameters such as calculator, pbc, etc.

#### FMPES

Using FMPES, the required parameters are pulling points, applied points and applied forces 
This can be done by using the keywords pp, ap, and pullforce respectively.

```
  pull.set_params(method='FMPES', pp=PullPoints, ap=AppPoints, pullforce = force)
```
In FMPES, the relative cartesian coordinates of pulling points and applied points are of utmost importance.
For example, a system described below where H0 and H1 are hydrogen atoms pulled towards points A and B respectively, 

```
      A <---  0 ------- 1  --->   B
```

with H2 xyz given as 
```
--> hydrogen.xyz
2
Hydrogen
 H                 0.000  0.000  0.000  
 H                 0.000  0.000  1.000

```
The pulling points and applied points could be given as
```
PullPoints = [[0.000,  0.000, -1.000],
              [0.000,  0.000,  2.000]]

AppPoints  = [0,1]
```
Care should be given that the position in pp and ap list corresponds with each other.  
That is, the first list in pp list is the direction where atom with index as the first element of ap list is pulled to.  

The magnitude of the applied force can be controlled using the pullforce key. 
Note that the force provided should be in atomic units and the pullforce is divided in the two force vectors equally.

#### EFEI

In comparison to FMPES, the EFEI pulls along the internal molecular coordinates.
Using this method, the pulling coordinate can be defined using only the applied points and pullforce.

```
  pull.set_params(method='EFEI', ap=AppPoints, pullforce = force)
```

```
	     H0 <-------> H1  
```

The ap list is same as the case above, which is
```
AppPoints  = [0,1]
```

Using EFEI formalism, the pull force is divided equally to the two atoms given by the ap list.

After parameterization of the LinearPull class, geometry optimization or molecular dynamics can be done.  
In the case of optimization, an example is given below.

```
pull.set_calculator(EMT())
dyn = BFGS(pull, trajectory='optimization.traj')
dyn.run(fmax=0.05)
```

Molecular dynamics can be done similarly. Users are suggested to visit ASE documentations for descriptions of parameters needed for MD and optimizations.

### WallPotential Class

Similar to LinearPull class, the wallpot inherits from the atoms class.
Hence, a prior instance of atoms should be provided.

```
mol = mc.WallPotential(slab)
```
The WallPotential class takes the parameters method, plane, height, and wallforce.

```
mol.set_params(method='linear', plane=atomplane , height=10, wallforce=force)
```
Method defines the type of interacting potential with the wall.
The plane is a list of len(3) which defines the equation of the plane.
The imaginary plane can then be displaced along the z direction with height key.
Similar to LinearPull class, the magnitude of force can be controlled by wallforce key.  

## Function Requests

For function request related to mechanochemistry simulations and force analysis tools please contact author.
