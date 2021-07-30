from ase import Atoms
from ase.io import read
from ase.optimize import BFGS
from ase.calculators.emt import EMT

import numpy as np
import OpenMechanochem as mc

mol = read('hydrogen.xyz')

"""
--> hydrogen.xyz
2
Hydrogen
 H                 0.000  0.000  0.000    |  0  |   
 H                 0.000  0.000  1.000    |  1  |

"""

PullPoints = [[0.000,  0.000, -2.000],
              [0.000,  0.000,  3.000]]

AppPoints  = [0,1]

pull = mc.LinearPull(mol)

"""   
Physical Model             

	     0 <-------> 1  
    	   
"""

for i in np.linspace(0,10,21):
  force = i #/ 1.60217733
  pull.set_params(method='EFEI', pp=None, ap=AppPoints, pullforce = force)
  pull.set_calculator(EMT())
  dyn = BFGS(pull, trajectory='1EFEI_Hydrogen_{}nN.traj'.format(force))
  dyn.run(fmax=0.05)


