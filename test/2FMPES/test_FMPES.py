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
 H                 0.000   0.000   0.000    |  0  |   
 H                 0.000   0.000   1.000    |  1  |

Pulling Points 
                   0.000   0.000 -20.000    |  A  |
                   0.000   0.000  30.000    |  B  |                                                      

"""

PullPoints = [[0.000,  0.000, -1.000],
              [0.000,  0.000,  2.000]]

AppPoints  = [0,1]

pull = mc.LinearPull(mol)

"""   
Physical Model             

   A <---  0 ------- 1  --->   B
         
"""

for i in np.linspace(0,10,21):
  force = i 
  pull.set_params(method='FMPES', pp=PullPoints, ap=AppPoints, pullforce = force)
  pull.set_calculator(EMT())
  dyn = BFGS(pull, trajectory='2FMPES_Hydrogen_{}nN.traj'.format(i))
  dyn.run(fmax=0.05)


