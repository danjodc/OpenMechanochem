import numpy as np
from ase.atoms import Atoms


class LinearPull(Atoms):

	""""Mechanochemical Simulation via Addition of Explicit Force
	    Parameters:
	    atoms: atoms object
	    method: string of method
	            * FMPES (Martinez) : J. Am. Chem. Soc. 2009, 131, 18, 6377–6379
		                   https://doi.org/10.1021/ja8095834
	            * EFEI (Marx)  : Angew. Chem. Int. Ed., 48, 4190 (2009)
		                   J. Am. Chem. Soc.132, 10609-10614 (2010)
		        * EGO (Baker) :   Molecular Physics,Vol 107, 22, (2009)
		                    Molecular Physicls, 1098, 14 (2010)  
	    pp: [[x1,y1,z1],[x2,y2,z2]]
	        list of xyz coordinates
	    ap: [int1,int2]
	        indices of atoms to be pulled
	        note: atom(int1) will be pulled to [x1,y1,z1]
	              and atom(int2) to [x2,y2,z2]
	    ego:[(int1, int2),(int1, int2),(int1, int2),...]
	         list of tuples
	         indices of atoms to be pulling
	         note: rotation and translation resulting to vector sum
	            should be considered.
	    force: float
	           applied force in nanoNewtons
	"""
	
	def set_params(self,method,pp,ap,ego = None, pullforce = 0.0):
		""" Initializes Linear Pulling Scheme """
		if method in ['FMPES', 'EFEI', 'EGO']:
			self.method = method
		else:
			raise NotImplementedError(method)

		if method =='FMPES' and len(pp) != 2:
			raise ValueError('Check Pulling Point Input')
		else:
			self.pp = pp

		if len(ap) != 2 and method != 'EGO':
			raise ValueError('Check Applied Point Input')
		else:
			self.ap = ap

		if method =='EGO':
			for elem in ego:
				if len(ego[elem]) != 2:
					raise ValueError('Check Pulling coordinates')
		else:
			self.ego = ego

		self.force = pullforce

	def get_forces(self, apply_constraint=True, md=False):
		""" Modified from ase.atoms.get_forces() """
		if self._calc is None:
			raise RuntimeError('Atoms object has no calculator.')
		forces = self._calc.get_forces(self) 

		if self.method == 'FMPES':
			forces = self.get_FMPESforces(forces)
		if self.method == 'EFEI':
			forces = self.get_EFEIforces(forces)
		else:
			forces = self.get_EGOforces(forces)

		if apply_constraint:
			for constraint in self.constraints:
				if md and hasattr(constraint,'redistribute_forces_md'):
					constraint.redistribute_forces_md(self, forces)
				if not md or hasattr(constraint, 'adjust_potential_energy'):
					constraint.adjust_forces(self, forces)
		return forces

	def get_FMPESforces(self,forces):
		""" Calculates Modified Forces according to FMPES Formalism """
		directionA = get_unitvector(self.get_positions()[self.ap[0]],
			                                    np.array(self.pp[0]))
		directionB = get_unitvector(self.get_positions()[self.ap[1]],
			                                    np.array(self.pp[1]))
		force_mag  = self.force / 2
		forceA     = force_mag * directionA
		forceB     = force_mag * directionB
		forces[self.ap[0]] = forces[self.ap[0]] + forceA
		forces[self.ap[1]] = forces[self.ap[1]] - forceB
		return forces

	def get_EFEIforces(self,forces):
		""" Calculates Modified Forces according to EFEI Formalism """
		direction = get_unitvector(self.get_positions()[self.ap[0]],
			                       self.get_positions()[self.ap[1]])
		force_mag = self.force / 2
		eff_force = force_mag * direction
		forces[self.ap[0]] = forces[self.ap[0]] - eff_force
		forces[self.ap[1]] = forces[self.ap[1]] + eff_force
		return forces

	def get_EGOforces(self,forces):
		""" Calculates Modified Forces according to EGO Formalism """
		force_mag = self.force / len(ego)
		for i in range(len(ego)):
			direction = get_unitvector(self.get_positions()[self.ego[i][0]],
				                       self.get_positions()[self.ego[i][1]])
			eff_force = force_mag * direction
			forces[self.ego[i][0]] = forces[self.ap[0]] - eff_force
			forces[self.ego[i][1]] = forces[self.ap[1]] + eff_force
		return forces


class WallPotential(Atoms):
	"""" Application of Invisible Wall Potential for Sandwich Pushing Simulations
	    Parameters:
	    atoms: atoms object
	    method: string of method
	                * linear      : r Potential     
		        * inverse     : 1 / r Potential
		        * lj          : 12-6 potential
       plane: [int1, int2, int3]
              Use to calculate a normal vector with respect to the plane
       height: float
               Plane height with respect to topmost atom
       force : float 
               linear  : equal to force
               inverse : multiplier to the -1/r^2
       params: list of len(2)
               epsilon
               sigma
	"""
	def set_params(self,method, plane, height, wallforce = None, params = None):
		""" Initializes Wall Potential """
		if method in ['linear', 'inverse', 'lj', 'repulsivelj']:
			self.method = method
		else:
			raise NotImplementedError(method)

		if len(plane) != 3:
			raise ValueError('Check Plane Definition')
		else:
			self.plane = plane

		self.height = height

		if (method =='linear' or  
			method == 'inverse' or 
			method == 'repulsivelj') and wallforce == None:
			raise ValueError('Check Force or Force Constant')
		else:
			self.wallforce = wallforce

		if method =='lj':
			raise NotImplementedError(method)


	def get_forces(self, apply_constraint=True, md=False):
		""" Modified from ase.atoms.get_forces() """
		if self._calc is None:
			raise RuntimeError('Atoms object has no calculator.')
		forces = self._calc.get_forces(self) 

		if self.method == 'linear':
			forces = self.get_linearforces(forces)
		if self.method == 'inverse':
			forces = self.get_inverseforces(forces)
		else:
			forces = self.get_ljforces(forces)

		if apply_constraint:
			for constraint in self.constraints:
				if md and hasattr(constraint,'redistribute_forces_md'):
					constraint.redistribute_forces_md(self, forces)
				if not md or hasattr(constraint, 'adjust_potential_energy'):
					constraint.adjust_forces(self, forces)
		return forces

	def get_linearforces(self,forces):
		normal = get_normalvector_plane(self.get_positions()[self.plane[0]],
			                            self.get_positions()[self.plane[1]],
			                            self.get_positions()[self.plane[2]])
		force_mag = self.wallforce * normal
		for item in range(len(forces)):
			forces[item] = forces[item] - force_mag
		return forces

	def get_inverseforces(self,forces):
		normal = get_normalvector_plane(self.get_positions()[self.plane[0]],
			                            self.get_positions()[self.plane[1]],
			                            self.get_positions()[self.plane[2]])
		for item in range(len(forces)):
			distance = get_planepoint_distance(normal,
				                               self.plane[0], 
				                               self.get_positions()[item])
			distance = distance + self.height
			force_mag = self.wallforce / distance**2 * normal
			forces[item] = forces[item] - force_mag
		return forces

	def get_ljforces(self,forces):
		normal = get_normalvector_plane(self.get_positions()[self.plane[0]],
			                            self.get_positions()[self.plane[1]],
			                            self.get_positions()[self.plane[2]])
		for item in range(len(forces)):
			distance = get_planepoint_distance(normal,self.plane[0], 
				                          self.get_positions()[item])
			distance = distance + self.height
			force_mag = 4 * ((-12/(distance**13))+(6/(distance**7)))
			forces[item] = forces[item] - force_mag
		return forces

def get_unitvector(atom1,atom2):
	""" Returns a unit vector of V(1->2) """
	dirvec = atom2-atom1
	length = np.linalg.norm(dirvec)
	dirvec = dirvec / length
	return dirvec	

def get_normalvector_plane(atom1, atom2, atom3):
	""" Returns a normal vector of plane containing atoms """
	dirvec1 = get_unitvector(atom1,atom2)
	dirvec2 = get_unitvector(atom1,atom3)
	normvec = np.cross(dirvec1, dirvec2)
	return normvec

def get_planepoint_distance(normal, planepoint, atompoint):
	""" Returns shortest distance between plane and point """
	val = atompoint - planepoint
	distance = np.dot(normal, val)
	return distance

