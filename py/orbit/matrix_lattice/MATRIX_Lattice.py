"""
The general Matrix lattice. The matrices track the bunch as 
the linear transport elements. It is a base class for TEAPOT_MATRIX_Lattice,
but it can be used by itself if user specifies the transport matrices by using 
addNode(BaseMATRIX). After adding nodes the user has to initialize() the lattice.
The initialization will create one turn matrix, but the user can ignore it if
the lattice is a linac lattice.
The  This class cannot calculate chromaticities.
"""
import os
import math

#import bunch
from bunch import Bunch

# import the function that creates multidimensional arrays
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import matrix class and generators
from orbit.teapot_base import MatrixGenerator
from orbit_utils import Matrix
from orbit_utils import PhaseVector

# import the AccNode implementation for a transport matrix
from BaseMATRIX import BaseMATRIX

class MATRIX_Lattice(AccLattice):
	"""
	The subclass of the AccLattice class. Shell class for the BaseMATRIX nodes collection.
	In the beginning the lattcie is empty.
	"""
	def __init__(self, name = None):
		AccLattice.__init__(self,name)
		self.oneTurmMatrix = Matrix(7,7)
		self.oneTurmMatrix.unit()		
		
	def initialize(self):
		"""
		Method. Initializes the matrix lattice, child node structures, and calculates 
		the one turn matrix.
		"""
		AccLattice.initialize(self)	
		self.makeOneTurnMatrix()
		
	def makeOneTurnMatrix(self):
		"""
		Calculates the one turn matrix.
		"""
		self.oneTurmMatrix.unit()
		for matrixNode in self.getNodes():
			if(isinstance(matrixNode,BaseMATRIX) == True):
				self.oneTurmMatrix = matrixNode.getMatrix().mult(self.oneTurmMatrix)
		return self.oneTurmMatrix

	def getOneTurnMatrix(self):
		"""
		Returns the one turn matrix.
		"""
		return self.oneTurmMatrix		
		
	def getRingParametersDict(self,momentum,mass):
		"""
		Returns the dictionary with different ring parametrs
		calculated from the one turn transport matrix.
		"""
		res_dict = {}
		Etotal = math.sqrt(momentum**2 + mass**2)
		beta = momentum/Etotal
		gamma = Etotal/mass
		Ekin = Etotal - mass
		res_dict["momentum [GeV/c]"] = momentum
		res_dict["mass [GeV]"] = mass
		res_dict["Ekin [GeV]"] = Ekin
		#longitudinal params
		c = 2.99792458e+8
		ring_length = self.getLength()
		T = ring_length/(beta*c)
		res_dict["period [sec]"] = T
		res_dict["frequency [Hz]"] = 1./T
		# the minus sign needed to take into account dz = - c*beta*dt
		eta_ring = - self.oneTurmMatrix.get(4,5)*beta*beta*(Ekin+mass)/ring_length
		res_dict["eta"] = eta_ring
		res_tmp = 1.0/(eta_ring + 1./(gamma*gamma))
		gamma_trans = math.sqrt(math.fabs(res_tmp))
		if(res_tmp < 0.): gamma_trans = complex(0.,gamma_trans)
		res_dict["gamma transition"] = gamma_trans
		res_dict["transition energy [GeV]"] = (gamma_trans - 1.0)*mass
		alpha_p = eta_ring + 1./(gamma*gamma)
		res_dict["momentum compaction"] = alpha_p
		#transverse twiss parameters
		mt = self.oneTurmMatrix
		res_dict["fractional tune x"] = None
		res_dict["fractional tune y"] = None
		res_dict["alpha x"] = None
		res_dict["alpha y"] = None		
		res_dict["beta x [m]"] = None
		res_dict["beta y [m]"] = None
		res_dict["gamma x [m^-1]"] = None
		res_dict["gamma y [m^-1]"] = None
		res_dict["dispersion x [m]"] = None
		res_dict["dispersion y [m]"] = None
		res_dict["dispersion prime x"] = None
		res_dict["dispersion prime y"] = None
		cos_phi_x = (mt.get(0,0)+mt.get(1,1))/2.0
		cos_phi_y = (mt.get(2,2)+mt.get(3,3))/2.0		
		if(abs(cos_phi_x) >= 1.0 or abs(cos_phi_x) >= 1.0):
			txt = "Class orbit.matrix_lattice.MATRIX_Lattice."+os.linesep
			txt = txt + "Method getRingParametersDict(momentum,mass):"+os.linesep
			txt = txt + "The one turn matrix is unstable: Sp(Mx) or Sp(My) > 2."
			orbitFinalize(txt)
			return res_dict
		sign_x = +1.0
		if(abs(mt.get(0,1)) != 0.): sign_x = mt.get(0,1)/abs(mt.get(0,1))
		sign_y = +1.0
		if(abs(mt.get(2,3)) != 0.): sign_y = mt.get(2,3)/abs(mt.get(2,3))
		sin_phi_x = math.sqrt(1. - cos_phi_x*cos_phi_x)*sign_x
		sin_phi_y = math.sqrt(1. - cos_phi_y*cos_phi_y)*sign_y
		nux = math.acos(cos_phi_x)/(2*math.pi)*sign_x
		nuy = math.acos(cos_phi_y)/(2*math.pi)	*sign_y
		res_dict["fractional tune x"] = nux
		res_dict["fractional tune y"] = nuy
		# alpha, beta, gamma
		beta_x = mt.get(0,1)/sin_phi_x 
		beta_y = mt.get(2,3)/sin_phi_y 
		alpha_x = (mt.get(0,0) - mt.get(1,1))/(2*sin_phi_x)
		alpha_y = (mt.get(2,2) - mt.get(3,3))/(2*sin_phi_y)
		gamma_x = -mt.get(1,0)/sin_phi_x
		gamma_y = -mt.get(3,2)/sin_phi_y
		# dispersion and dispersion prime
		m_coeff =  momentum*momentum/(mass + Ekin)
		disp_x = m_coeff*(mt.get(0,5)*(1-mt.get(1,1))+mt.get(0,1)*mt.get(1,5))/(2-mt.get(0,0)-mt.get(1,1)) 
		disp_y = m_coeff*(mt.get(2,5)*(1-mt.get(3,3))+mt.get(2,3)*mt.get(3,5))/(2-mt.get(2,2)-mt.get(3,3)) 
		disp_pr_x = m_coeff*(mt.get(1,0)*mt.get(0,5)+mt.get(1,5)*(1-mt.get(0,0)))/(2-mt.get(0,0)-mt.get(1,1))
		disp_pr_y = m_coeff*(mt.get(3,2)*mt.get(2,5)+mt.get(3,5)*(1-mt.get(2,2)))/(2-mt.get(2,2)-mt.get(3,3))
		res_dict["alpha x"] = alpha_x
		res_dict["alpha y"] = alpha_y
		res_dict["beta x [m]"] = beta_x
		res_dict["beta y [m]"] = beta_y
		res_dict["gamma x [m^-1]"] = gamma_x
		res_dict["gamma y [m^-1]"] = gamma_y
		res_dict["dispersion x [m]"] = disp_x
		res_dict["dispersion y [m]"] = disp_y 
		res_dict["dispersion prime x"] = disp_pr_x
		res_dict["dispersion prime y"] = disp_pr_y
		return res_dict

	def getRingTwissDataX(self,momentum,mass):
		"""
		Returns the tuple (tuneX, [(position, alphaX),...],[(position,betaX),...] ). 
		"""
		res_dict = self.getRingParametersDict(momentum,mass)
		alpha_x = res_dict["alpha x"]
		beta_x = res_dict["beta x [m]"]
		return self.trackTwissData(alpha_x,beta_x,"x")
	
	def getRingTwissDataY(self,momentum,mass):
		"""
		Returns the tuple (tuneY, [(position, alphaY),...],[(position,betaY),...] ).
		"""
		res_dict = self.getRingParametersDict(momentum,mass)
		alpha_y = res_dict["alpha y"]
		beta_y = res_dict["beta y [m]"]
		return self.trackTwissData(alpha_y,beta_y,"y")

	def trackTwissData(self, alpha, beta, direction = "x"):
		"""
		Returns the tuple (tune, [(position, alpha),...],[(position,beta),...] ). 
		The tracking starts from the values specified as the initial parameters. 
		The possible values for direction parameter "x" or "y".
		"""
		if(direction.lower() != "x" and direction.lower() != "y"):
			orbitFinalize("Class orbit.matrix_lattice.MATRIX_Lattice, method trackTwissData(...): direction should be x or y.")
		#track twiss 
		eps_length = 0.00001 # 10^-6 meter
		dir_ind = 0
		if(direction.lower() == "y"):
			dir_ind = 2
		gamma = (1.0+alpha*alpha)/beta
		track_m = Matrix(3,3)
		track_m.unit()
		track_v = PhaseVector(3)
		track_v.set(0,alpha)
		track_v.set(1,beta)
		track_v.set(2,gamma)
		position = 0.
		pos_arr = []
		alpha_arr = []
		beta_arr = []
		#phi is for tune accumulation
		phi = 0.
		#put in array the initial twiss
		pos_arr.append(position)
		alpha_arr.append(track_v.get(0))
		beta_arr.append(track_v.get(1))
		position_old = position
		beta_old = beta
		#count = 0
		for matrixNode in self.getNodes():
			if(isinstance(matrixNode,BaseMATRIX) == True):
				beta = track_v.get(1)
				if(abs(position_old-position) > eps_length or abs(beta_old - beta) > eps_length):
					pos_arr.append(position)
					alpha_arr.append(track_v.get(0))
					beta_arr.append(track_v.get(1))
				mt = matrixNode.getMatrix()
				ind0 = 0+dir_ind
				ind1 = 1+dir_ind
				#if(count < 5):
				#	print "debug =========count =",count, " position=",position
				#	print "debug a00=",mt.get(ind0,ind0)," a01=",mt.get(ind0,ind1)
				#	print "debug a10=",mt.get(ind1,ind0)," a11=",mt.get(ind1,ind1)
				#	print "debug OLD aplha =",track_v.get(0)," beta=",track_v.get(1)," gamma=",track_v.get(2)
				track_m.set(0,0,mt.get(ind0,ind0)*mt.get(ind1,ind1)+mt.get(ind0,ind1)*mt.get(ind1,ind0))
				track_m.set(0,1,-mt.get(ind0,ind0)*mt.get(ind1,ind0))
				track_m.set(0,2,-mt.get(ind0,ind1)*mt.get(ind1,ind1))
				track_m.set(1,0,-2*mt.get(ind0,ind0)*mt.get(ind0,ind1))
				track_m.set(1,1,mt.get(ind0,ind0)*mt.get(ind0,ind0))
				track_m.set(1,2,mt.get(ind0,ind1)*mt.get(ind0,ind1))
				track_m.set(2,0,-2*mt.get(ind1,ind0)*mt.get(ind1,ind1))
				track_m.set(2,1,mt.get(ind1,ind0)*mt.get(ind1,ind0))
				track_m.set(2,2,mt.get(ind1,ind1)*mt.get(ind1,ind1))
				alpha_0 = track_v.get(0)
				beta_0 = track_v.get(1)
				delta_phi = math.atan(	mt.get(ind0,ind1)/(	beta_0*mt.get(ind0,ind0) - alpha_0*mt.get(ind0,ind1)))
				phi = phi + delta_phi
				track_v = track_m.mult(track_v)	
				#if(count < 5):
				#	print "debug NEW aplha =",track_v.get(0)," beta=",track_v.get(1)," gamma=",track_v.get(2)			
				position_old = position
				beta_old = beta_0
				position = position + matrixNode.getLength()
				#count = count + 1
		pos_arr.append(position)
		alpha_arr.append(track_v.get(0))
		beta_arr.append(track_v.get(1))
		#pack the resulting tuple
		tune = phi/(2*math.pi)	
		graph_alpha_arr = []
		graph_beta_arr = []
		for i in range(len(pos_arr)):
			graph_alpha_arr.append((pos_arr[i],alpha_arr[i]))
			graph_beta_arr.append((pos_arr[i],beta_arr[i]))
		return (tune,graph_alpha_arr,graph_beta_arr)

	def getRingDispersionDataX(self,momentum,mass):
		"""
		Returns the tuple  ([(position, dispX),...],[(position,disp_pX),...] ). 
		"""
		res_dict = self.getRingParametersDict(momentum,mass)
		disp = res_dict["dispersion x [m]"]
		disp_p = res_dict["dispersion prime x"]
		return self.trackDispersionData(momentum,mass,disp, disp_p,"x")

	def getRingDispersionDataY(self,momentum,mass):
		"""
		Returns the tuple  ([(position, dispY),...],[(position,disp_pY),...] ). 
		"""
		res_dict = self.getRingParametersDict(momentum,mass)
		disp = res_dict["dispersion y [m]"]
		disp_p = res_dict["dispersion prime y"]
		return self.trackDispersionData(momentum,mass,disp, disp_p,"y")

	def trackDispersionData(self,momentum,mass, disp, disp_p, direction = "x"):
		"""
		Returns the tuple ([(position, disp),...],[(position,disp_p),...] ). 
		The tracking starts from the values specified as the initial parameters. 
		The possible values for direction parameter "x" or "y".
		"""
		if(direction.lower() != "x" and direction.lower() != "y"):
			orbitFinalize("Class orbit.matrix_lattice.MATRIX_Lattice, method trackDispersionData(...): direction should be x or y.")
		#track dispersion
		eps_length = 0.00001 # 10^-6 meter
		dir_ind = 0
		if(direction.lower() == "y"):
			dir_ind = 2
		track_m = Matrix(3,3)
		track_m.unit()
		track_m.set(2,0,0.)
		track_m.set(2,1,0.)			
		track_m.set(2,2,1.)			
		track_v = PhaseVector(3)
		#kinematics coefficient calculation
		Etotal = math.sqrt(momentum**2 + mass**2)
		beta = momentum/Etotal
		gamma = Etotal/mass
		Ekin = Etotal - mass		
		m_coeff =  momentum*momentum/(mass + Ekin)		
		track_v.set(0,disp)
		track_v.set(1,disp_p)
		track_v.set(2,1.)
		position = 0.
		pos_arr = []
		disp_arr = []
		disp_p_arr = []
		#put in array the initial dispersions
		pos_arr.append(position)
		disp_arr.append(track_v.get(0))
		disp_p_arr.append(track_v.get(1))
		position_old = position
		disp_old = disp
		#count = 0
		for matrixNode in self.getNodes():
			if(isinstance(matrixNode,BaseMATRIX) == True):
				disp = track_v.get(0)
				if(abs(position_old-position) > eps_length or abs(disp_old - disp) > eps_length):
					pos_arr.append(position)
					disp_arr.append(track_v.get(0))
					disp_p_arr.append(track_v.get(1))
				disp_old = disp
				position_old = position
				mt = matrixNode.getMatrix()
				ind0 = 0+dir_ind
				ind1 = 1+dir_ind
				#if(count < 200):
				#	print "debug =========count =",count, " position=",position
				#	print "debug a00=",mt.get(ind0,ind0)," a01=",mt.get(ind0,ind1)," a02=",mt.get(ind0,5)
				#	print "debug a10=",mt.get(ind1,ind0)," a11=",mt.get(ind1,ind1)," a12=",mt.get(ind0,5)
				#	print "debug OLD disp =",track_v.get(0)," disp_prime=",track_v.get(1)
				track_m.set(0,0,mt.get(ind0,ind0))
				track_m.set(0,1,mt.get(ind0,ind1))
				track_m.set(1,0,mt.get(ind1,ind0))
				track_m.set(1,1,mt.get(ind1,ind1))		
				track_m.set(0,2,mt.get(ind0,5)*m_coeff)
				track_m.set(1,2,mt.get(ind1,5)*m_coeff)
				track_v = track_m.mult(track_v)	
				#if(count < 200):
				#		print "debug NEW disp =",track_v.get(0)," disp_prime=",track_v.get(1)
				position = position + matrixNode.getLength()
				#count = count + 1
		pos_arr.append(position)
		disp_arr.append(track_v.get(0))
		disp_p_arr.append(track_v.get(1))
		#pack the resulting tuple
		graph_disp_arr = []
		graph_disp_p_arr = []
		for i in range(len(pos_arr)):
			graph_disp_arr.append((pos_arr[i],disp_arr[i]))
			graph_disp_p_arr.append((pos_arr[i],disp_p_arr[i]))
		return (graph_disp_arr,graph_disp_p_arr)
	
	def getSubLattice(self, index_start = -1, index_stop = -1,):
		"""
		It returns the new MATRIX_Lattice with children with indexes 
		between index_start and index_stop inclusive
		"""
		return self._getSubLattice(MATRIX_Lattice(),index_start,index_stop)

	def trackBunch(self, bunch, paramsDict = {}, actionContainer = None):
		"""
		It tracks the bunch through the lattice.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Bunch Tracking")
		paramsDict["bunch"] = bunch
		
		def track(paramsDict):
			node = paramsDict["node"]
			node.track(paramsDict)
			
		actionContainer.addAction(track, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict)
		actionContainer.removeAction(track, AccActionsContainer.BODY)


