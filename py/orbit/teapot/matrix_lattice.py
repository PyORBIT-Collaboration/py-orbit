"""
Module. Includes classes for a Matrix lattice. The Matrix lattice is generated
by using the TEAPOT lattice. The matrices are the linear part of the TEAPOT elements
tracking. The number of transport matrices in the lattice is equal to the sum of 
all parts of TEAPOT elements. The RF cavities in the Matrix lattice are the 
TEAPOT RF Cavity class instances.
"""
import os
import math

#import bunch
from bunch import Bunch

# import the function that creates multidimensional arrays
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer

# import matrix class and generators
from orbit.teapot_base import MatrixGenerator
from orbit_utils import Matrix
from orbit_utils import PhaseVector

# import the MAD parser to construct lattices of TEAPOT elements.
from orbit.teapot import TEAPOT_Lattice, RingRFTEAPOT, BaseTEAPOT

class MATRIX_Lattice(AccLattice):
	"""
	The subclass of the AccLattice class. Shell class for the BaseMATRIX nodes collection.
	The TEAPOT_Lattice instance is needed to create the MATRIX_Lattice instance. The nodes for
	RF elements are the usual TEAPOT RingRFTEAPOT nodes. The Bunch instance is needed for
	the MATRIX_Lattice constructor to specify the particle's energy and other parameters.
	"""
	def __init__(self, teapot_lattice, bunch, name = None):
		AccLattice.__init__(self,name)
		self.oneTurmMatrix = Matrix(7,7)
		self.oneTurmMatrix.unit()		
		if(isinstance(teapot_lattice,TEAPOT_Lattice) != True):
			orbitFinalize("Constructor orbit.teapot.MATRIX_Lattice needs  the TEAPOT_Lattice instance.")
		#memorize the TEAPOT LATTICE
		if(name == None):
			name = teapot_lattice.getName()
		self.setName(name)
		self.teapot_lattice = teapot_lattice
		self.bunch = Bunch()
		bunch.copyEmptyBunchTo(self.bunch)
		self.matrixGenerator = MatrixGenerator()
		#----------make MATRIX lattice from TEAPOT
		self.index = 0
		self.matrixArr = []
		
		def twissAction(paramsDict):
			node = paramsDict["node"]
			bunch = paramsDict["bunch"]
			active_index = node.getActivePartIndex()
			n_parts = node.getnParts()
			length = node.getLength(active_index)
			if(isinstance(node, BaseTEAPOT) == True and isinstance(node,RingRFTEAPOT) == False): 
				self.matrixGenerator.initBunch(bunch)
				node.track(paramsDict)
				#bunch is ready
				matrixNode = BaseMATRIX(node.getName()+"_"+str(active_index))
				matrixNode.addParam("matrix_parent_node",node)
				matrixNode.addParam("matrix_parent_node_type",node.getType())
				matrixNode.addParam("matrix_parent_node_n_nodes",n_parts)
				matrixNode.addParam("matrix_parent_node_active_index",active_index)
				matrixNode.addParam("matrix_index",self.index)
				matrixNode.setLength(length)
				self.matrixGenerator.calculateMatrix(bunch,matrixNode.getMatrix())
				#print "debug i=",self.index," name=",matrixNode.getName(),
				#print " type=",matrixNode.getParam("matrix_parent_node_type"),			
				#print " L=",matrixNode.getLength()	
				self.addNode(matrixNode)
				self.matrixArr.append(matrixNode)
				self.index = self.index + 1				
			if(isinstance(node,RingRFTEAPOT) == True):
				rf_node = RingRFTEAPOT(node.getName())
				rf_node.setParamsDict(node.getParamsDict().copy())
				rf_node.addParam("matrix_index",self.index)
				self.addNode(rf_node)
				self.index = self.index + 1	
				
		accContainer = AccActionsContainer()
		accContainer.addAction(twissAction,AccActionsContainer.BODY)
		paramsDict = {}
		paramsDict["bunch"] = self.bunch
		paramsDict["position"] = 0.
		self.index = 0
		self.teapot_lattice.trackActions(accContainer,paramsDict)		
		self.makeOneTurnMatrix()
		self.initialize()		
					
	def getKinEnergy(self):
		return self.bunch.getSyncParticle().kinEnergy()

	def rebuild(self, Ekin = -1.0):
		if(Ekin > 0.):
			self.bunch.getSyncParticle().kinEnergy(Ekin)
		for matrixNode in self.matrixArr:
			node = matrixNode.getParam("matrix_parent_node")
			active_index = matrixNode.getParam("matrix_parent_node_active_index")
			n_parts = matrixNode.getParam("matrix_parent_node_n_nodes")
			if(n_parts != node.getnParts()):
				msg = " orbit.teapot.MATRIX_Lattice class" + os.linesep
				msg = msg + "  rebuild(Ekin = -1.0) method" + os.linesep
				msg = msg + "  TEAPOT node="+node.getName() + os.linesep
				msg = msg + "  has been changed!" + os.linesep
				msg = msg + "  Stop!" + os.linesep
				orbitFinalize(msg)
			self.matrixGenerator.initBunch(self.bunch)
			paramsDict = {}
			paramsDict["bunch"] = self.bunch
			paramsDict["node"] = node
			node.setActivePartIndex(active_index)
			node.track(paramsDict)
			self.matrixGenerator.calculateMatrix(self.bunch,matrixNode.getMatrix())
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
		
	def getRingParametersDict(self):
		"""
		Returns the dictionary with different ring parametrs
		calculated from the one turn transport matrix.
		"""
		res_dict = {}
		momentum = self.bunch.getSyncParticle().momentum()
		beta = self.bunch.getSyncParticle().beta()
		gamma = self.bunch.getSyncParticle().gamma()
		mass = self.bunch.getSyncParticle().mass()
		Ekin = self.bunch.getSyncParticle().kinEnergy()
		res_dict["momentum [GeV/c]"] = momentum
		res_dict["mass [GeV]"] = mass
		res_dict["Ekin [GeV]"] = Ekin
		#longitudinal params
		c = 2.99792458e+8
		ring_length = self.getLength()
		T = ring_length/(beta*c)
		res_dict["period [sec]"] = T
		res_dict["frequency [Hz]"] = 1./T
		eta_ring = self.oneTurmMatrix	.get(4,5)*beta*beta*(Ekin+mass)/ring_length
		res_dict["eta"] = eta_ring
		gamma_trans =1.0/math.sqrt( eta_ring + 1./(gamma*gamma))
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

	def getRingTwissDataX(self):
		"""
		Returns the tuple (tuneX, [(position, alphaX),...],[(position,betaX),...] ). 
		"""
		res_dict = self.getRingParametersDict()
		alpha_x = res_dict["alpha x"]
		beta_x = res_dict["beta x [m]"]
		return self.trackTwissData(alpha_x,beta_x,"x")
	
	def getRingTwissDataY(self):
		"""
		Returns the tuple (tuneY, [(position, alphaY),...],[(position,betaY),...] ).
		"""
		res_dict = self.getRingParametersDict()
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
			orbitFinalize("Class orbit.teapot.MATRIX_Lattice, method trackTwissData(...): direction should be x or y.")
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
		for matrixNode in self.matrixArr:
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

	def getRingDispersionDataX(self):
		"""
		Returns the tuple  ([(position, dispX),...],[(position,disp_pX),...] ). 
		"""
		res_dict = self.getRingParametersDict()
		disp = res_dict["dispersion x [m]"]
		disp_p = res_dict["dispersion prime x"]
		return self.trackDispersionData(disp, disp_p,"x")

	def getRingDispersionDataY(self):
		"""
		Returns the tuple  ([(position, dispY),...],[(position,disp_pY),...] ). 
		"""
		res_dict = self.getRingParametersDict()
		disp = res_dict["dispersion y [m]"]
		disp_p = res_dict["dispersion prime y"]
		return self.trackDispersionData(disp, disp_p,"y")

	def trackDispersionData(self, disp, disp_p, direction = "x"):
		"""
		Returns the tuple ([(position, disp),...],[(position,disp_p),...] ). 
		The tracking starts from the values specified as the initial parameters. 
		The possible values for direction parameter "x" or "y".
		"""
		if(direction.lower() != "x" and direction.lower() != "y"):
			orbitFinalize("Class orbit.teapot.MATRIX_Lattice, method trackDispersionData(...): direction should be x or y.")
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
		momentum = self.bunch.getSyncParticle().momentum()
		mass = self.bunch.getSyncParticle().mass()
		Ekin = self.bunch.getSyncParticle().kinEnergy()		
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
		for matrixNode in self.matrixArr:
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

	def getChromaticitiesXY(self):
		"""
		Calculates chromaticities for X,Y planes for the whole ring
		"""
		(tuneX,tmp0,tmp1) = self.getRingTwissDataX()
		(tuneY,tmp0,tmp1) = self.getRingTwissDataY()
		self.matrixGenerator.initBunchChromCoeff(self.bunch)
		#track bunch through the TEAPOT nodes
		def twissAction(paramsDict):
			node = paramsDict["node"]
			if(isinstance(node, BaseTEAPOT) == True and isinstance(node,RingRFTEAPOT) == False): 
				node.track(paramsDict)
		
		accContainer = AccActionsContainer()
		accContainer.addAction(twissAction,AccActionsContainer.BODY)
		paramsDict = {}
		paramsDict["bunch"] = self.bunch
		self.teapot_lattice.trackActions(accContainer,paramsDict)	
		
		res_coeff = self.matrixGenerator.calcChromCoeff(self.bunch)
		(coeff_x_dE,coeff_xp_dE,coeff_y_dE,coeff_yp_dE) = res_coeff
		momentum = self.bunch.getSyncParticle().momentum()
		mass = self.bunch.getSyncParticle().mass()
		Ekin = self.bunch.getSyncParticle().kinEnergy()
		chromX = - (momentum/(2*math.sin(2*math.pi*tuneX)))*(coeff_x_dE+coeff_xp_dE)
		chromX = (momentum/(mass+Ekin))*chromX
		chromY = - (momentum/(2*math.sin(2*math.pi*tuneY)))*(coeff_y_dE+coeff_yp_dE)
		chromY = (momentum/(mass+Ekin))*chromY
		return (chromX/(2*math.pi),chromY/(2*math.pi))
		

	def initialize(self):
		"""
		Initializes the lattice.
		"""
		AccLattice.initialize(self)
	
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

class BaseMATRIX(AccNode):
	""" The base abstract class of the MATRIX accelerator elements hierarchy. """
	def __init__(self, name = "no name"):
		"""
		Constructor. Creates the base MATRIX element.
		"""
		AccNode.__init__(self,name)
		self.setType("base matrix")
		self.matrix = Matrix(7,7)
		
	def getMatrix(self):
		"""
		Returns the (7,7) Matrix for this transport AccNode.
		"""
		return self.matrix
		
	def trackBunch(self, bunch, paramsDict = {}, actionContainer = None):
		"""
		It tracks the bunch through the BaseMATRIX instance.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Bunch Tracking")
		paramsDict["bunch"] = bunch
		
		def track(paramsDict):
			node = paramsDict["node"]
			node.track(paramsDict)
			
		actionContainer.addAction(track, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict)
		actionContainer.removeAction(track, AccActionsContainer.BODY)		
		
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the element.
		"""
		bunch = paramsDict["bunch"]	
		self.matrix.track(bunch)

