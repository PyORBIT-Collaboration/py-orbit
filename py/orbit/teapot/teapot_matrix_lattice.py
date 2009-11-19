"""
The TEAPOT MATRIX_Lattice is a subclass of a MAXTRIX_Lattice class. The Matrix lattice is generated
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

# import C++ matrix, phase vector, and generator classes
from orbit.teapot_base import MatrixGenerator

from orbit.matrix_lattice import MATRIX_Lattice, BaseMATRIX

# import the MAD parser to construct lattices of TEAPOT elements.
from orbit.teapot import TEAPOT_Lattice, RingRFTEAPOT, BaseTEAPOT

class TEAPOT_MATRIX_Lattice(MATRIX_Lattice):
	"""
	The subclass of the MATRIX_Lattice class. Shell class for the BaseMATRIX nodes collection.
	The TEAPOT_Lattice instance is needed to create the TEAPOT_MATRIX_Lattice instance. 
	The nodes for RF elements are the usual TEAPOT RingRFTEAPOT nodes. The Bunch instance 
	is needed for the MATRIX_Lattice constructor to specify the particle's energy and 
	other parameters.
	"""
	def __init__(self, teapot_lattice, bunch, name = None):
		MATRIX_Lattice.__init__(self,name)
		if(isinstance(teapot_lattice,TEAPOT_Lattice) != True):
			orbitFinalize("Constructor orbit.teapot.TEAPOT_MATRIX_Lattice needs  the TEAPOT_Lattice instance.")
		#memorize the TEAPOT LATTICE
		if(name == None):
			name = teapot_lattice.getName()
		self.setName(name)
		self.teapot_lattice = teapot_lattice
		self.bunch = Bunch()
		bunch.copyEmptyBunchTo(self.bunch)
		self.matrixGenerator = MatrixGenerator()
		#----------make MATRIX lattice from TEAPOT		
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
				matrixNode.setLength(length)
				self.matrixGenerator.calculateMatrix(bunch,matrixNode.getMatrix())
				#print "============= name=",matrixNode.getName(),
				#print " type=",matrixNode.getParam("matrix_parent_node_type"),			
				#print " L=",matrixNode.getLength()	
				self.addNode(matrixNode)		
			if(isinstance(node,RingRFTEAPOT) == True):
				rf_node = RingRFTEAPOT(node.getName())
				rf_node.setParamsDict(node.getParamsDict().copy())
				self.addNode(rf_node)
				
		accContainer = AccActionsContainer()
		accContainer.addAction(twissAction,AccActionsContainer.BODY)
		paramsDict = {}
		paramsDict["bunch"] = self.bunch
		paramsDict["position"] = 0.
		self.teapot_lattice.trackActions(accContainer,paramsDict)		
		self.makeOneTurnMatrix()
		self.initialize()		
					
	def getKinEnergy(self):
		return self.bunch.getSyncParticle().kinEnergy()

	def rebuild(self, Ekin = -1.0):
		if(Ekin > 0.):
			self.bunch.getSyncParticle().kinEnergy(Ekin)
		for matrixNode in self.getNodes():
			if(isinstance(matrixNode,BaseMATRIX) == True):
				node = matrixNode.getParam("matrix_parent_node")
				active_index = matrixNode.getParam("matrix_parent_node_active_index")
				n_parts = matrixNode.getParam("matrix_parent_node_n_nodes")
				if(n_parts != node.getnParts()):
					msg = " orbit.teapot.TEAPOT_MATRIX_Lattice class" + os.linesep
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
				
	def getRingParametersDict(self):
		"""
		Returns the dictionary with different ring parametrs
		calculated from the one turn transport matrix. It overloads the 
		getRingParametersDict(p,m) method from the parent MATRIX_Lattice 
		class.
		"""
		momentum = self.bunch.getSyncParticle().momentum()
		mass = self.bunch.getSyncParticle().mass()	
		return MATRIX_Lattice.getRingParametersDict(self, momentum, mass)

	def getRingTwissDataX(self):
		"""
		Returns the tuple (tuneX, [(position, alphaX),...],[(position,betaX),...] ). 
		It overloads the getRingTwissDataX(p,m) method from the parent MATRIX_Lattice 
		class.
		"""
		res_dict = self.getRingParametersDict()
		alpha_x = res_dict["alpha x"]
		beta_x = res_dict["beta x [m]"]
		return self.trackTwissData(alpha_x,beta_x,"x")
	
	def getRingTwissDataY(self):
		"""
		Returns the tuple (tuneY, [(position, alphaY),...],[(position,betaY),...] ).
		It overloads the getRingTwissDataY(p,m) method from the parent MATRIX_Lattice 
		class.
		"""
		res_dict = self.getRingParametersDict()
		alpha_y = res_dict["alpha y"]
		beta_y = res_dict["beta y [m]"]
		return self.trackTwissData(alpha_y,beta_y,"y")

	def getRingDispersionDataX(self):
		"""
		Returns the tuple  ([(position, dispX),...],[(position,disp_pX),...] ). 
		It overloads the getRingDispersionDataX(p,m) method from the parent MATRIX_Lattice 
		class.
		"""
		res_dict = self.getRingParametersDict()
		disp = res_dict["dispersion x [m]"]
		disp_p = res_dict["dispersion prime x"]
		momentum = res_dict["momentum [GeV/c]"]
		mass = res_dict["mass [GeV]"]			
		return self.trackDispersionData(momentum, mass, disp, disp_p,"x")

	def getRingDispersionDataY(self):
		"""
		Returns the tuple  ([(position, dispY),...],[(position,disp_pY),...] ). 
		It overloads the getRingDispersionDataY(p,m) method from the parent MATRIX_Lattice 
		class.
		"""
		res_dict = self.getRingParametersDict()
		disp = res_dict["dispersion y [m]"]
		disp_p = res_dict["dispersion prime y"]
		momentum = res_dict["momentum [GeV/c]"]
		mass = res_dict["mass [GeV]"]		
		return self.trackDispersionData(momentum, mass, disp, disp_p,"y")

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

