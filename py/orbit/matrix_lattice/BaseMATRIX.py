"""
The AccNode subclass that represents a transport 7x7 matrix in a MATRIX_Lattice instance. 
"""

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer

# import matrix class
from orbit_utils import Matrix

class BaseMATRIX(AccNode):
	""" The base abstract class of the BaseMATRIX accelerator lattice elements hierarchy. """
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
		It is tracking the parameter dictionary (with bunch as a "bunch") through the BaseMATRIX element.
		"""
		bunch = paramsDict["bunch"]	
		self.matrix.track(bunch)

