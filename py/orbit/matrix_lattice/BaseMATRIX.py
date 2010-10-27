"""
The AccNode subclass that represents a transport 7x7 matrix in a MATRIX_Lattice instance. 
"""

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import matrix class
from orbit_utils import Matrix

class BaseMATRIX(AccNodeBunchTracker):
	""" The base abstract class of the BaseMATRIX accelerator lattice elements hierarchy. """
	def __init__(self, name = "no name"):
		"""
		Constructor. Creates the base MATRIX element.
		"""
		AccNodeBunchTracker.__init__(self,name)
		self.setType("base matrix")
		self.matrix = Matrix(7,7)
		
	def getMatrix(self):
		"""
		Returns the (7,7) Matrix for this transport AccNode.
		"""
		return self.matrix
				
	def track(self, paramsDict):
		"""
		It is tracking the parameter dictionary (with bunch as a "bunch") through the BaseMATRIX element.
		"""
		bunch = paramsDict["bunch"]	
		self.matrix.track(bunch)

