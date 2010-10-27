from orbit.lattice import AccActionsContainer, AccNode

class AccNodeBunchTracker(AccNode):
	"""
	Class. Base class of the accelerator nodes that track the pyORBIT bunch.
	"""

	def __init__(self, name = "no name", type_in = "bunch tracker"):
		"""
		Constructor. Creates an empty bunch tracker accelerator node.
		"""
		AccNode.__init__(self, name, type_in)
		
	def trackBunch(self, bunch, paramsDict = {}, actionContainer = None):
		"""
		It tracks the bunch through the AccNodeBunchTracker instance.
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
		It is tracking the bunch through the element. Each element 
		should implement this method.
		"""
		pass	
