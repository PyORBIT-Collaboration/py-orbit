"""
The module includes classes for an acceleartor lattice.
"""

import sys
import os

from pyORBIT_utils import orbitFinalize

class AccActionsConatainer:
	""" The class of the container of accelerator actions. """
	def __init__(self, name = "no name container" ):
		"""
		Creates an empty accelerator actions container.
		"""
		self.__name = name
		self.__type = "action container"
		self.__entranceActions = []
		self.__bodyActions = []
		self.__exitActions = []
		self.__shouldStop = False	

	def performEntranceActions(self, paramDict):
		"""
		Method. It performs the set of actions at the entrance of accelerator node.
		"""
		for action in self.__entranceActions:
			if(self.__shouldStop):
				return
			action(paramDict)

	def performExitActions(self, paramDict):
		"""
		Method. It performs the set of actions 
		at the exit of accelerator node.
		"""
		for action in self.__exitActions:
			if(self.__shouldStop):
				return
			action(paramDict)

	def performActions(self, paramDict):
		"""
		Method. It performs the set of actions for the 
		body of accelerator node.
		"""
		for action in self.__bodyActions:
			if(self.__shouldStop):
				return
			action(paramDict)
	
	def addEntranceAction(self, action):
		"""
		Method. It adds an entarnce action to 
		the beginning of the container.
		"""
		if(self.__entranceActions.count(action) == 0):
			self.__entranceActions.insert(0,action)

	def getNumbOfEntranceActions(self):
		"""
		Method. It returns the number of entrance actions 
		in the container.
		"""
		return len(self.__entranceActions)

	def getEntranceActions(self):
		"""
		Method. It returns the entrance actions as a list.
		"""
		return self.__entranceActions

	def removeEntranceAction(self,action):
		"""
		Method. It removes the action from the container.
		"""
		return self.__entranceActions.remove(action)

	def addExitAction(self, action):
		"""
		Method. It adds an exit action to 
		the beginning of the container.
		"""
		if(self.__exitActions.count(action) == 0):
			self.__exitActions.insert(0,action)

	def getNumbOfExitActions(self):
		"""
		Method. It returns the number of exit actions 
		in the container.
		"""
		return len(self.__exitActions)

	def getExitActions(self):
		"""
		Method. It returns the exit actions as a list.
		"""
		return self.__exitActions

	def removeExitAction(self,action):
		"""
		Method. It removes the action from the container.
		"""
		return self.__exitActions.remove(action)

	def addAction(self, action):
		"""
		Method. It adds an action to the beginning 
		of the container.
		"""
		if(self.__bodyActions.count(action) == 0):
			self.__bodyActions.insert(0,action)

	def getNumbOfActions(self):
		"""
		Method. It returns the number of  actions 
		in the container.
		"""
		return len(self.__bodyActions)

	def getActions(self):
		"""
		Method. It returns the actions as a list.
		"""
		return self.__bodyActions

	def removeAction(self,action):
		"""
		Method. It removes the action from the container.
		"""
		return self.__bodyActions.remove(action)

	def setShouldStop(self, shouldStop = True):
		"""
		Method. It sets a should_stop property of 
		the actions container.
		"""
		self.__shouldStop = shouldStop 

	def getShouldStop(self):
		"""
		Method. It returns a should_stop property 
		of the actions container.
		"""
		return self.__shouldStop

	def setName(self,name = "no_name"):
		"""
		Method. It sets a name of the actions container.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. It returns a name of the actions container.
		"""
		return self.__name

class AccElement:
	""" The base class of the accelerator elements hierarchy. """
	def __init__(self, name = "no name"):
		"""
		nParts - number of parts to which an elements was devided.
		length the physical length of the elements in m.
		"""
		#------------------------------------------------
		# there is no position parameter, 
		# because element could be in different lattices 
		#------------------------------------------------
		self.__name = name
		self.__type = "generic"
		self.__nParts = 1
		self.__lengthArr = [0.]
		self.__length = 0.
		self.__activePartIndex = 0
		self.__isInitialized = False
		self.__params = {}
		# the elements that should be used at the start,
		#     inside, and at the finish of this element
		#     __insideChildNodesDic - dictinary with indexes
		#             as keys and lists of nodes as elements
		self.__startChildNodes = []
		self.__insideChildNodesDic = {0:[]}
		self.__finishChildNodes = []

	def initialize(self, paramDict):
		"""
		Method. It initializes the necessary structures of the node.
		This method will call  _initialize(probe) method of the subclass before all
		children initialization. The subclasses should implement __initialize(probe)
		method.
		"""
		#----initialization self---------
		self._initialize(paramDict)
		self.setInitialized(True)

	def _initialize(self, paramDict):
		"""
		Abstract Method. The subclasses should implement _initialize(probe) method.
		"""
		pass

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. It is tracking the actions through the accelerator element.
		"""
		paramsDict["node"] = self
		parentNode = paramsDict["parentNode"]
		self.__activePartIndex = 0

		if(actionsContainer.getShouldStop()):
			return		

		actionsContainer.performEntranceActions(paramsDict)

		if(actionsContainer.getShouldStop()):
			return		

		for node in self.__startChildNodes:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self
			
			if(actionsContainer.getShouldStop()):
				return		

			node.trackActions(actionsContainer, paramsDict)

		for i in xrange(self.__nParts):
			paramsDict["node"] = self
			paramsDict["parentNode"] = parentNode
			self.__activePartIndex = i
			
			if(actionsContainer.getShouldStop()):
				return					
			
			actionsContainer.performActions(paramsDict)

			nodes = self.__insideChildNodesDic[i]
			for node in nodes:
					paramsDict["node"] = node
					paramsDict["parentNode"] = self

					if(actionsContainer.getShouldStop()):
						return							
					
					node.trackActions(actionsContainer, paramsDict)

		for node in self.__finishChildNodes:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self

			if(actionsContainer.getShouldStop()):
				return					
			
			node.trackActions(actionsContainer, paramsDict)

		paramsDict["node"] = self
		paramsDict["parentNode"] = 	parentNode

		if(actionsContainer.getShouldStop()):
			return				
		
		actionsContainer.performExitActions(paramsDict)
		self.__activePartIndex = 0

	def track(self, paramDict):
		"""
		Abstract Method. It is tracking the dictionary with parameters through
		the node. The subclasses should implement this method.
		"""
		pass

	def addChildNodeAtStart(self, node, index = -1):
		"""
		Method. It adds the child node to the beginning of the this node.
		The third parameter is an index of node among the child nodes
		at the start of the node.
		"""
		if(isinstance(node,AccElement) != True):
			msg = "The child of AccElemen can be its subclass only!"
			msg = msg + "method addChildNodeAtStart(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		if(index <= 0):
			self.__startChildNodes.insert(0,node)
		else:
			if(index < len(self.__startChildNodes)):
				self.__startChildNodes.insert(index,node)
			else:
				self.__startChildNodes.append(node)

	def addChildNodeAtFinish(self, node, index = -1):
		"""
		Method. It adds the child node to the end of the this node.
		The third parameter is an index of node among the child nodes
		at the end of the node.
		"""
		if(isinstance(node,AccElement) != True):
			msg = "The child of AccElemen can be its subclass only!"
			msg = msg + "method addChildNodeAtFinish(self, node):"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		if(index <= 0):
			self.__finishChildNodes.insert(0,node)
		else:
			if(index < len(self.__finishChildNodes)):
				self.__finishChildNodes.insert(index,node)
			else:
				self.__finishChildNodes.append(node)

	def addChildNode(self, node, index = 0):
		"""
		Method. It adds the child node to the sub-element
		with a certain index.
		"""
		if(isinstance(node,AccElement) != True):
			msg = "The child of AccElemen can be its subclass only!"
			msg = msg + "method addChildNode(self, node, index = 0)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		chArr = self.__insideChildNodesDic[index]
		chArr.insert(0,node)

	def getChildrenNodesAll(self):
		"""
		Method. It returns all children of this node as a list.
		"""
		res = []
		for node in self.__startChildNodes:
			res.append(node)
		for i in xrange(self.__nParts):
			nodes = self.__insideChildNodesDic[i]
			for node in nodes:
				res.append(node)
		for node in self.__finishChildNodes:
			res.append(node)
		return res

	def getChildrenNodesAtStart(self):
		"""
		Method. It returns all children of this node at START of the node.
		"""
		res = []
		for node in self.__startChildNodes:
			res.append(node)
		return res

	def getChildrenNodesAtFinish(self):
		"""
		Method. It returns all children of this node at FINISH of the node.
		"""
		res = []
		for node in self.__finishChildNodes:
			res.append(node)
		return res
	
	def getChildrenNodesInMiddle(self):
		"""
		Method. It returns all children of this node that are in the middle.
		"""
		res = []
		for i in xrange(self.__nParts):
			nodes = self.__insideChildNodesDic[i]
			for node in nodes:
				res.append(node)
		return res	

	def setName(self,name = "no_name"):
		"""
		Method. It sets a name of the node.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. It returns a name of the node.
		"""
		return self.__name

	def setType(self, tp = "generic"):
		"""
		Method. It sets a type of the node.
		"""
		self.__type = tp

	def getType(self):
		"""
		Method. It returns a type of the node.
		"""
		return self.__type

	def setnParts(self, n = 1):
		"""
		Method. It sets a number of inner parts of the node. 
		This method does not call the initilize() method.
		"""
		self.__nParts = n
		self.__lengthArr = []
		self.__insideChildNodesDic = {}		
		for i in xrange(self.__nParts):
			self.__lengthArr.append(0.)
			self.__insideChildNodesDic[i] = []

	def getnParts(self):
		"""
		Method. It returns a number of inner parts of the node.
		"""
		return self.__nParts

	def getActivePartIndex(self):
		"""
		Method. It returns an active part index of the node.
		"""
		return self.__activePartIndex

	def setLength(self, L = 0.):
		"""
		Method. It sets a physical length of the node.
		"""
		if(abs(L) < 1.0e-9): L = 0.
		self.__length = L

	def getLength(self):
		"""
		Method. It returns a physical length of the node.
		"""
		return self.__length

	def setPartLength(self, index, L = 0.):
		"""
		Method. It sets a physical length of the part
		of the node with particular index.
		"""
		self.__lengthArr[index] = L

	def getPartLength(self, index):
		"""
		Method. It returns a physical length of the part
		of the node with particular index.
		"""
		return self.__lengthArr[index]

	def addParam(self, key, value):
		"""
		Sets the parameter of the node
		"""
		self.__params[key] = value

	def addParams(self, params):
		"""
		Adds the parameters from the external 
		dictionary to the parameters of the node.
		The existing keys' values will be replaced 
		by the new ones from the external dictionary.
		"""		
		self.__params.update(params)
		
	def setParams(self, params):
		"""
		Sets the external dictionary as a parameters dictionary for this element.
		"""
		self.__params = params
		
	def getParam(self, key):
		"""
		Returns the parameter of the node
		"""
		if(not self.hasParam(key)):
			msg = "The node does not have a prameter for the key you asked for!"
			msg = msg + "method getParam(self, key)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "key=" + str(key)
			orbitFinalize(msg)
		return self.__params[key]

	def hasParam(self, key):
		"""
		Returns true if the node has a parameter for this key.
		"""
		return self.__params.has_key(key)

	def getParams(self):
		"""
		Returns the whole parameters dictionary.
		"""
		return self.__params

	def isInitialized(self):
		"""
		Method. It returns true or false about an initialization status.
		"""
		return self.__isInitialized

	def setInitialized(self, status):
		"""
		Method. It sets true or false about an initialization status.
		"""
		if(type(status) != type(True)):
			orbitFinalize("setInitialized(boolean) needs boolean")
		self.__isInitialized = status


class AccLine:
	"""
	The accelerator line class.
	A Line could include other lines or elements.
	"""
	def __init__(self, name = "no name"):
		"""
		A constructor for an accelerator line class instance.
		"""
		#------------------------------------------------
		# there is no position parameter, 
		# because AccLine could be in different lattices 
		#------------------------------------------------		
		self.__name = name
		self.__type = "line"
		self.__length = 0.
		self.__isInitialized = False
		self.__children = []

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. It is tracking the actions through the accelerator line.
		"""
		paramsDict["node"] = self
		parentNode = paramsDict["parentNode"]

		if(actionsContainer.getShouldStop()):
			return		

		actionsContainer.performEntranceActions(paramsDict)

		for node in self.__children:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self

			if(actionsContainer.getShouldStop()):
				return					
			
			node.trackActions(actionsContainer, paramsDict)

		paramsDict["node"] = self
		paramsDict["parentNode"] = 	parentNode

		if(actionsContainer.getShouldStop()):
			return		
		actionsContainer.performExitActions(paramsDict)


	def addChildNode(self, node):
		"""
		Method. It adds a child node at the end of the line.
		"""
		if(isinstance(node,AccElement) != True and isinstance(node,AccLine) != True ):
			msg = "The child of AccLine can be AccElement or AccLine) only!"
			msg = msg + "method addChildNode(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		self.__children.append(node)

	def insertChildNode(self, node, index = 0):
		"""
		Method. It inserts a child node at the position with certain index.
		"""
		if(isinstance(node,AccElement) != True and isinstance(node,AccLine) != True ):
			msg = "The child of AccLine can be AccElement or AccLine) only!"
			msg = msg + "method insertChildNode(self, node, index)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		self.__children.insert(index,node)

	def removeChildNode(self, index = 0):
		"""
		Method. It removes a child node at the position with certain index.
		"""
		if(index > (len(self.__children)-1)):
			msg = "You try to remove the child node with wrong index!"
			msg = msg + "method removeChildNode(self, index = 0)"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "index=" + str(index)
			msg = msg + os.linesep
			msg = msg + "number of children=" + str(self.getChildrenNodesCount())
			orbitFinalize(msg)
		del self.__children[index:(index+1)]

	def removeAllChildren(self):
		"""
		Method. It removes all children nodes.
		"""
		self.__children = []

	def getChildrenNodes(self):
		"""
		Method. It returns all children of this node as a list.
		"""
		return self.__children

	def setName(self,name = "no_name"):
		"""
		Method. It sets a name of the line.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. It returns a name of the line.
		"""
		return self.__name

	def getType(self):
		"""
		Method. It returns a type of the line.
		"""
		return self.__type

	def setLength(self, L = 0.):
		"""
		Method. It sets a physical length of the line.
		"""
		if(abs(L) < 1.0e-9): L = 0.
		self.__length = L

	def getLength(self):
		"""
		Method. It returns a physical length of the line.
		"""
		return self.__length

	def isInitialized(self):
		"""
		Method. It returns true or false about an initialization status.
		"""
		return self.__isInitialized

	def setInitialized(self, initialized = True):
		"""
		Method. Sets an initialization status.
		"""
		self.__isInitialized = initialized

class AccLattice:
	""" The class of the accelerator lattice. """
	def __init__(self,name = "no name"):
		"""
		Method. It creates an empty accelerator lattice.
		"""
		self.__name = name
		self.__type = "lattice"
		self.__length = 0.
		self.__isInitialized = False
		self.__children = []

	def initialize(self, actions = None, paramsDict = {}):
		"""
		Method. It initializes the necessary structures of the lattice and child nodes.
		"""
		if(actions == None):
			actions = AccActionsConatainer()
		d = {"position":0,"position_start_line":0}
		d["position"] = 0.
		d["position_start_line"] = []

		def accElemEntrance(paramsDict):
			node = paramsDict["node"]
			if(isinstance(node,AccElement)):
				node.initialize(paramsDict)
			if(isinstance(node,AccLine)):
				node.setLength(0.)
				d["position_start_line"].append(d["position"])

		def accElemExit(paramsDict):
			node = paramsDict["node"]
			d["position"] += node. getLength()
			if(isinstance(node,AccLine)):
				line_ind = len(d["position_start_line"]) - 1
				start_pos = d["position_start_line"][line_ind]
				del d["position_start_line"][line_ind]
				node.setLength(d["position"] - start_pos)

		actions.addEntranceAction(accElemEntrance)
		actions.addExitAction(accElemExit)

		self.trackActions(actions,paramsDict)
		self.__length = d["position"]
		self.__isInitialized = True

	def trackActions(self, actionsContainer, paramsDict = {}):
		"""
		Method. It is tracking the actions through the all nodes.
		"""
		paramsDict["lattice"] = self
		paramsDict["actions"] = actionsContainer
		paramsDict["node"] = self
		paramsDict["parentNode"] = None

		actionsContainer.performEntranceActions(paramsDict)
		
		for node in self.__children:
			paramsDict["node"] = node
			paramsDict["parentNode"] = self

			if(actionsContainer.getShouldStop()):
				return	
			
			node.trackActions(actionsContainer, paramsDict)

			if(actionsContainer.getShouldStop()):
				return			

		paramsDict["node"] = self
		paramsDict["parentNode"] = None
		actionsContainer.performExitActions(paramsDict)

	def addChildNode(self, node):
		"""
		Method. It adds the child node to the end of the this node.
		"""
		if(isinstance(node,AccElement) != True and isinstance(node,AccLine) != True):
			msg = "The child of AccLattice can be  AccElement or AccLine only!"
			msg = msg + "method addChildNode(self, node)"
			msg = msg + os.linesep
			msg = msg + "Name of lattice=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "Child node=" + str(node)
			orbitFinalize(msg)
		self.__children.append(node)

	def setName(self,name = "no_name"):
		"""
		Method. It sets a name of the lattice.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. It returns a name of the lattice.
		"""
		return self.__name

	def getType(self):
		"""
		Method. It returns a type of the lattice.
		"""
		return self.__type

	def getLength(self):
		"""
		Method. It returns a physical length of the node.
		"""
		return self.__length

	def isInitialized(self):
		"""
		Method. It returns true or false about an initialization status.
		"""
		return self.__isInitialized

	def setInitialized(self, initialized = True):
		"""
		Method. Sets an initialization status.
		"""
		actions = AccActionsConatainer()

		def accElemExit(paramsDict):
			node = paramsDict["node"]
			node.setInitialized(initialized)

		actions.addExitAction(accElemExit)
		self.trackActions(actions,paramsDict)		
		self.__isInitialized = initialized

	def trackBunch(self, bunch, actions = None, paramsDict = {}):
		"""
		Method. It track the bunch trough the lattice.
		"""
		paramsDict["bunch"] = bunch
		if(actions == None):
			actions = AccActionsConatainer()

		def track(paramsDict):
			node = paramsDict["node"]
			if(isinstance(node,AccElement)):
				node.track(paramsDict)

		actions.addAction(track)
		self.trackActions(actions,paramsDict)