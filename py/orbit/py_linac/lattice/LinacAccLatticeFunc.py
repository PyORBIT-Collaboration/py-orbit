"""
This package is a collection of the lattice functions.
Some of these functions could be lattice class methods,
but I think they are too specific. So they will be functions.
"""

import os
import math
import sys

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from orbit.py_linac.lattice import BaseLinacNode, Drift, Quad
from orbit.py_linac.lattice import LinacMagnetNode
from orbit.py_linac.lattice import AxisFieldRF_Gap
from orbit.py_linac.lattice import AxisField_and_Quad_RF_Gap

# import acc. nodes
from LinacAccNodes import Quad
from LinacRfGapNodes import AxisFieldRF_Gap

from LinacFieldOverlappingNodes import AxisField_and_Quad_RF_Gap
from LinacFieldOverlappingNodes import OverlappingQuadsNode

def GetGlobalQuadGradient(accLattice,z):
	"""
	The service function for the overlapping fields package.
	Returns the quad field for certain position in the lattice that
	has usual and overlapping quads.
	"""
	node_pos_dict = accLattice.getNodePositionsDict()
	nodes = accLattice.getNodes()
	G = 0.
	(node,index,posBefore,posAfter) = accLattice.getNodeForPosition(z)
	if(isinstance(node,Quad)):
		return node.getParam("dB/dr")
	if(isinstance(node,OverlappingQuadsNode)):
		G = node.getTotalField(z - (posBefore+posAfter)/2)			
		return G
	if(isinstance(node,AxisField_and_Quad_RF_Gap)):
		(z_min,z_max) = node.getZ_Min_Max()
		G = node.getTotalField((z - posBefore)+z_min)			
		return G		
	return G

def GetGlobalRF_AxisField(accLattice,z):
	"""
	The service function for the overlapping RF fields package.
	Returns the RF field on the axis of the RF cavities 
	for certain position in the lattice. If we have 
	the BaseRF_Gap instance we will get 0, because it is 
	an element with zero length.
	"""	
	node_pos_dict = accLattice.getNodePositionsDict()
	nodes =accLattice.getNodes()
	Ez = 0.
	(node,index,posBefore,posAfter) = accLattice.getNodeForPosition(z)
	if(isinstance(node,AxisField_and_Quad_RF_Gap) or isinstance(node,AxisFieldRF_Gap)):
		(z_min,z_max) = node.getZ_Min_Max()
		Ez = node.getEzFiled(z - posBefore + z_min)
		modePhase = node.getParam("mode")*math.pi
		Ez = Ez*math.cos(modePhase)
	return Ez
	
def getNodeForNameFromWholeLattice(accLattice,name):
	"""
	Returns the accelerator node or an array of nodes with the same name.
	This function could be replaced later by the method in the AccLattice class.
	"""
	
	paramsDict = {}
	actions = AccActionsContainer()
	nodes = []

	def accNodeExitAction(paramsDict):
		"""
		Non-bound function. Finds the node in the lattice
		with the specified name.
		positions. This is a closure (well, maybe not
		exactly). It uses external objects.
		"""
		node = paramsDict["node"]
		if(node.getName() == name):
			nodes.append(node)
		
	actions.addAction(accNodeExitAction, AccNode.EXIT)
	accLattice.trackActions(actions, paramsDict)
	if(len(nodes) == 1):
		return nodes[0]
	elif(len(nodes) == 0):
		return None
	else:
		return nodes
		
def getNodePosDictForWholeLattice(accLattice):
	"""
	Returns the dict[node] = (posStart,posEnd) for all nodes (not only for the firts level).
	This function could be replaced later by the method in the AccLattice class.
	"""
	paramsDict = {}
	actions = AccActionsContainer()
	posStartDict = {}
	posStopDict = {}

	def accNodeEntranceAction(paramsDict):
		"""
		Non-bound function. Sets node's end positions. 
		This is a closure (well, maybe not exactly). .
		"""
		node = paramsDict["node"]
		pos = paramsDict["path_length"]
		posStartDict[node] = pos

	def accNodeExitAction(paramsDict):
		"""
		Non-bound function. Sets node's end positions. 
		This is a closure (well, maybe not exactly). .
		"""
		node = paramsDict["node"]
		pos = paramsDict["path_length"]
		posStopDict[node] = pos

	actions.addAction(accNodeEntranceAction, AccNode.ENTRANCE)	
	actions.addAction(accNodeExitAction, AccNode.EXIT)
	accLattice.trackActions(actions,paramsDict)
	
	posStartStopDict = {}
	for key in posStartDict:
		pos_start = posStartDict[key]
		pos_end = posStopDict[key]
		posStartStopDict[key] = (pos_start,pos_end)

	return posStartStopDict

def getAllNodesInLattice(accLattice):
	"""
	Returns the array with all nodes on all levels (even sub-child).
	"""
	
	nodes = []
	paramsDict = {}
	actions = AccActionsContainer()

	def accNodeEntranceAction(paramsDict):
		"""
		Non-bound function. Add the node to the nodes array
		at the entrance inside the node.
		This is a closure (well, maybe not
		exactly). It uses external objects.
		"""
		node = paramsDict["node"]
		nodes.append(node)
		
	actions.addAction(accNodeEntranceAction, AccNode.ENTRANCE)
	accLattice.trackActions(actions, paramsDict)
	return nodes

def getAllMagnetsInLattice(accLattice):
	"""
	Returns the array with all magnets including magnets on all levels (e.g correctors)
	"""
	
	magnets = []
	paramsDict = {}
	actions = AccActionsContainer()

	def accNodeExitAction(paramsDict):
		"""
		Non-bound function. Finds the node in the lattice
		which is a magnet.
		This is a closure (well, maybe not
		exactly). It uses external objects.
		"""
		node = paramsDict["node"]
		if(isinstance(node,LinacMagnetNode)):
			magnets.append(node)
		
	actions.addAction(accNodeExitAction, AccNode.EXIT)
	accLattice.trackActions(actions, paramsDict)
	return 	magnets
