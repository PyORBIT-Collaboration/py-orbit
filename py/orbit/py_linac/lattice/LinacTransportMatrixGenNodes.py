#!/usr/bin/env python

#--------------------------------------------------------
# Linac Accelerator Nodes for Transport Matrices generation
# These nodes are using the Initial Coordinates particles Attributes.
# Each node (if it is not the first one) calculates the transport matrix
# between the previous node and itself. 
# The matrix is a 7x7 matrix that transforms the initial particles 
# coordinates to the final ones that are in the bunch.
#--------------------------------------------------------

import math
import sys
import os

# import the finalization function 
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from orbit.py_linac.lattice import MarkerLinacNode

import orbit_utils
from orbit_utils import bunch_utils_functions
from bunch_utils_functions import copyCoordsToInitCoordsAttr
from bunch_utils_functions import transportMtrxFromInitCoords

from orbit_utils import Matrix

class LinacTrMatrixGenNode(MarkerLinacNode):
	"""
	Linac Accelerator Nodes for Transport Matrices generation.
	These nodes are using thethe Initial Coordinates particles Attrubutes.
	Each node (if it is not the first one) calculates the transport matrix
	between the previous node and itself. 
	The matrix is a 7x7 matrix that transforms the initial particles 
	coordinates to the final ones that are in the bu	
	"""
	def __init__(self, trMatricesController, name = "TrMatrixGen"):
		if(name == "TrMatrixGen"):
		 name += name+":"+str(trMatricesController.getCount())
		MarkerLinacNode.__init__(self,name)
		self.trMatricesController = trMatricesController
		self.trMtrxNode_ind = trMatricesController.getCount()
		self.use_twiss_weight_x = 0
		self.use_twiss_weight_y = 0
		self.use_twiss_weight_z = 0
		self.relativistic_beta = 0.
		self.relativistic_gamma = 0.
		#--------------------------------------
		self.trMtrx = Matrix(7,7)
		#--------------------------------------
		self.trMatricesController.addNode(self)		

	def setInternalIndex(self,ind):
		"""
		Sets the index of the TrMatrxGenNode in the controller
		"""
		self.trMtrxNode_ind = ind

	def getTrMatricesController(self):
		"""
		Returns the LinacTrMatricesContrioller that keeps the references to the TrMatrxGenNodes.
		"""
		return self.trMatricesController

	def getTwissWeightUse(self):
		"""
		Returns (use_x,use,use_z) tuple where use_{} == 1 means the Twiss weights will be used.
		"""
		res_arr = [True,True,True]
		if(self.use_twiss_weight_x == 0): res_arr[0] = False
		if(self.use_twiss_weight_y == 0): res_arr[1] = False
		if(self.use_twiss_weight_z == 0): res_arr[2] = False
		return tuple(res_arr)
		
	def setTwissWeightUse(self,use_twiss_weight_x,use_twiss_weight_y,use_twiss_weight_z):
		"""
		Sets (use_x,use,use_z) tuple where use_{} == 1 means the Twiss weights will be used.
		"""
		self.use_twiss_weight_x = 0
		self.use_twiss_weight_y = 0
		self.use_twiss_weight_z = 0		
		if(use_twiss_weight_x == True): self.use_twiss_weight_x = 1
		if(use_twiss_weight_y == True): self.use_twiss_weight_y = 1
		if(use_twiss_weight_z == True): self.use_twiss_weight_z = 1

	def track(self, paramsDict):
		bunch = paramsDict["bunch"]
		self.relativistic_beta = bunch.getSyncParticle().beta()
		self.relativistic_gamma = bunch.getSyncParticle().gamma()
		if(self.trMtrxNode_ind == 0):
			self.trMtrx.unit()
			copyCoordsToInitCoordsAttr(bunch)
		else:
			transportMtrxFromInitCoords(bunch,self.trMtrx,self.use_twiss_weight_x,self.use_twiss_weight_y,self.use_twiss_weight_z)

	def trackDesign(self, paramsDict):
		"""
		This method does nothing for the aperture case.
		"""
		pass
	
	def getBeta(self):
		"""
		Returns relativistic beta at this node.
		"""
		return self.relativistic_beta
		
	def getGamma(self):
		"""
		Returns relativistic gamma at this node.
		"""
		return self.relativistic_gamma
		
	def getTransportMatrix(self):
		"""
		Return transport matrix (7x7).
		"""
		return self.trMtrx
		
	def getDetXYZ(self, trMtrx = None):
		"""
		Returns the determinants of the transformations in (x,y,z) directions.
		"""
		if(trMtrx == None): trMtrx = self.trMtrx
		det_x = trMtrx.get(0,0)*trMtrx.get(1,1) - trMtrx.get(1,0)*trMtrx.get(0,1)
		det_y = trMtrx.get(0+2,0+2)*trMtrx.get(1+2,1+2) - trMtrx.get(1+2,0+2)*trMtrx.get(0+2,1+2)
		det_z = trMtrx.get(0+4,0+4)*trMtrx.get(1+4,1+4) - trMtrx.get(1+4,0+4)*trMtrx.get(0+4,1+4)
		return (det_x,det_y,det_z)
		
	def getNormDetXYZ(self):
		"""
		Returns the normalized determinants of the transformations in (x,y,z) directions.
		"""
		(node0,node1) = self.getTwoNodes()
		beta_in = node0.getBeta()
		beta_out = node1.getBeta()
		gamma_in = node0.getGamma()
		gamma_out = node1.getGamma()
		gb_in = beta_in*gamma_in
		gb_out = beta_out*gamma_out
		(det_x,det_y,det_z) = self.getDetXYZ(self.trMtrx)
		return ((gb_out/gb_in)*det_x,(gb_out/gb_in)*det_y,(beta_in/beta_out)*det_z)
		
	def getTwoNodes(self):
		"""
		Returns two LinacTrMatrixGenNode nodes. The transport matrix is between these nodes.
		"""
		node0 = self
		if(self.trMtrxNode_ind > 0):
			node0 = self.trMatricesController.getNode(0)
		node1 = self.trMatricesController.getNode(self.trMtrxNode_ind)
		return (node0,node1)

	def printMatrix(self):
		"""
		Print the matrix.
		"""
		name0 = "None"
		if(self.trMtrxNode_ind > 0):
			name0 = self.trMatricesController.getNode(self.trMtrxNode_ind-1).getName()
		name1 = self.trMatricesController.getNode(self.trMtrxNode_ind).getName()
		print "----Transport matrix--- from name0=",name0," to name1=",name1
		m = self.trMtrx
		for i in xrange(m.size()[0]):
			for j in xrange(m.size()[1]):
				print ("m(" + str(i) + "," + str(j)+")="+"%12.5g"%m.get(i,j) + " "),
			print " "		
	

class LinacTrMatricesContrioller:
	"""
	LinacTrMatricesContrioller keeps the references to the LinacTrMatrixGenNode
	instances.
	"""
	def __init__(self):
		self.trMatrxNodes = []
		
	def getCount(self):
		return len(self.trMatrxNodes)
		
	def getNode(self,ind):
		return self.trMatrxNodes[ind]
		
	def addNode(self,trMatrxNode):
		self.trMatrxNodes.append(trMatrxNode)
		
	def getNodes(self):
		return self.trMatrxNodes
		
	def init(self):
		#--- place nodes in the right order
		nodes = self.trMatrxNodes
		self.trMatrxNodes = sorted(nodes, key = lambda x: x.getPosition(), reverse = False)
		for node_ind in range(len(self.trMatrxNodes)):
			node = self.trMatrxNodes[node_ind]
			node.setInternalIndex(node_ind)
			
	def addTrMatrxGenNodes(self, accLattice, node_or_nodes, place = MarkerLinacNode.ENTRANCE):
		"""
		Adds the LinacTrMatrixGenNode to the nodes as child nodes.
		"""
		nodes = []
		if(type(node_or_nodes) in [tuple,list]):
			for node in node_or_nodes:
				nodes.append(node)
		else:
			nodes.append(node_or_nodes)
		#-----------------------------
		for node in nodes:
			trMatrxGenNode = LinacTrMatrixGenNode(self,node.getName()+":trMatrx")	
			node.addChildNode(trMatrxGenNode,place)
		#----- set up the position of the TrMatrix nodes
		actions = AccActionsContainer()
		
		def accNodeExitAction(paramsDict):
			"""
			Nonbound function. Sets the position of the TrMatrix nodes.
			"""
			node = paramsDict["node"]
			if(isinstance(node,LinacTrMatrixGenNode)):
				pos = paramsDict["path_length"]
				node.setPosition(pos)		
				
		actions.addAction(accNodeExitAction, AccNode.EXIT)
		accLattice.trackActions(actions)
		self.init()
		return self.trMatrxNodes
			
	def addTrMatrxGenNodesAtEntrance(self, accLattice, node_or_nodes):
		"""
		Adds the LinacTrMatrixGenNode to the nodes as child nodes at the entrance.
		"""
		return self.addTrMatrxGenNodes(accLattice, node_or_nodes, MarkerLinacNode.ENTRANCE)
		
	def addTrMatrxGenNodesAtExit(self, accLattice, node_or_nodes):
		"""
		Adds the LinacTrMatrixGenNode to the nodes as child nodes at the exit.
		"""		
		return self.addTrMatrxGenNodes(accLattice, node_or_nodes, MarkerLinacNode.EXIT)	
			
