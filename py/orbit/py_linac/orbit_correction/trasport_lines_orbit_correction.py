#!/usr/bin/env python

import math
import sys

from orbit.py_linac.lattice import Quad
from orbit.py_linac.lattice import DCorrectorH, DCorrectorV
from orbit.py_linac.lattice import MarkerLinacNode
from orbit.py_linac.lattice import BaseLinacNode

class TransverseBPM(BaseLinacNode):
	"""
	BPM node for transverse position report
	"""
	def __init__(self,bpm_marker_node):
		name = bpm_marker_node.getName()+"_diag"
		BaseLinacNode.__init__(self,name)
		self.bpm_marker_node = bpm_marker_node
		self.x = 0.
		self.y = 0.
		self.xp = 0.
		self.yp = 0.		

	def track(self, paramsDict):
		"""
		It is tracking the bunch through this node.
		"""
		bunch = paramsDict["bunch"]
		nParts = bunch.getSize()
		if(nParts != 1):
			print "debug TransverseBPM class nParts=",nParts
			print "debug It should be 1"
			print "debug Stop"
			sys.exit(1)
		self.x = bunch.x(0)
		self.y = bunch.y(0)
		self.xp = bunch.xp(0)
		self.yp = bunch.yp(0)
		
	def getGeCoordinates(self):
		"""
		returns coordinates of the particle
		"""
		return (self.x,self.xp,self.y,self.yp)

class TrajectoryCorrection:
	"""
	Class to perform correction of the trajectory in the linac type lattice
	or in the part of it. To correct trajectories the Linac DCorrectorH and 
	DCorrectorV objects will be used. The trajectory is defined by BPMs nodes
	that specifically defined for this trajectory correction. The correction
	algorithm will build response matrix for correctors (assuming no tilt) and
	apply the correctors field to provide the best correction to the target
	values for each BPM. By default the goal for all BPMs is 0., but user can 
	redefine these values and set of BPMs. User also can redefine what 
	dipole correctors will be used.
	"""
	def __init__(self, lattice, start_node = None, stop_node = None):
		self.lattice = lattice
		self.start_node = start_node
		self.stop_node = stop_node
		self.bpm_node_arr = []
		
	def _getStartStopIndexes(self):
		"""
		Returns the start and stop indexes of the start and stop nodes
		"""
		start_ind = -1
		if(self.start_node != None):
			start_ind = self.lattice.getNodeIndex(self.start_node)
		stop_ind = -1
		if(self.stop_node != None):
			stop_ind = self.lattice.getNodeIndex(self.stop_node)
		return (start_ind,stop_ind)
		
	def _returnFilteredNodes(self,nodes):
		"""
		Returns list of nodes between start and stop nodes.
		"""
		(start_ind,stop_ind) = self._getStartStopIndexes()
		nodes_tmp = []
		for node in nodes:
			ind = self.lattice.getNodeIndex(node)
			res = True
			if(start_ind >= 0):
				if(ind < start_ind):
					res = False
			if(stop_ind >= 0):
				if(ind > stop_ind):
					res = False
			if(res):
				nodes_tmp.append(node)
		return nodes_tmp
		
	def _updateBPM_Nodes(self):
		"""
		Updates BPM nodes between start_node and stop_node
		"""
		self.bpm_node_arr = []
		markers = self.lattice.getNodesOfClass(MarkerLinacNode)
		bpm_nodes = []
		for node in markers:
			if(node.getName().find("BPM") >= 0):
				bpm_nodes.append(node)
		self.bpm_node_arr = self._returnFilteredNodes(bpm_nodes)

			
		
				
		

		