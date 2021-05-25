#!/usr/bin/env python

import math
import sys

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

from orbit.py_linac.lattice import Quad
from orbit.py_linac.lattice import DCorrectorH, DCorrectorV
from orbit.py_linac.lattice import MarkerLinacNode
from orbit.py_linac.lattice import BaseLinacNode

class TransverseBPM(BaseLinacNode):
	"""
	BPM node for transverse position report
	"""
	def __init__(self,bpm_marker_node, sfx = ""):
		name = bpm_marker_node.getName()+"_diag" + sfx
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
		self.transvBPM_arr = []
		self.dch_node_arr = []
		self.dcv_node_arr = []
		#-------------------------
		self.quad_node_arr = []
		self.quad_transvBPM_arr = []
		#-------------------------
		self._updateBPM_Nodes()
		self._updateDC_Nodes(None,DCorrectorH)
		self._updateDC_Nodes(None,DCorrectorV)
		self._updateQuad_Nodes()
		
	def setStartStopNodes(self,start_node = None, stop_node = None):
		self.start_node = start_node
		self.stop_node = stop_node
		self._updateBPM_Nodes()
		self._updateDC_Nodes(None,DCorrectorH)
		self._updateDC_Nodes(None,DCorrectorV)
		self._updateQuad_Nodes()
		
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
		
	def _updateBPM_Nodes(self, bpms = None):
		"""
		Updates BPM nodes between start_node and stop_node.
		Returns array of bpm nodes.
		"""
		self.cleanBPM_Nodes()
		self.bpm_node_arr = []
		bpm_nodes = []
		if(bpms == None):
			markers = self.lattice.getNodesOfClass(MarkerLinacNode)
			for node in markers:
				if(node.getName().find("BPM") >= 0):
					bpm_nodes.append(node)
		else:
			bpm_nodes = bpms
		self.bpm_node_arr = self._returnFilteredNodes(bpm_nodes)
		for bpm in self.bpm_node_arr:
			transvBPM = TransverseBPM(bpm)
			bpm.addChildNode(transvBPM,AccNode.ENTRANCE)
			self.transvBPM_arr.append(transvBPM)
		return self.bpm_node_arr

	def _updateDC_Nodes(self, nodes = None, class_type = None):
		"""
		Updates DCorrector nodes that we use for orbit correction
		"""
		if(class_type == None): return
		dc_node_arr = None
		if(class_type == DCorrectorH):
			dc_node_arr = self.dch_node_arr
		if(class_type == DCorrectorV):
			dc_node_arr = self.dcv_node_arr
		if(dc_node_arr == None): return
		del dc_node_arr[:]
		if(nodes == None):
			dc_node_arr += self.lattice.getNodesOfClass(class_type)
		else:
			for node in nodes:
				if(isinstance(node,class_type)):
					dc_node_arr.append(node)
		node_arr = self._returnFilteredNodes(dc_node_arr)
		del dc_node_arr[:]
		dc_node_arr += node_arr
		return dc_node_arr
		
	def _updateQuad_Nodes(self, nodes = None):
		"""
		Updates Quad nodes with TransverseBPM instances
		"""
		self.cleanQuad_Nodes()
		self.quad_node_arr = []
		quad_arr = []
		if(nodes == None):
			quad_arr = self.lattice.getNodesOfClass(Quad)
		else:
			quad_arr = nodes
		quad_arr = self._returnFilteredNodes(quad_arr)
		for quad in quad_arr:
			transvBPM = TransverseBPM(quad,"_entrance")
			quad.addChildNode(transvBPM,AccNode.ENTRANCE)
			self.quad_transvBPM_arr.append(transvBPM)
			transvBPM = TransverseBPM(quad,"_exit")
			quad.addChildNode(transvBPM,AccNode.EXIT)
			self.quad_transvBPM_arr.append(transvBPM)
		self.quad_node_arr = quad_arr
		return self.quad_node_arr		
			
	def setBPMs(self,bpms):
		"""
		Sets custom bpm array for analysis.
		In reality they can be any nodes.
		Returns array of bpm nodes.
		"""
		return self._updateBPM_Nodes(bpms)
		
	def getBPMs(self):
		"""
		Returns array of bpm nodes.
		"""
		return self.bpm_node_arr
		
	def setDCHs(self, dchs):
		"""
		Sets the DCorrectorHs 
		"""
		self._updateDC_Nodes(dchs,DCorrectorH)
		return self.dch_node_arr
		
	def getDCHs(self):
		"""
		Returns DCorrectorHs array
		"""
		return self.dch_node_arr
		
	def setDCHVs(self, dchs):
		"""
		Sets the DCorrectorVs 
		"""
		self._updateDC_Nodes(dchs,DCorrectorV)
		return self.dcv_node_arr
		
	def getDCVs(self):
		"""
		Returns DCorrectorVs array
		"""
		return self.dcv_node_arr	
		
	def setQuads(self, quads):
		"""
		Sets Quads
		"""
		self._updateQuad_Nodes(quads)
		return self.quad_node_arr
		
	def getQuads(self):
		"""
		Returns Quads array
		"""
		return self.quad_node_arr			
		
	def getTransverseBPMs(self):
		return self.transvBPM_arr
		
	def getQuadTransverseBPMs(self):
		return self.quad_transvBPM_arr
		
	def getTransverseBPMforBPM(self,bpm):
		"""
		Retuns the TransverseBPM instance for particular BPM.
		"""
		for child in bpm.getChildNodes(AccNode.ENTRANCE):
			if(isinstance(child,TransverseBPM)):
				return child
		return None
	
	def cleanBPM_Nodes(self):
		"""
		Removes TransverseBPM child nodes from BPM nodes
		"""
		for bpm in self.bpm_node_arr:
			child_arr = bpm.getChildNodes(AccNode.ENTRANCE)
			transvBPM = None
			for child in child_arr:
				if(isinstance(child,TransverseBPM)):
					transvBPM = child
					break
			if(transvBPM != None):
				child_arr.remove(transvBPM)
		self.transvBPM_arr = []
		
	def cleanQuad_Nodes(self):
		"""
		Removes TransverseBPM child nodes from Quad nodes
		"""
		for node in self.quad_node_arr:
			for place in [AccNode.ENTRANCE,AccNode.EXIT]:
				child_arr = node.getChildNodes(place)
				transvBPM = None
				for child in child_arr:
					if(isinstance(child,TransverseBPM)):
						transvBPM = child
						break
				if(transvBPM != None):
					child_arr.remove(transvBPM)
		self.quad_transvBPM_arr = []
		
	
		

			
		
				
		

		