#!/usr/bin/env python

import math
import sys
import os

# import the utilities
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

from bunch import Bunch
from bunch import BunchTwissAnalysis

from orbit.py_linac.lattice import Quad
from orbit.py_linac.lattice import DCorrectorH, DCorrectorV
from orbit.py_linac.lattice import MarkerLinacNode
from orbit.py_linac.lattice import BaseLinacNode

from orbit_utils import Matrix, PhaseVector

def printM(m):
	print "----matrix--- size=",m.size()
	for i in xrange(m.size()[0]):
		for j in xrange(m.size()[1]):
			print ("m(" + str(i) + "," + str(j)+")= %10.3g "%m.get(i,j) + " "),
		print ""


class TransverseBPM(BaseLinacNode):
	"""
	BPM node for transverse position report
	"""
	def __init__(self,trajCorrection,bpm_marker_node, sfx = ""):
		name = bpm_marker_node.getName()+"_diag" + sfx
		BaseLinacNode.__init__(self,name)
		#---- the trajectory correction instance that created this node
		self.trajCorrection = trajCorrection
		self.bpm_marker_node = bpm_marker_node
		self.x = 0.
		self.y = 0.
		self.xp = 0.
		self.yp = 0.

	def trackDesign(self, paramsDict):
		"""
		Nothing should happen here
		"""
		pass

	def track(self, paramsDict):
		"""
		It is tracking the bunch through this node.
		"""
		bunch = paramsDict["bunch"]
		nParts = bunch.getSize()
		if(nParts < 1):
			msg = "TransverseBPM class:"
			msg = msg + os.linesep
			msg = msg + "Node name="+self.getName()
			msg = msg + os.linesep
			msg = msg + "Number of macro-particles in bunch nParts=" + str(nParts)
			msg = msg + os.linesep
			msg = msg + "It should be more than 0"
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)
		elif(nParts > 1):
			bunch = self.trajCorrection.makeOneParticleBunch(bunch)
		self.x = bunch.x(0)
		self.y = bunch.y(0)
		self.xp = bunch.xp(0)
		self.yp = bunch.yp(0)

	def getCoordinates(self):
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
		#---- change in the corrector field
		self.deltaB = 0.001
		#----------------------------------
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
		#-------------------------
		self.twiss_analysis = BunchTwissAnalysis()

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
		node_pos_dict = self.lattice.getNodePositionsDict()
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
			(start_pos,end_pos) = node_pos_dict[bpm]
			transvBPM = TransverseBPM(self,bpm)
			transvBPM.setPosition(start_pos)
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
			nodes_tmp = self.lattice.getNodes()
			if(len(nodes_tmp) == 0): return
			(start_ind,stop_ind) = self._getStartStopIndexes()
			if(start_ind < 0): start_ind = 0
			if(stop_ind < 0): stop_ind = len(nodes_tmp)
			nodes_tmp = nodes_tmp[start_ind:stop_ind]
			for node in nodes_tmp:
				if(isinstance(node,class_type)):
					dc_node_arr.append(node)
				else:
					child_nodes = node.getAllChildren()
					for child_node in child_nodes:
						if(isinstance(child_node,class_type)):
							dc_node_arr.append(child_node)						
		else:
			for node in nodes:
				if(isinstance(node,class_type)):
					dc_node_arr.append(node)
				else:
					child_nodes = node.getAllChildren()
					for child_node in child_nodes:
						if(isinstance(child_node,class_type)):
							dc_node_arr.append(child_node)
		return dc_node_arr

	def _updateQuad_Nodes(self, nodes = None):
		"""
		Updates Quad nodes with TransverseBPM instances
		"""
		node_pos_dict = self.lattice.getNodePositionsDict()
		self.cleanQuad_Nodes()
		self.quad_node_arr = []
		quad_arr = []
		if(nodes == None):
			quad_arr = self.lattice.getNodesOfClass(Quad)
		else:
			quad_arr = nodes
		quad_arr = self._returnFilteredNodes(quad_arr)
		for quad in quad_arr:
			(start_pos,end_pos) = node_pos_dict[quad]
			transvBPM = TransverseBPM(self,quad,"_entrance")
			transvBPM.setPosition(start_pos)
			quad.addChildNode(transvBPM,AccNode.ENTRANCE)
			self.quad_transvBPM_arr.append(transvBPM)
			transvBPM = TransverseBPM(self,quad,"_exit")
			transvBPM.setPosition(end_pos)
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

	def setDCVs(self, dcvs):
		"""
		Sets the DCorrectorVs
		"""
		self._updateDC_Nodes(dcvs,DCorrectorV)
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

	def correctTrajectory(self,bunch_initial):
		"""
		This method will calculate and applyes fields to dipole correctors
		to achive minimal beam deviation from the center at all BPMs.
		LSQM method is used.
		"""
		bunch_in = self.makeOneParticleBunch(bunch_initial)
		#--------------------------
		bunch_init = Bunch()
		bunch_in.copyEmptyBunchTo(bunch_init)
		#---- add one particle to the bunch_init
		x0  = 0.
		xp0 = 0.
		y0  = 0.
		yp0 = 0.
		z0  = 0.
		dE  = 0.
		bunch_init.addParticle(x0,xp0,y0,yp0,z0,dE)
		#----------------------------------
		n_bpms = len(self.bpm_node_arr)
		n_corrs_hor = len(self.dch_node_arr)
		n_corrs_ver = len(self.dcv_node_arr)
		horResponceMtrx = Matrix(n_bpms,n_corrs_hor)
		verResponceMtrx = Matrix(n_bpms,n_corrs_ver)
		self._calculatBPM_Matrix(bunch_init,horResponceMtrx,self.dch_node_arr,axis = 0)
		self._calculatBPM_Matrix(bunch_init,verResponceMtrx,self.dcv_node_arr,axis = 1)
		#printM(horResponceMtrx)
		horResponceMtrxTr = horResponceMtrx.copy()
		horResponceMtrxTr.transpose()
		verResponceMtrxTr = verResponceMtrx.copy()
		verResponceMtrxTr.transpose()
		horAmtrx = (horResponceMtrxTr.mult(horResponceMtrx)).invert()
		verAmtrx = (verResponceMtrxTr.mult(verResponceMtrx)).invert()
		#printM(horAmtrx)
		#print "det(horAmtrx) = ",horAmtrx.det()
		horATmtrx = horAmtrx.mult(horResponceMtrxTr)
		verATmtrx = verAmtrx.mult(verResponceMtrxTr)
		#---- matrices for LSQM are ready - now use initial bunch to get
		#---- the trajectory that should be flattened
		bunch_init = Bunch()
		bunch_in.copyBunchTo(bunch_init)
		(bpm_value_hor_arr,bpm_value_ver_arr) = self.calculateTrajectory(bunch_init)
		horBPM_V = PhaseVector(len(bpm_value_hor_arr))
		for ind in range(len(bpm_value_hor_arr)):
			horBPM_V.set(ind,bpm_value_hor_arr[ind])
		verBPM_V = PhaseVector(len(bpm_value_ver_arr))
		for ind in range(len(bpm_value_ver_arr)):
			verBPM_V.set(ind,bpm_value_ver_arr[ind])
		dch_val_V = horATmtrx.mult(horBPM_V)
		dcv_val_V = verATmtrx.mult(verBPM_V)
		#print "debug =========== final DCH fields ================"
		for ind in range(dch_val_V.size()):
			dc = self.dch_node_arr[ind]
			delta_field = dch_val_V.get(ind)
			#print "debug corr =",dc.getName()," init field [T] =",dc.getParam("B")," delta [T] =",delta_field
			dc.setParam("B",dc.getParam("B") - delta_field)
		#print "debug =========== final DCV fields ================"
		for ind in range(dcv_val_V.size()):
			dc = self.dcv_node_arr[ind]
			delta_field = dcv_val_V.get(ind)
			#print "debug corr =",dc.getName()," init field [T] =",dc.getParam("B")," delta [T] =",delta_field
			dc.setParam("B",dc.getParam("B") - delta_field)		

	def _calculatBPM_Matrix(self,bunch_init,responceMtrx,dc_node_arr,axis = None):
		corr_field_arr = []
		for dc_node in dc_node_arr:
			corr_field_arr.append(dc_node.getParam("B"))
		#---------------------------------------------
		bpm_value_init_arr = self.calculateTrajectory(bunch_init)[axis]
		for dc_node_ind in range(len(dc_node_arr)):
			dc_node = dc_node_arr[dc_node_ind]
			field = corr_field_arr[dc_node_ind] + self.deltaB
			dc_node.setParam("B",field)
			bpm_value_arr = self.calculateTrajectory(bunch_init)[axis]
			for bpm_ind in range(len(self.transvBPM_arr)):
				deltaVal = bpm_value_arr[bpm_ind] - bpm_value_init_arr[bpm_ind]
				derivative = deltaVal/self.deltaB
				responceMtrx.set(bpm_ind,dc_node_ind,derivative)
			dc_node.setParam("B",corr_field_arr[dc_node_ind])

	def calculateTrajectory(self,bunch_init, print_info = False):
		"""
		It tracks bunch with one particle through the lattice and
		fills out the BPM data in self.transvBPM_arr
		Returns arrays of x,y corrdinates of the bunch in
		(pm_value_hor_arr,bpm_value_ver_arr) tuple.
		"""
		bunch = Bunch()
		bunch_init.copyBunchTo(bunch)
		(start_ind,stop_ind) = self._getStartStopIndexes()
		self.lattice.trackDesignBunch(bunch,None,None,start_ind,stop_ind)
		self.lattice.trackBunch(bunch,None,None,start_ind,stop_ind)
		bpm_value_hor_arr = []
		bpm_value_ver_arr = []
		for transvBPM in self.transvBPM_arr:
			val_hor = transvBPM .getCoordinates()[0]
			val_ver = transvBPM .getCoordinates()[2]
			bpm_value_hor_arr.append(val_hor)
			bpm_value_ver_arr.append(val_ver)
		#---- Start trajectory printing
		if(print_info):
			print " index name position[m]  x[mm] y[mm] "
			for ind in range(len(bpm_value_hor_arr)):
				hor_val = bpm_value_hor_arr[ind]*1000.
				ver_val = bpm_value_ver_arr[ind]*1000.
				bpm_node = self.bpm_node_arr[ind]
				transvBPM = self.transvBPM_arr[ind]
				pos = transvBPM.getPosition()
				print " %2d %20s  %8.3f   %+6.2f %+6.2f "%(ind,bpm_node.getName(),pos,hor_val,ver_val)
		#---- Stop trajectory printing
		return (bpm_value_hor_arr,bpm_value_ver_arr)

	def makeOneParticleBunch(self,bunch_init):
		"""
		Returns the bunch with one particle which coordinates
		are equal to the average coordinates in the initial
		bunch.
		"""
		bunch_out = Bunch()
		bunch_init.copyEmptyBunchTo(bunch_out)
		x_avg  = 0.
		xp_avg = 0.
		y_avg  = 0.
		yp_avg = 0.
		z_avg  = 0.
		dE_avg = 0.
		nParts = bunch_init.getSizeGlobal()
		if(nParts > 0):
			self.twiss_analysis.analyzeBunch(bunch_init)
			x_avg  = self.twiss_analysis.getAverage(0)
			xp_avg = self.twiss_analysis.getAverage(1)
			y_avg  = self.twiss_analysis.getAverage(2)
			yp_avg = self.twiss_analysis.getAverage(3)
			z_avg  = self.twiss_analysis.getAverage(4)
			dE_avg = self.twiss_analysis.getAverage(5)
		#----------------------------------------
		#---- add one particle to the bunch_out
		bunch_out.addParticle(x_avg,xp_avg,y_avg,yp_avg,z_avg,dE_avg)
		return bunch_out














