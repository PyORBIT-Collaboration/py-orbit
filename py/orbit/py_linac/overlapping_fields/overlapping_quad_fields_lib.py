#!/usr/bin/env python

#--------------------------------------------------------
# This is a collection of classes for describing quad 
# overlapping fields.
# Here we have only the model for overlapping quadrupole 
# magnetic fields.
#--------------------------------------------------------
# The classes for the quads with the overlapping fields uses
# the field distribution from the following paper:
#    M.Berz, B. Erdelyn, K.Makino
#    Fringe Field Effects in Small Rings of Large Acceptance
#    Phys. Rev STAB, V3, 124001(2000)
#
#--------------------------------------------------------

import math
import sys
import os

from orbit.py_linac.lattice import BaseLinacNode, Drift, Quad
from orbit.py_linac.lattice import AxisFieldRF_Gap
from overlapping_rf_and_quad_fields_lib import AxisField_and_Quad_RF_Gap

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

from orbit_utils import Function


class AbstractQuadFieldSourceFunction:
	"""
	It is an abstract class describing the quadrupole magnetic 
	field as a function of the longitudinal coordinate.
	"""
	def __init__(self):
		pass
	
	def getLimitsZ(self):
		pass
	
	def getFuncValue(self,z):
		""" 
		Returns the quad's normalized field distribution 
		at the distance z from the center.
		the distribution should be normalized to 1 over the
		integration along the longitudinal coordinate.
		"""
		pass
	
class SimpleQuadFieldFunc(AbstractQuadFieldSourceFunction):
	"""
	It is an implementation of the QuadFieldSourceFunction class for
	a simple quad with constant field between (-L/2;+L/2).
	"""
	def __init__(self,quad):
		self.quad = quad
		
	def getLimitsZ(self):
		"""
		Returns (+L/2,-L/2) as longitudinal limits of the quad field.
		"""
		L = self.quad.getLength()
		return (-L/2,L/2)
		
	def getFuncValue(self,z):
		""" 
		Returns the quad's normalized field distribution at the distance z from the center 
		"""
		L = self.quad.getLength()
		if(abs(z) <= L/2):
			return 1./L
		return 0.
		
class EngeFunction(AbstractQuadFieldSourceFunction):
	""" 
	The Enge function with parameters from Berz's paper 
	M.Berz, B. Erdelyn, K.Makino
  Fringe Field Effects in Small Rings of Large Acceptance
  Phys. Rev STAB, V3, 124001(2000)	
	"""
	def __init__(self, length_param, acceptance_diameter_param, cutoff_level = 0.01):
		self.length = length_param
		self.acceptance_diameter = acceptance_diameter_param
		self.a_arr = [0.296471,4.533219,-2.270982,1.068627,-0.036391,0.022261]	
		self.normalization= 1.0
		self.n_func_points = 500
		#-----find cut-off z value
		self.cutoff_z = self.acceptance_diameter		
		step = self.acceptance_diameter
		self.cutoff_level = cutoff_level
		self.cutoff_z = self._findCutOff(step, cutoff_level)
		#------------------------------------------------------
		self.func = Function()
		self._normalize()
		
	def setEngeCoefficients(self,a_arr):
		"""
		Sets new values for Enge function's coeffients.
		"""
		self.a_arr = a_arr
		step = self.length/2
		self.cutoff_z = self._findCutOff(step, self.cutoff_level)	
		self._normalize()
		
	def setCutOffLevel(self,cutoff_level):
		""" Sets the cutoff level for quad's  field """
		step = self.length/2
		self.cutoff_level = cutoff_level
		self.cutoff_z = self._findCutOff(step, cutoff_level)
		self._normalize()
		
	def setCutOffZ(self,cutoff_z):
		""" Sets the cutoff distance from the center of quad's field"""
		self.cutoff_z = cutoff_z
		self._normalize()
		
	def setLength(self,length):
		""" Sets the length of quad's field"""
		self.length = length
		step = self.length/2.0
		self.cutoff_z = self._findCutOff(step, self.cutoff_level)
		self._normalize()
		
	def setAcceptanceDiameter(self,acceptance_diameter):
		""" Sets the acceptance diameter of the quad """
		self.acceptance_diameter = acceptance_diameter
		step = self.length/2.0
		self.cutoff_z = self._findCutOff(step, self.cutoff_level)
		self._normalize()
		
	def setNumberOfPoints(self,n_func_points):
		""" Sets the number of points in the field function """
		self.n_func_points = n_func_points
		step = self.length/2.0
		self.cutoff_z = self._findCutOff(step, self.cutoff_level)
		self._normalize()
		
	def getCuttOffZ(self):
		""" Returns the cutoff distance from the center of quad's field"""
		return self.cutoff_z 
		
	def getNumberOfPoints(self):
		""" Returns the number of points in the field function """
		return self.n_func_points
		
	def _getTrueEngeFunc(self, x):
		""" Returns the quad's field at the distance x from the center """
		# x is the distance from the center of the magnet with the iron length l """
		x = (math.fabs(x) - self.length/2.0)/self.acceptance_diameter
		sum_exp = self.a_arr[0]
		x0 = x
		for i in range(1,len(self.a_arr)):
			sum_exp += self.a_arr[i]*x0
			x0 *= x
		if(abs(sum_exp) > 30.): sum_exp = 30.0*sum_exp/abs(sum_exp)
		return self.normalization/(1.0+math.exp(sum_exp))

	def _findCutOff(self,step, cutoff_level):
		""" Finds the distance from the center where the field is less than cutoff level """
		self.normalization = 1.0
		init_val = self._getTrueEngeFunc(0.)
		z = step
		val = self._getTrueEngeFunc(z)/init_val
		if(val <= cutoff_level):
			return z
		while(val > cutoff_level):
			z += step
			val = self._getTrueEngeFunc(z)/init_val
		z0 = z - step
		z1 = z
		step_z = step/self.n_func_points
		val0 =  self._getTrueEngeFunc(z0)/init_val
		val1 =  self._getTrueEngeFunc(z1)/init_val		
		while(abs(z0-z1) > step_z):
			z_new = (z0+z1)/2.0
			val_new = self._getTrueEngeFunc(z_new)/init_val
			if(val_new <= cutoff_level):
				z1 = z_new
				val1 = val_new
			else:
				z0 = z_new
				val0 = val_new			
		self.cutoff_z = (z0+z1)/2.0
		return self.cutoff_z
				
	def _normalize(self):
		""" Normalizes the quad field function to the integral of 1 """
		self.normalization = 1.0
		step = self.cutoff_z/(self.n_func_points - 1)
		self.func.clean()
		sum_int = 0.
		for ind in range(self.n_func_points):
			z = step*ind
			val = self._getTrueEngeFunc(z)
			self.func.add(z,val)
			sum_int += val
		sum_int -= (self._getTrueEngeFunc(0.) + self._getTrueEngeFunc(step*(self.n_func_points - 1)))/2.0
		sum_int *= 2.0*step
		self.normalization = 1.0/sum_int
		self.func.setConstStep(1)
		
	def getFuncValue(self,z):
		""" Returns the quad's field at the distance z from the center """
		if(abs(z) >= self.func.getMaxX()): return 0.
		return self.normalization*self.func.getY(abs(z))
		
	def getLimitsZ(self):
		""" Returns the tuple with min and max Z value for this field """
		z_max = self.func.getMaxX()
		return (-z_max,z_max)
		
class OverlappingQuadsNode(BaseLinacNode):
	"""
	The set of quads with the overlapping fields.
	"""
	def __init__(self, name = "OverlappingQuads"):
		"""
		Constructor. Creates the OverlappingQuadsNode instance.
		"""
		BaseLinacNode.__init__(self,name)
		self.setType("OVRLPQ")
		self.setnParts(1)
		#----quads_fields_arr is an array of [quad, fieldFunc, z_center_of_field]
		self.quads_fields_arr = []
		#---- current z position
		self.z_value = 0.
		#---- z-step - the step in longitudinal direction during the tracking in [m]
		self.z_step = 0.01 
		#---- the initial length is 0.
		self.setLength(0.)
		
	def addQuad(self, quad, fieldFunc, z_center_of_field):
		"""
		Adds the quad with the field function and the position.
		The position of the quad is relative to the beginning of this OverlappingQuadsNode.
		"""
		self.quads_fields_arr.append((quad, fieldFunc, z_center_of_field))
				
	def setZ_Step(self,z_step):
		"""
		Sets the longitudinal step for the tracking along the node.
		"""
		self.z_step = z_step
		
	def getZ_Step(self):
		"""
		Returns the longitudinal step for the tracking along the node.
		"""		
		return self.z_step
		
	def getZ_Min_Max(self):
		"""
		Returns the tuple (z_min,z_max) with the limits of z coordinate from the center.
		These parameters define the length of the node. The center of the node
		is at 0.
		"""
		L2 = self.getLength()/2
		return (-L2,L2)	
		
	def getQuads(self):
		"""
		Returns the list of quads in this node.
		"""
		quads = []
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			quads.append(quad)
		return quads
	
	def getCentersOfField(self):
		"""
		Returns the array of centers of the quads in this node.
		"""
		centers_arr = []
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			centers_arr.append(z_center_of_field)
		return centers_arr

	def initialize(self):
		"""
		The OverlappingQuadsNode class implementation
		of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		length = self.getLength()
		lengthStep = length/nParts
		for i in xrange(nParts):
			self.setLength(lengthStep,i)

	def track(self, paramsDict):
		"""
		The  OverlappingQuadsNode class implementation of the AccNode class track(paramDict) method.
		"""
		index = self.getActivePartIndex()	
		length = self.getLength(index)
		if(index == 0): self.z_value = - self.getLength()/2
		bunch = paramsDict["bunch"]
		momentum = bunch.getSyncParticle().momentum()		
		n_steps = int(length/self.z_step)+1
		z_step = length/n_steps
		for z_ind in range(n_steps):
			z = self.z_value + z_step*(z_ind+0.5)
			G = self.getTotalField(z)
			kq = G/(3.335640952*momentum)
			if(abs(kq) == 0.):
				TPB.drift(bunch,z_step)
				continue
			#------- track through a combined quad
			TPB.quad1(bunch,z_step/4.0, kq)
			TPB.quad2(bunch,z_step/2.0)
			TPB.quad1(bunch,z_step/2.0, kq)
			TPB.quad2(bunch,z_step/2.0)
			TPB.quad1(bunch,z_step/4.0, kq)
		self.z_value += length
		
	def getTotalField(self,z_from_center):
		"""
		Returns the combined field of all overlapping quads.
		z_from_center - is a distance from the center of the node.
		z - is a distance from the beginning of the node.
		"""
		z = z_from_center + self.getLength()/2
		G = 0.
		if(z < 0. or z > self.getLength()): return G
		for [quad, fieldFunc, z_center_of_field] in self.quads_fields_arr:
			if(fieldFunc != None):
				gl = quad.getParam("dB/dr")*quad.getLength()
				G += gl*fieldFunc.getFuncValue(z - z_center_of_field)
			else:
				G += quad.getTotalField(z - z_center_of_field)
		return G
	
class OverlappingQuadsController:
	"""
	A container for the set of OverlappingQuadsNode instances.
	It performs the lattice modification and allows control
	over the OverlappingQuadsNode nodes.
	This controller is not the part of the lattice.
	"""
	def __init__(self, name = "OverlappingQuadsCntrl"):
		"""
		Constructor. Creates the OverlappingQuadsController instance.
		"""
		self.name = name
		#--- The array with OverlappingQuadsNode instances
		self.overlapping_quads = []
		#--- The array with [quad, fieldFunc, z_center_of_field] lists
		self.quads_fields_arr = []
		#--- Total length of all nodes
		self.length = 0.
		#--- The longitudinal step over the magnetic field during the tracking
		self.z_step = 0.005
		#--- The maximal distance between parts of the node
		self.z_max_dst = 0.01

	def getName(self):
		return self.name
		
	def setName(self, name):
		self.name = name
	
	def getLength(self):
		return self.length
	
	def addQuadsAndFields(self,accLattice,quads_fields_arr):
		"""
		Adds the OverlappingQuadsNodes to the Acc. Lattice as a replacement 
		for the quads in the quads_fields_arr [[quad,fieldFunc]] 
		with the combination of the fieldFuncs. We should take into account 
		the existing drifts and other nodes with zero length. We cannot have 
		non-zero nodes in the area covered by overlapping quads.
		"""
		n_quads = len(quads_fields_arr)
		if(n_quads < 1):
			print "Class OverlappingQuadsController. We have to have one or more quads here!"
			print "Stop"
			sys.exit(1)
		node_pos_dict = accLattice.getNodePositionsDict()
		#print "debug dict=",node_pos_dict
		self.quads_fields_arr = []
		quads = []
		for [quad, fieldFunc] in quads_fields_arr:
			#print "debug quad=",quad.getName()
			quads.append(quad)
			(posBefore, posAfter) = node_pos_dict[quad]
			z_center_of_field = (posBefore + posAfter)/2.0
			self.quads_fields_arr.append([quad, fieldFunc, z_center_of_field])
		z_min = self.quads_fields_arr[0][1].getLimitsZ()[0]+self.quads_fields_arr[0][2]
		z_max = self.quads_fields_arr[n_quads-1][1].getLimitsZ()[1]+self.quads_fields_arr[n_quads-1][2]
		self.length = z_max - z_min
		#print "debug length=",self.length
		#---- let's create the array of the AccNodes that are inside or touching the region 
		#---- between z_min and z_max
		#---- Thin nodes should be preserved
		thin_nodes = []
		#---- the two drifts: before the overlapping region and after [drift,length,ind]
		drift_nodes = []
		for node_ind in range(len(accLattice.getNodes())):
			node = accLattice.getNodes()[node_ind]
			(posBefore, posAfter) = node_pos_dict[node]
			if(posBefore < z_min and posAfter > z_min):
				drift_nodes.append([node,z_min - posBefore,node_ind])	
			if(posBefore < z_max and posAfter > z_max):
				drift_nodes.append([node,posAfter - z_max,node_ind])
			if((posBefore >= z_min and posBefore <= z_max) or (posAfter >= z_min and posAfter <= z_max)):
				length = node.getLength()
				if(length > 0.):
					if(not ((node in quads) or isinstance(node,Drift))):
						print "Class OverlappingQuadsController. There is a node with L > 0 (not quad)!"
						print "It cannot be included in the overlapping node!"
						print "This node should be Quad or Drift."
						print "node =",node.getName()
						print "type =",node.getType()
						print "L=",length
						sys.exit(1)
				else:
					thin_nodes.append([(posBefore+posAfter)/2.0,node])
		if(len(drift_nodes) != 2 and accLattice.getNodes()[drift_nodes[0][2]].getSequence() != accLattice.getNodes()[drift_nodes[1][2]].getSequence()):
			print "Class OverlappingQuadsController. We cannot put the fields inside the lattice!"
			print "Fields' z_min=",z_min
			print "Fields' z_max=",z_max
			print "Seq0 = ",accLattice.getNodes()[drift_nodes[0][2]].getSequence().getName()
			print "Seq1 = ",accLattice.getNodes()[drift_nodes[1][2]].getSequence().getName()
			print "edge drift nodes=",drift_nodes
			print "Stop"
			sys.exit(1)	
		accSeq = accLattice.getNodes()[drift_nodes[0][2]].getSequence()
		new_lattice_nodes_before = []
		for node_ind in range(drift_nodes[0][2]):
			node = accLattice.getNodes()[node_ind]
			new_lattice_nodes_before.append(node)
		drift_before = Drift(drift_nodes[0][0].getName())
		drift_before.setLength(drift_nodes[0][1])
		drift_before.setSequence(accSeq)
		drift_before.setName(drift_nodes[0][0].getName())
		new_lattice_nodes_before.append(drift_before)
		new_lattice_nodes_after = []
		drift_after = Drift(drift_nodes[1][0].getName())
		drift_after.setLength(drift_nodes[1][1])
		drift_after.setSequence(accSeq)
		drift_after.setName(drift_nodes[1][0].getName())
		new_lattice_nodes_after.append(drift_after)
		for node_ind in range(drift_nodes[1][2]+1,len(accLattice.getNodes())):
			node = accLattice.getNodes()[node_ind]
			new_lattice_nodes_after.append(node)
		#------- the OverlappingQuadsNodes and possible markers
		#------- self.overlapping_quads include only OverlappingQuadsNode type
		#------- allNewNodes includes all nodes (even possible zero length nodes)
		self.overlapping_quads = []
		allNewNodes = []
		z_start = z_min
		z_end = z_max
		for ind in range(len(thin_nodes)+1):
			overlappingNode = OverlappingQuadsNode()
			overlappingNode.setSequence(accSeq)
			self.overlapping_quads.append(overlappingNode)
			allNewNodes.append(overlappingNode)
			z_end = z_max
			if(len(thin_nodes) > 0 and ind < len(thin_nodes)):
				z_end = thin_nodes[ind][0]
				allNewNodes.append(thin_nodes[ind][1])
			overlapping_length =  z_end - z_start
			if(overlapping_length < 0.):
				print "Class OverlappingQuadsController. Something wrong about OverlappingQuadsNode!"
				print "negative length=",overlapping_length
				print "before the thin node=",thin_nodes[ind][1].getName()
				print "Stop"
				sys.exit(1)					
			for quad_ind in range(len(self.quads_fields_arr)):
				[quad, fieldFunc, z_center_of_field] = self.quads_fields_arr[quad_ind]
				z_pos = z_center_of_field - z_start
				overlappingNode.addQuad(quad, fieldFunc, z_pos)
			overlappingNode.setLength(overlapping_length)
			overlappingNode.initialize()
			z_start = z_end
		if(len(self.overlapping_quads) > 1):
			count = 0
			for overlappingNode in self.overlapping_quads:
				overlappingNode.setName(self.getName()+"_p"+str(count))
				count += 1
		else:
			self.overlapping_quads[0].setName(self.getName())
		allNodes = new_lattice_nodes_before + allNewNodes + new_lattice_nodes_after
		accLattice.setNodes(allNodes)
		accLattice.initialize()
		
	def setZ_Step(self,z_step):
		"""
		Sets the longitudinal step for the tracking along the overlapping nodes.
		"""
		self.z_step = z_step
		for overlappingNode in self.overlapping_quads:
			overlappingNode.setZ_Step(z_step)
		
	def getZ_Step(self):
		"""
		Returns the longitudinal step for the tracking along the overlapping nodes.
		"""		
		return self.z_step
		
	def setMaxPartsDistance(self,z_max_dst):
		"""
		Sets the number of parts in the overlapping nodes, so
		the max distance between parts is less than the z_max_dst.
		We need this for the space charge nodes settings.
		"""
		self.z_max_dst = z_max_dst
		for overlappingNode in self.overlapping_quads:
			length = overlappingNode.getLength()
			nParts = int(length/z_max_dst)+1
			overlappingNode.setnParts(nParts)
	
	def getMaxPartsDistance(self):
		"""
		Returns the number of parts in the overlapping nodes, so
		the max distance between parts is less than the z_max_dst.
		We need this for the space charge nodes settings.
		"""		
		return self.z_max_dst
		
	def getOverlappingNodes(self):
		"""
		Returns the overlapping nodes.
		"""
		return self.overlapping_quads
		
def GetGlobalQuadGradient(accLattice,z):
	"""
	The service function for the overlapping fields package.
	Returns the quad field for certain position in the lattice that
	has usual and overlapping quads.
	"""
	node_pos_dict = accLattice.getNodePositionsDict()
	nodes = accLattice.getNodes()
	G = 0.
	(node,posBefore,posAfter) = GetNodeForPosition(accLattice,z)
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
	(node,posBefore,posAfter) = GetNodeForPosition(accLattice,z)
	if(isinstance(node,AxisField_and_Quad_RF_Gap) or isinstance(node,AxisFieldRF_Gap)):
		(z_min,z_max) = node.getZ_Min_Max()
		Ez = node.getEzFiled(z - posBefore + z_min)
	return Ez

def GetNodeForPosition(accLattice,z):
	"""
	It is a local convinience function. It returns the node which
	coordsinates cover the z-position.
	"""
	node_pos_dict = accLattice.getNodePositionsDict()
	nodes = accLattice.getNodes()	
	index0 = 0
	index1 = len(nodes) - 1
	(posBefore0, posAfter0) = node_pos_dict[nodes[index0]]
	(posBefore1, posAfter1) = node_pos_dict[nodes[index1]]
	index = 0
	while(index0 != index1):			
		index = (index0 + index1)/2
		#print "debug z=",z," index0=",index0," index1=",index1," index=",index," (posBefore0, posAfter0)=",(posBefore0, posAfter0)," (posBefore1, posAfter1)=",(posBefore1, posAfter1)
		(posBefore, posAfter) = node_pos_dict[nodes[index]]
		if(z < posBefore):
			index1 = index
			(posBefore1, posAfter1) = node_pos_dict[nodes[index1]]
		elif(z > posAfter):
			index0 = index
			(posBefore0, posAfter0) = node_pos_dict[nodes[index0]]
		elif(z >= posBefore and z <= posAfter):
			break
		if(z >= posBefore0 and z <= posAfter0):
			index =index0 
			break
		if(z >= posBefore1 and z <= posAfter1):
			index =index1
			break		
	#-------------------------------------------
	node = nodes[index]
	(posBefore, posAfter) = node_pos_dict[node]
	return (node,posBefore,posAfter)

#-----------------------------------------------------------------------		
#-----Test of the Enge Function ----------------	
#-----------------------------------------------------------------------
if __name__ == "__main__":	
	#---- MEBT quads ----
	length_param = 0.066
	acceptance_diameter_param = 0.0363
	#---- DTL Permanent Quad
	length_param = 0.035
	acceptance_diameter_param = 0.025
	#--------------------
	cutoff_level = 0.01	
	func = EngeFunction(length_param,acceptance_diameter_param,cutoff_level)
	z_max = func.getCuttOffZ()
	np = func.getNumberOfPoints()
	np = 50
	step = z_max/(np-1)
	for ind in range(2*np-1):
		z = ind*step-z_max
		val = func.getFuncValue(z)
		print " %12.5e  %12.5e "%(z*1000.,val)
	print "z limits=",func.getLimitsZ()
	func.setCutOffZ(0.5*func.getLimitsZ()[1])
	print "new z limits=",func.getLimitsZ()



