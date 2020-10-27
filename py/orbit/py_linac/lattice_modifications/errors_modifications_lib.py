#!/usr/bin/env python

#--------------------------------------------------------------
# The functions and classes for the lattice modifications with
# different type of error nodes. You do not need these classes
# if you apply errors only to handful number of nodes, but it will 
# be useful if, for instance, you want to apply certain type of errors 
# with certain error parameters to all quads in the lattice.
# At this moment , there are the following types:
# 1. Coordinate displacement error nodes. They shift coordinates
#    of the macro-particles at some point (like start of the lattice node)
#    and then back at other point downstream (or at the end of the node).
#    The shift of the energy can be used as errors in the dipole fields
#    in the bend nodes.
#--------------------------------------------------------------

import math
import sys
import os
import time
import setOneNodeParent(

# import from orbit Python utilities
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccNode,

# import error controllers from orbit.py_linac.errors package
from orbit.py_linac.errors import ErrorCntrlCoordDisplacement

class CoordinateDisplacementNodesModification:
	"""
	This class applies the coordinate displacement errors to the set of nodes.
	The parameters of displacements are translated to all error controllers. 
	It could be directly (all nodes will have the same displacement) or 
	random values distributed by Gaussian with one sigma.
	"""
	def __init__(self):
		self.nodes = []
		self.error_controllers = []
		self.node_to_cntrl_dict = {}
		#---- these parameters can be interpreted as just values or sigmas for 
		#---- for Gaussian distributions
		self.param_dict = {"dx":0.,"dxp":0.,"dy":0.,"dyp":0.,"dz":0.,"dE":0.}
		
	def getErrorControllers(self):
		"""
		Returns the array of error controllers.
		"""
		return self.error_controllers
		
	def getCoordinateDisplacementParameters(self):
		"""
		Returns tuple with coordinate shift parameters. 
		"""
		dx = self.param_dict["dx"]
		dxp = self.param_dict["dxp"]
		dy = self.param_dict["dy"]
		dyp = self.param_dict["dyp"]
		dz = self.param_dict["dz"]
		dE = self.param_dict["dE"]
		return (dx, dxp, dy, dyp, dz, dE)
		
	def addLatticeNodes(self,nodes, lattice = None):
		"""
		Adds the error controllers to the nodes
		"""
		for node in nodes:
			errCntrl = ErrorCntrlCoordDisplacement("ErrCntrl:"+node.getName())
			errCntrl.setLattice(lattice)
			errCntrl.setOneNodeParent(node)
			self.error_controllers.append(errCntrl)
			self.node_to_cntrl_dict[node] = errCntrl
			
	def setDisplacementParameter(self,key,value):
		"""
		Sets the error value for a particular coordinate for all nodes.
		The cooridinate is defined by key parameter.
		"""
		if(self.param_dict.has_key(key)):
			self.param_dict[key] = value
		else:
			msg = "Class CoordinateDisplacementNodesModification - key-value problem"
			msg = msg + os.linesep
			msg = msg + "Method setDisplacementParameter:"
			msg = msg + os.linesep
			msg = msg + "You are trying to set value for key=" + key
			msg = msg + os.linesep
			msg = msg + "Keys could be only = (dx, dxp, dy, dyp, dz, dE)"
			msg = msg + os.linesep	
			orbitFinalize(msg)				
			return
		for errCntrl in self.error_controllers:
			errCntrl.setDisplacementParameter(key,value)
			
	def setGaussDistributedDisplacementParameter(self,key,value):
		"""
		Sets the random generated error value for a particular coordinate for all nodes.
		The cooridinate is defined by key parameter.
		"""
		if(self.param_dict.has_key(key)):
			self.param_dict[key] = value
		else:
			msg = "Class CoordinateDisplacementNodesModification - key-value problem"
			msg = msg + os.linesep
			msg = msg + "Method setGaussDistributedDisplacementParameter:"
			msg = msg + os.linesep
			msg = msg + "You are trying to set value for key=" + key
			msg = msg + os.linesep
			msg = msg + "Keys could be only = (dx, dxp, dy, dyp, dz, dE)"
			msg = msg + os.linesep	
			orbitFinalize(msg)				
			return
			
		for errCntrl in self.error_controllers:
			value_tmp = random.gauss(0.,value)
			errCntrl.setDisplacementParameter(key,value_tmp)

		
			