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
import random

#---- we need MPI for Gaussian distribution errors to be sure the lattices 
#---- are the same across all relevant node (the same communicator)
import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

# import from orbit Python utilities
from orbit.utils import orbitFinalize
from orbit.utils   import NamedObject
from orbit.utils   import TypedObject

# import general accelerator elements and lattice
from orbit.lattice import AccNode
from orbit.py_linac.lattice import Quad
from orbit.py_linac.lattice import DCorrectorH, DCorrectorV
from orbit.py_linac.lattice import Bend

# import error controllers from orbit.py_linac.errors package
from orbit.py_linac.errors import ErrorCntrlCoordDisplacement
from orbit.py_linac.errors import ErrorCntrlBendField
from orbit.py_linac.errors import ErrorCntrlLongitudinalDisplacement
from orbit.py_linac.errors import ErrorCntrlStraightRotationX
from orbit.py_linac.errors import ErrorCntrlStraightRotationY
from orbit.py_linac.errors import ErrorCntrlStraightRotationZ


class ErrorForNodesModification(NamedObject,TypedObject):
	"""
	The base abstract class for set of separate nodes modification 
	with two error nodes: one at the entrance and one at the exit 
	of the lattice node. The lattice is specified for possible
	needs in the future when we will introduce errors for some
	section of the lattice as a whole. Right now we are working 
	only with nodes.
	"""
	def __init__(self, name = "no_name", type_in = "ErrorForNodesModification"):
		NamedObject.__init__(self, name)
		TypedObject.__init__(self, type_in)
		self.nodes = []
		self.error_controllers = []
		self.node_to_cntrl_dict = {}
		self.lattice = None
		
	def _getInstanceOfErrorController(self):
		"""
		This is abstract method. The subclasses should implement this method
		according the particular error type.
		"""
		return None
		
	def updateErrorParameters(self):
		"""
		This is abstract method. The sub-classes should implement this method
		according the particular error type. It updates the error defining 
		parameters for all registered error controllers.
		"""
		pass	
		
	def getLatticeNodes(self):
		"""
		Returns the list of lattice nodes to which we applied errors. 
		"""
		return self.nodes
		
	def getErrorControllers(self):
		"""
		Returns the array of error controllers.
		"""
		return self.error_controllers
		
	def addLatticeNodes(self,nodes, lattice = None):
		"""
		Adds the error controllers to the nodes
		"""
		self.lattice = lattice
		self.nodes += nodes
		for node in self.nodes:
			errCntrl = self._getInstanceOfErrorController()
			errCntrl.setName("ErrCntrl:" + errCntrl.getShortTypeName() + ":" + node.getName())
			errCntrl.setLattice(lattice)
			errCntrl.setOneNodeParent(node)
			self.error_controllers.append(errCntrl)
			self.node_to_cntrl_dict[node] = errCntrl
		self.updateErrorParameters()
		
class CoordinateDisplacementNodesModification(ErrorForNodesModification):
	"""
	This class applies the coordinate displacement errors to the set of nodes.
	The parameters of displacements are translated to all error controllers. 
	It could be fixed (all nodes will have the same displacement) or 
	random values distributed around 0. by Gaussian with one sigma.
	"""
	def __init__(self, name = "no_name"):
		ErrorForNodesModification.__init__(self,name,"CoordinateDisplacementNodesModification")
		#---- these parameters can be interpreted as just values or sigmas for 
		#---- for Gaussian distributions
		self.param_dict = {"dx":0.,"dxp":0.,"dy":0.,"dyp":0.,"dz":0.,"dE":0.}
		
	def _getInstanceOfErrorController(self):
		"""
		Returns the instance of ErrorCntrlCoordDisplacement error controller.
		"""
		return ErrorCntrlCoordDisplacement()
		
	def updateErrorParameters(self):
		"""
		This is an implementation of parent class abstract method. 
		It updates the error defining parameters for all registered 
		error controllers.
		"""
		for key in 	self.param_dict.keys():
			for errCntrl in self.error_controllers:
				errCntrl.setDisplacementParameter(key,self.param_dict[key])			

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
		self.updateErrorParameters()
			
	def setGaussDistributedDisplacementParameter(self,key,value,cut_off_level = 3.0, comm = mpi_comm.MPI_COMM_WORLD):
		"""
		Sets the random generated error value for a particular coordinate for all nodes.
		The cooridinate is defined by key parameter.
		"""
		value = abs(value)
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
		#---------------------------------------------------------------------
		for errCntrl in self.error_controllers:
			value_tmp = random.gauss(0.,value)
			while(abs(value_tmp) > abs(value)*cut_off_level):
				value_tmp = random.gauss(0.,value)
			main_rank = 0
			value_tmp = orbit_mpi.MPI_Bcast(value_tmp,mpi_datatype.MPI_DOUBLE,main_rank,comm)
			errCntrl.setDisplacementParameter(key,value_tmp)

class LongitudinalDisplacementNodesModification(ErrorForNodesModification):
	"""
	This class shifts lattice nodes in longitudinal directions.
	The shift length is sent to all error controllers. 
	It could be fixed (all nodes will have the same displacement) or 
	random values distributed around 0. by Gaussian with one sigma.
	"""
	def __init__(self, name = "no_name"):
		ErrorForNodesModification.__init__(self,name,"LongitudinalDisplacementNodesModification")
		#---- this parameter can be interpreted as just a value or a sigma for 
		#---- for Gaussian distributions
		self.shift_length = shift_length
		
	def _getInstanceOfErrorController(self):
		"""
		Returns the instance of ErrorCntrlCoordDisplacement error controller.
		"""
		return ErrorCntrlLongitudinalDisplacement()

	def updateErrorParameters(self):
		"""
		This is an implementation of parent class abstract method. 
		It updates the error defining parameters for all registered 
		error controllers.
		"""
		for errCntrl in self.error_controllers:
			errCntrl.setShiftLength(self.shift_length)

	def getShiftLength(self):
		"""
		Returns the shift-length parameter.
		"""
		return self.shift_length	

	def setShiftLength(self,shift_length):
		"""
		Sets the shift-length parameter.
		"""
		self.shift_length = shift_length
		self.updateErrorParameters()

	def setGaussDistributedShiftLength(self,shift_length,cut_off_level = 3.0, comm = mpi_comm.MPI_COMM_WORLD):
		"""
		Sets the random generated error shift_length for all nodes.
		"""
		for errCntrl in self.error_controllers:
			shift_length_tmp = random.gauss(0.,shift_length)
			while(abs(shift_length_tmp) > abs(shift_length)*cut_off_level):
				shift_length_tmp = random.gauss(0.,shift_length)
			main_rank = 0
			shift_length_tmp = orbit_mpi.MPI_Bcast(shift_length,mpi_datatype.MPI_DOUBLE,main_rank,comm)
			errCntrl.setShiftLength(shift_length_tmp)
			
class StraightRotationZ_NodesModification(ErrorForNodesModification):
	"""
	This class rotate lattice nodes around z-axis directions.
	The rotating angle is sent to all error controllers. 
	It could be fixed (all nodes will have the same angle) or 
	random values distributed around 0. by Gaussian with one sigma.
	"""
	def __init__(self, name = "no_name"):
		ErrorForNodesModification.__init__(self,name,"StraightRotationZ_NodesModification")
		#---- this parameter can be interpreted as just a value or a sigma for 
		#---- for Gaussian distributions
		self.angle = 0.
		
	def _getInstanceOfErrorController(self):
		"""
		Returns the instance of ErrorCntrlStraightRotationZ error controller.
		"""
		return ErrorCntrlStraightRotationZ()

	def updateErrorParameters(self):
		"""
		This is an implementation of parent class abstract method. 
		It updates the error defining parameters for all registered 
		error controllers.
		"""
		for errCntrl in self.error_controllers:
			errCntrl.setRotationAngle(self.angle)

	def getRotationAngle(self):
		"""
		Returns the angle parameter.
		"""
		return self.angle	

	def setRotationAngle(self,angle):
		"""
		Sets the angle parameter.
		"""
		self.angle = angle
		self.updateErrorParameters()

	def setGaussDistributedAngle(self,angle,cut_off_level = 3.0, comm = mpi_comm.MPI_COMM_WORLD):
		"""
		Sets the random generated error angle for all nodes.
		"""
		for errCntrl in self.error_controllers:
			angle_tmp = random.gauss(0.,angle)
			while(abs(angle_tmp) > abs(angle)*cut_off_level):
				angle_tmp = random.gauss(0.,angle)
			main_rank = 0
			angle_tmp = orbit_mpi.MPI_Bcast(angle,mpi_datatype.MPI_DOUBLE,main_rank,comm)
			errCntrl.setRotationAngle(angle_tmp)			
			
class StraightRotationX_NodesModification(ErrorForNodesModification):
	"""
	This class rotate lattice nodes around x-axis directions.
	The rotating angle is sent to all error controllers. 
	It could be fixed (all nodes will have the same angle) or 
	random values distributed around 0. by Gaussian with one sigma.
	"""
	def __init__(self, name = "no_name"):
		ErrorForNodesModification.__init__(self,name,"StraightRotationX_NodesModification")
		#---- this parameter can be interpreted as just a value or a sigma for 
		#---- for Gaussian distributions
		self.angle = 0.
		
	def _getInstanceOfErrorController(self):
		"""
		Returns the instance of ErrorCntrlStraightRotationZ error controller.
		"""
		return ErrorCntrlStraightRotationX()
		
	def updateErrorParameters(self):
		"""
		This is an implementation of parent class abstract method. 
		It updates the error defining parameters for all registered 
		error controllers.
		"""
		for errCntrl in self.error_controllers:
			entr_parent_node = errCntrl.getEntanceNodeParent()
			exit_parent_node = errCntrl.getExitNodeParent()
			if(entr_parent_node != exit_parent_node):
				msg = "Class StraightRotationX_NodesModification update parameters problem!"
				msg = msg + os.linesep
				msg = msg + "For this type of Error Node Controllers parent node is one, not two!"
				msg = msg + os.linesep
				msg = msg + "Entrance parent node="+entr_parent_node.getName()
				msg = msg + os.linesep
				msg = msg + "Exit     parent node="+exit_parent_node.getName()
				msg = msg + os.linesep
				msg = msg + "Stop"
				msg = msg + os.linesep
				orbitFinalize(msg)
				return
			errCntrl.setRotationAngle(self.angle)
			errCntrl.setBaseLength(entr_parent_node.getLength())

	def getRotationAngle(self):
		"""
		Returns the angle parameter.
		"""
		return self.angle	

	def setRotationAngle(self,angle):
		"""
		Sets the angle parameter.
		"""
		self.angle = angle
		self.updateErrorParameters()

	def setGaussDistributedAngle(self,angle,cut_off_level = 3.0, comm = mpi_comm.MPI_COMM_WORLD):
		"""
		Sets the random generated error angle for all nodes.
		"""
		for errCntrl in self.error_controllers:
			angle_tmp = random.gauss(0.,angle)
			while(abs(angle_tmp) > abs(angle)*cut_off_level):
				angle_tmp = random.gauss(0.,angle)
			main_rank = 0
			angle_tmp = orbit_mpi.MPI_Bcast(angle,mpi_datatype.MPI_DOUBLE,main_rank,comm)
			errCntrl.setRotationAngle(angle_tmp)
			
class StraightRotationY_NodesModification(ErrorForNodesModification):
	"""
	This class rotate lattice nodes around y-axis directions.
	The rotating angle is sent to all error controllers. 
	It could be fixed (all nodes will have the same angle) or 
	random values distributed around 0. by Gaussian with one sigma.
	"""
	def __init__(self, name = "no_name"):
		ErrorForNodesModification.__init__(self,name,"StraightRotationY_NodesModification")
		#---- this parameter can be interpreted as just a value or a sigma for 
		#---- for Gaussian distributions
		self.angle = 0.
		
	def _getInstanceOfErrorController(self):
		"""
		Returns the instance of ErrorCntrlStraightRotationZ error controller.
		"""
		return ErrorCntrlStraightRotationY()
		
	def updateErrorParameters(self):
		"""
		This is an implementation of parent class abstract method. 
		It updates the error defining parameters for all registered 
		error controllers.
		"""
		for errCntrl in self.error_controllers:
			entr_parent_node = errCntrl.getEntanceNodeParent()
			exit_parent_node = errCntrl.getExitNodeParent()
			if(entr_parent_node != exit_parent_node):
				msg = "Class StraightRotationY_NodesModification update parameters problem!"
				msg = msg + os.linesep
				msg = msg + "For this type of Error Node Controllers parent node is one, not two!"
				msg = msg + os.linesep
				msg = msg + "Entrance parent node="+entr_parent_node.getName()
				msg = msg + os.linesep
				msg = msg + "Exit     parent node="+exit_parent_node.getName()
				msg = msg + os.linesep
				msg = msg + "Stop"
				msg = msg + os.linesep
				orbitFinalize(msg)
				return
			errCntrl.setRotationAngle(self.angle)
			errCntrl.setBaseLength(entr_parent_node.getLength())

	def getRotationAngle(self):
		"""
		Returns the angle parameter.
		"""
		return self.angle	

	def setRotationAngle(self,angle):
		"""
		Sets the angle parameter.
		"""
		self.angle = angle
		self.updateErrorParameters()

	def setGaussDistributedAngle(self,angle,cut_off_level = 3.0,comm = mpi_comm.MPI_COMM_WORLD):
		"""
		Sets the random generated error angle for all nodes.
		"""
		for errCntrl in self.error_controllers:
			angle_tmp = random.gauss(0.,angle)
			while(abs(angle_tmp) > abs(angle)*cut_off_level):
				angle_tmp = random.gauss(0.,angle)
			main_rank = 0
			angle_tmp = orbit_mpi.MPI_Bcast(angle,mpi_datatype.MPI_DOUBLE,main_rank,comm)
			errCntrl.setRotationAngle(angle_tmp)
			
#------------------------------------------------------------------
# The magnet field errors application classes
#------------------------------------------------------------------

class QuadFieldsErrorsDeployment(NamedObject,TypedObject):
	"""
	Class will apply the errors to the fields of the quads
	"""
	def __init__(self, name = "no_name", type_in = "QuadFieldsErrorsDeployment"):
		NamedObject.__init__(self, name)
		TypedObject.__init__(self, type_in)
		#---- self.quad_and_field_arr[ind] = [[quad,field_init],...]
		self.quad_and_field_arr = []
		
	def addQuads(self,quads):
		"""
		Add quads to the inner array of quads.
		"""
		for quad in quads:
			if(isinstance(quad,Quad)):
				self.quad_and_field_arr.append([quad,quad.getParam("dB/dr")])
			
	def addQuad(self,quad):
		"""
		Add one quad to the inner array of quads.
		"""	
		if(isinstance(quad,Quad)):
			self.quad_and_field_arr.append([quad,quad.getParam("dB/dr")])
		
	def restoreFields(self):
		for [quad,field_init] in self.quad_and_field_arr:
			quad.setParam("dB/dr",field_init)
	
	def setGaussDistributedRealtiveErrors(self,relative_error,cut_off_level = 3.0,comm = mpi_comm.MPI_COMM_WORLD):
		"""
		Sets the random generated error field for all quads.
		"""
		for [quad,field_init] in self.quad_and_field_arr:
			rel_err = random.gauss(0.,relative_error)
			while(abs(rel_err) > abs(relative_error)*cut_off_level):
				rel_err = random.gauss(0.,relative_error)
			main_rank = 0
			rel_err = orbit_mpi.MPI_Bcast(rel_err,mpi_datatype.MPI_DOUBLE,main_rank,comm)
			field = field_init*(1.0 + rel_err)
			quad.setParam("dB/dr",field)
			
class BendFieldNodesModification(ErrorForNodesModification):
	"""
	This class will apply the errors to the fields of the bends using energy shift
	"""
	def __init__(self, name = "no_name"):
		ErrorForNodesModification.__init__(self,name,"BendFieldNodesModification")
		#---- defines a relative change in the bend field
		self.relative_field_change = 0.
		
	def _getInstanceOfErrorController(self):
		"""
		Returns the instance of ErrorCntrlCoordDisplacement error controller.
		"""
		return ErrorCntrlBendField()
		
	def updateErrorParameters(self):
		"""
		This is an implementation of parent class abstract method. 
		It updates the realtive field change parameter for all registered 
		error controllers.
		"""
		for errCntrl in self.error_controllers:
			errCntrl.setRelativeFieldChange(self.relative_field_change)
			
	def addBends(self,bends):
		"""
		Adds up bend nodes array.
		"""
		self.addLatticeNodes(bends)
		
	def addBend(self,bend):
		"""
		Adds up bend node.
		"""
		self.addLatticeNodes([bend,])		
		
	def setRelativeFieldChange(self,relative_field_change):
		"""
		Sets the realtive change in the bend field.
		"""
		self.relative_field_change = relative_field_change
		self.updateErrorParameters()
		
	def getRelativeFieldChange(self):
		"""
		Returns the realtive change in the bend field.
		"""
		return self.relative_field_change

	def setGaussDistributedRelativeFieldError(self,relative_error,cut_off_level = 3.0, comm = mpi_comm.MPI_COMM_WORLD):
		"""
		Sets the random generated field error for all bends. The same value for all registered bends.
		"""
		rel_err = random.gauss(0.,relative_error)
		while(abs(rel_err) > abs(relative_error)*cut_off_level):
			rel_err = random.gauss(0.,relative_error)
		main_rank = 0
		rel_err = orbit_mpi.MPI_Bcast(rel_err,mpi_datatype.MPI_DOUBLE,main_rank,comm)
		self.relative_field_change = rel_err
		self.updateErrorParameters()
