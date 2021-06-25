"""
This module is a collection of Error Nodes and Controller Classes. Controller are not
subclasses of AccNode, and they are not part of the lattice. They resemble the 
Linac Cavities as controlling structure for RF gaps. The gaps are real nodes 
in the lattice, and Cavities perform synchronisation of time arrivals and 
field amplitudes in RF gaps. The Error Controller Classes keep references for 
two error nodes: first at the entrance of error applied part of the lattice (or node)
and the second at the exit. For example, if we apply offset to a part of the lattice,
we will shift particles' coordinates transversly by some amount at the beginning
of this part and back at the exit.

The elementary operations are performed by functions from error_base 
implemented in C++.

/scr/orbit/Errors/errorbase.cc

"""

import os
import math
import sys

# import the finalization function 
from orbit.utils import orbitFinalize
from orbit.utils   import NamedObject
from orbit.utils   import TypedObject
from orbit.utils   import ParamsDictObject

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

# import error packages from C++
from error_base import CoordDisplacement
from error_base import LongDisplacement
from error_base import StraightRotationXY
from error_base import StraightRotationXSI
from error_base import StraightRotationXSF
from error_base import StraightRotationYSI
from error_base import StraightRotationYSF
from error_base import BendFieldI
from error_base import BendFieldF
from error_base import BendDisplacementXI
from error_base import BendDisplacementXF
from error_base import BendDisplacementYI
from error_base import BendDisplacementYF
from error_base import BendDisplacementLI
from error_base import BendDisplacementLF
from error_base import RotationI
from error_base import RotationF
from error_base import DipoleKickerOsc
from error_base import QuadKicker
from error_base import QuadKickerOsc


#------------------------------------------------------------------------
#     Error Nodes classes. All of them should be sub-classes
#     of AccErrorNode class. User will not create these nodes
#     directly, They will created by sub-classes of BaseErrorController
#     class.
#------------------------------------------------------------------------

class AccErrorNode(AccNodeBunchTracker):
	"""
	This is an abstract base class for Error nodes hierarchy. The sub-classes
	should implement track(self, paramsDict) method. Usually it will be 
	call to /scr/orbit/Errors/errorbase.cc functions with parameters 
	provided by the controller class through 
	self.error_controller_params_func call.
	"""
	def __init__(self,name = "no_name", type_in = "Base_Error_Controller"):
		AccNodeBunchTracker.__init__(self,name,type_in)
		self.error_controller_params_func = None
	
	def setErrorControllerParamFunc(self,error_controller_params_func):
		"""
		Sets the callable object error_controller_params_func from the
		error node controller.
		"""
		self.error_controller_params_func = error_controller_params_func

	def track(self, paramsDict):
		"""
		It is tracking the bunch through the node. Each subclass 
		should implement this method. errorCntrl will provide
		parameters for tracking.
		"""
		error_controller_params_func = self.error_controller_params_func
		pass	
	
	def trackDesign(self, bunch, paramsDict = {}, actionContainer = None):
		"""
		It tracks the design bunch through the AccErrorNode instance - does nothing.
		The existence of this method allows to use this class in the Linac type
		lattice in addition to the Teapot lattices.
		"""	
		pass

class ErrorLongitudinalDisplacementNode(AccErrorNode):
	"""
	The subclass of AccErrorNode. It can shift the bunch longitudinally by using a drift
	TEAPOT function. This node is used to simulate longitudinal shift of lattice
	elements. The shift length is provided by the controller's 
	error_controller_params_func. This function return shift_length for the entrance
	node and -shift_length for the exit one.
	"""
	def __init__(self, name = "no_name"):
		AccErrorNode.__init__(self,name,"ErrorLongitudinalDisplacementNode")
		
	def track(self, paramsDict):
		"""
		Performs shift_length TEAPOT drift tracking.
		"""
		bunch = paramsDict["bunch"]
		shift_length = self.error_controller_params_func()
		LongDisplacement(bunc,shift_length)
		
class ErrorCoordDisplacementNode(AccErrorNode):
	"""
	The subclass of AccErrorNode. It can shift any of 6 coordinates by specified 
	value. The parameters are provided by the controller's error_controller_params_func
	callable method reference.
	The shift could be performed for (dx, dxp, dy, dyp, dz, dE), default values are 0.
	"""
	def __init__(self, name = "no_name"):
		AccErrorNode.__init__(self,name,"ErrorCoordDisplacementNode")
		
	def track(self, paramsDict):
		"""
		Performs (dx, dxp, dy, dyp, dz, dE) shift with CoordDisplacement function.
		"""
		bunch = paramsDict["bunch"]
		(dx, dxp, dy, dyp, dz, dE) = self.error_controller_params_func()
		CoordDisplacement(bunch, dx, dxp, dy, dyp, dz, dE)		

class ErrorBendFieldNode(AccErrorNode):
	"""
	The subclass of AccErrorNode. It changes of the energy of each particle 
	in the bunch to emulate the relative change of the magnetic field of the
	bend. deltaB/B = deltaP/P = (1/beta)^2*deltaE/E
	"""
	def __init__(self, name = "no_name"):
		AccErrorNode.__init__(self,name,"ErrorBendFieldNode")
		
	def track(self, paramsDict):
		"""
		Performs deltaE energy shift with CoordDisplacement function.
		"""
		bunch = paramsDict["bunch"]
		syncPart = bunch.getSyncParticle()
		e_total = (syncPart.kinEnergy() + syncPart.mass())
		beta = syncPart.beta()
		relative_field_change = self.error_controller_params_func()
		deltaE = relative_field_change*beta**2*e_total
		CoordDisplacement(bunch, 0., 0., 0., 0., 0., deltaE)

class ErrorStraightRotationZNode(AccErrorNode):
	"""
	The subclass of AccErrorNode. It describes the rotation around z-axis 
	by a certain angle. The angle parameter is provided by the controller's 
	error_controller_params_func. This function return angle for the entrance
	node and -angle for the exit one.
	"""
	def __init__(self, name = "no_name"):
		AccErrorNode.__init__(self,name,"ErrorStraightRotationZNode")
		
	def track(self, paramsDict):
		"""
		Performs rotation by TEAPOT error_base StraightRotationZ function.
		"""
		bunch = paramsDict["bunch"]
		angle = self.error_controller_params_func()	
		StraightRotationXY(bunch, angle)

class ErrorStraightRotationXNode(AccErrorNode):
	"""
	The subclass of AccErrorNode. It describes the rotation around x-axis 
	by a certain angle. The angle parameter is provided by the controller's 
	error_controller_params_func. This function return angle,length, and direction
	of the node/nodes that are rotated for the entrance node and for the exit one.
	The rotation is performed by StraightRotationYSI and StraightRotationYSF
	error_base TEAPOT functions for the entrance (Init) and exit (Final) 
	error nodes.
	"""
	def __init__(self, name = "no_name"):
		AccErrorNode.__init__(self,name,"ErrorStraightRotationXNode")
		
	def track(self, paramsDict):
		"""
		Performs rotation by TEAPOT error_base StraightRotationYSI and
		StraightRotationYSF functions.
		"""
		bunch = paramsDict["bunch"]
		(angle,length,direction) = self.error_controller_params_func()
		if(direction > 0):
			StraightRotationYSI(bunch,angle,length)
			return
		if(direction < 0):
			StraightRotationYSF(bunch,angle,length)
			return
		msg = "ErrorStraightRotationXNode - no direction parameter for the rotation!"
		msg = msg + os.linesep
		msg = msg + "direction="+str(direction)
		msg = msg + os.linesep
		msg = msg + "Stop."
		msg = msg + os.linesep
		orbitFinalize(msg)	

class ErrorStraightRotationYNode(AccErrorNode):
	"""
	The subclass of AccErrorNode. It describes the rotation around y-axis 
	by a certain angle. The angle parameter is provided by the controller's 
	error_controller_params_func. This function return angle,length, and direction
	of the node/nodes that are rotated for the entrance node and for the exit one.
	The rotation is performed by StraightRotationXSI and StraightRotationXSF
	error_base TEAPOT functions for the entrance (Init) and exit (Final) 
	error nodes.
	"""
	def __init__(self, name = "no_name"):
		AccErrorNode.__init__(self,name,"ErrorStraightRotationYNode")
		
	def track(self, paramsDict):
		"""
		Performs rotation by TEAPOT error_base StraightRotationXSI and
		StraightRotationXSF functions.
		"""
		bunch = paramsDict["bunch"]
		(angle,length,direction) = self.error_controller_params_func()
		if(direction > 0):
			StraightRotationXSI(bunch,angle,length)
			return
		if(direction < 0):
			StraightRotationXSF(bunch,angle,length)
			return
		msg = "ErrorStraightRotationYNode - no direction parameter for the rotation!"
		msg = msg + os.linesep
		msg = msg + "direction="+str(direction)
		msg = msg + os.linesep
		msg = msg + "Stop."
		msg = msg + os.linesep
		orbitFinalize(msg)			
		
#------------------------------------------------------------------------
#     Error Control classes. All of them should be sub-classes
#     of BaseErrorController class. The Error Controller Classes keep 
#     references for two error nodes: first at the entrance of error 
#     applied part of the lattice (or node) and the second at the exit.
#     The Error Control classes implement two functions
#     getEntanceNodeParameters() and getExitNodeParameters() that return
#     parameters for Error Nodes. The Error Control classes also keep
#     the reference to the parent Node/Nodes (real lattice nodes - like 
#     quads,bends, RF gaps etc.) of the Error nodes. 
#------------------------------------------------------------------------

class BaseErrorController(NamedObject, TypedObject, ParamsDictObject):
	"""
	The base class for Error Controllers. It keeps the reference to start and stop
	lattice nodes. The transformation function is defined by  self.entranceErrorAccNode,
	and backward transformation is defined by  self.exitErrorAccNode which are set up as
	children of two (or just one) entrance and exit lattice nodes.
	"""
	def __init__(self, name = "no_name", type_in = "Base_Error_Controller"):
		NamedObject.__init__(self, name)
		TypedObject.__init__(self, type_in)
		ParamsDictObject.__init__(self)
		#--------------------------------------
		#---- entrance and exit Error AccNodes
		self.short_type_name = "None"
		self.entranceErrorAccNode = None
		self.exitErrorAccNode = None
		self.entranceErrorAccNodeParent = None
		self.exitErrorAccNodeParent = None
		self.accLattice = None
		
	def getShortTypeName(self):
		"""
		Returns the short type name that will be included into the error node name.
		"""
		return self.short_type_name
		
	def getEntanceNodeParameters(self):
		"""
		Returns tuple with parameters for entranceAccNode.
		This method should be implemented by sub-classes.
		"""
		pass 
	
	def getExitNodeParameters(self):
		"""
		Returns tuple with parameters for exitAccNode.
		This method should be implemented by sub-classes.
		"""
		pass 
		
	def setLattice(self,accLattice):
		"""
		Sets the lattice.
		"""
		self.accLattice = accLattice
		
	def getLattice(self):
		"""
		Returns the reference to the lattice.
		"""
		return self.accLattice
		
	def getEntanceNode(self):
		"""
		Returns the lattice entrance (start) node for errors application.
		"""
		return self.entranceErrorAccNode
		
	def getExitNode(self):
		"""
		Returns the lattice exit (stop) node for errors application.
		"""
		return self.exitErrorAccNode
		
	def setEntanceNodeParent(self,entranceAccNodeParent):
		"""
		Set up the entrance (start) lattice node for error effects.
		"""
		self.entranceErrorAccNodeParent = entranceAccNodeParent
		self.entranceErrorAccNode.setName("ErrNode:"+self.short_type_name+":Entr:"+entranceAccNodeParent.getName())
		self.entranceErrorAccNodeParent.getChildNodes(AccNode.ENTRANCE).insert(0,self.entranceErrorAccNode)
		
	def setExitNodeParent(self,exitAccNodeParent):
		"""
		Set up the exit (stop) lattice node for error effects.
		"""
		self.exitErrorAccNodeParent = exitAccNodeParent
		self.exitErrorAccNode.setName("ErrNode:"+self.short_type_name+":Exit:"+exitAccNodeParent.getName())
		self.exitErrorAccNodeParent.addChildNode(self.exitErrorAccNode,AccNode.EXIT)
		
	def setOneNodeParent(self,accNodeParent):
		"""
		If we apply error transformation only for one lattice node.		
		"""
		self.setEntanceNodeParent(accNodeParent)
		self.setExitNodeParent(accNodeParent)
		
	def getEntanceNodeParent(self):
		"""
		Returns the entrance (start) lattice node for error effects.
		"""
		return self.entranceErrorAccNodeParent
		
	def getExitNodeParent(self):
		"""
		Returns the exit (stop) lattice node for error effects.
		"""
		return self.exitErrorAccNodeParent
		
	def cleanParentNodes(self):
		"""
		Removes children representing error nodes from the parent node or nodes.
		"""
		if(self.entranceErrorAccNodeParent != None and self.entranceErrorAccNode != None):
			self.entranceErrorAccNodeParent.getChildNodes(AccNode.ENTRANCE).remove(self.entranceErrorAccNode)
		if(self.exitErrorAccNodeParent != None and self.exitErrorAccNode != None):
			self.exitErrorAccNodeParent.getChildNodes(AccNode.EXIT).remove(self.exitErrorAccNode)

class ErrorCntrlLongitudinalDisplacement(BaseErrorController):
	"""
	Subclass of BaseErrorController. It defines the error controller for 
	longitudinal shift Error node. 
	"""
	def __init__(self,name = "no_name"):
		BaseErrorController.__init__(self,name,"ErrorCntrlLongitudinalDisplacement")
		self.short_type_name = "LongDisp"
		self.entranceErrorAccNode = ErrorLongitudinalDisplacementNode()
		self.exitErrorAccNode = ErrorLongitudinalDisplacementNode()
		self.entranceErrorAccNode.setErrorControllerParamFunc(self.getEntanceNodeShiftLength)
		self.exitErrorAccNode.setErrorControllerParamFunc(self.getExitNodeShiftLength)
		self.shift_length = 0.
		
	def getEntanceNodeShiftLength(self):
		"""
		Returns shift-length parameter for entrance Error node
		"""
		return self.shift_length
		
	def getExitNodeShiftLength(self):
		"""
		Returns shift-length parameter for exit Error node
		"""
		return self.shift_length		
		
	def setShiftLength(self,shift_length):
		"""
		Sets the shift-length parameter.
		"""
		self.shift_length = shift_length
		
	def getShiftLength(self):
		"""
		Returns the shift-length parameter.
		"""
		return self.shift_length	

class ErrorCntrlCoordDisplacement(BaseErrorController):
	"""
	Subclass of BaseErrorController. It defines the error controller for 
	shift transformations. User can shift any of 6 coordiantes 
	(dx, dxp, dy, dyp, dz, dE). The values of shfing parameters are defined by 
	setDisplacementParameter(self,key,value) method with "dx","dxp","dy","dyp",
	"dz","dE" keys.
	"""
	def __init__(self,name = "no_name"):
		BaseErrorController.__init__(self,name,"ErrorCntrlCoordDisplacement")
		self.short_type_name = "CoordDisp"
		self.entranceErrorAccNode = ErrorCoordDisplacementNode()
		self.exitErrorAccNode = ErrorCoordDisplacementNode()
		self.entranceErrorAccNode.setErrorControllerParamFunc(self.getEntanceNodeParameters)
		self.exitErrorAccNode.setErrorControllerParamFunc(self.getExitNodeParameters)
		self.param_dict = {"dx":0.,"dxp":0.,"dy":0.,"dyp":0.,"dz":0.,"dE":0.}
		
	def getEntanceNodeParameters(self):
		"""
		Returns tuple with parameters for ErrorCoordDisplacementNode tracking
		at the entrance.
		"""
		dx = self.param_dict["dx"]
		dxp = self.param_dict["dxp"]
		dy = self.param_dict["dy"]
		dyp = self.param_dict["dyp"]
		dz = self.param_dict["dz"]
		dE = self.param_dict["dE"]
		return (dx, dxp, dy, dyp, dz, dE)
	
	def getExitNodeParameters(self):
		"""
		Returns tuple with parameters for ErrorCoordDisplacementNode tracking
		at the exit.
		"""
		dx = self.param_dict["dx"]
		dxp = self.param_dict["dxp"]
		dy = self.param_dict["dy"]
		dyp = self.param_dict["dyp"]
		dz = self.param_dict["dz"]
		dE = self.param_dict["dE"]
		return (-dx, -dxp, -dy, -dyp, -dz, -dE)
		
	def setDisplacementParameter(self,key,value):
		"""
		Sets one of the parameters for keys "dx","dxp","dy","dyp","dz","dE".
		"""
		self.param_dict[key] = value
		
	def getDisplacementParameters(self):
		"""
		Returns the (dx, dxp, dy, dyp, dz, dE) parameters.
		"""
		return self.getEntanceNodeParameters()

class ErrorCntrlBendField(BaseErrorController):
	"""
	Subclass of BaseErrorController. It defines the error controller for 
	bend field. User specifies the relative change in the bend magnetic
	field.
	"""
	def __init__(self,name = "no_name"):
		BaseErrorController.__init__(self,name,"ErrorCntrlBendField")
		self.short_type_name = "BendFieldNoise"
		self.entranceErrorAccNode = ErrorBendFieldNode()
		self.exitErrorAccNode = ErrorBendFieldNode()
		self.entranceErrorAccNode.setErrorControllerParamFunc(self.getEntanceNodeParameters)
		self.exitErrorAccNode.setErrorControllerParamFunc(self.getExitNodeParameters)
		self.relative_field_change = 0.
		
	def getEntanceNodeParameters(self):
		"""
		Returns the realtive change in the bend field at the entrance node.
		"""
		return self.relative_field_change
	
	def getExitNodeParameters(self):
		"""
		Returns the realtive change in the bend field at the entrance node.
		"""
		return -self.relative_field_change
		
	def setRelativeFieldChange(self,relative_field_change):
		"""
		Sets the realtive change in the bend field.
		"""
		self.relative_field_change = relative_field_change
		
	def getRelativeFieldChange(self):
		"""
		Returns the realtive change in the bend field.
		"""
		return self.getEntanceNodeParameters()

class ErrorCntrlStraightRotationZ(BaseErrorController):
	"""
	Subclass of BaseErrorController. It defines the error controller for 
	the rotation around z-axis Error node. 
	"""
	def __init__(self,name = "no_name"):
		BaseErrorController.__init__(self,name,"ErrorCntrlStraightRotationZ")
		self.short_type_name = "RotZ"
		self.entranceErrorAccNode = ErrorStraightRotationZNode()
		self.exitErrorAccNode = ErrorStraightRotationZNode()
		self.entranceErrorAccNode.setErrorControllerParamFunc(self.getEntanceNodeRotationAngle)
		self.exitErrorAccNode.setErrorControllerParamFunc(self.getExitNodeRotationAngle)
		self.angle = 0.
		
	def getEntanceNodeRotationAngle(self):
		"""
		Returns angle parameter for entrance Error node
		"""
		return self.angle
		
	def getExitNodeRotationAngle(self):
		"""
		Returns angle parameter for exit Error node
		"""
		return -self.angle

	def setRotationAngle(self,angle):
		"""
		Sets the angle parameter.
		"""
		self.angle = angle
		
	def getRotationAngle(self):
		"""
		Returns the angle parameter.
		"""
		return self.angle			

class ErrorCntrlStraightRotationX(BaseErrorController):
	"""
	Subclass of BaseErrorController. It defines the error controller for 
	the rotation around x-axis Error node. 
	"""
	def __init__(self,name = "no_name"):
		BaseErrorController.__init__(self,name,"ErrorCntrlStraightRotationX")
		self.short_type_name = "RotX"
		self.entranceErrorAccNode = ErrorStraightRotationXNode()
		self.exitErrorAccNode = ErrorStraightRotationXNode()
		self.entranceErrorAccNode.setErrorControllerParamFunc(self.getEntanceNodeRotationParams)
		self.exitErrorAccNode.setErrorControllerParamFunc(self.getExitNodeRotationParams)
		self.angle = 0.
		self.base_length = 0.
		
	def getEntanceNodeRotationParams(self):
		"""
		Returns angle, base_length, and direction parameters for entrance Error node
		"""
		return (self.angle,self.base_length,+1)
		
	def getExitNodeRotationParams(self):
		"""
		Returns angle, base_length, and direction parameters for exit Error node
		"""
		return (self.angle,self.base_length,-1)

	def setRotationAngle(self,angle):
		"""
		Sets the angle parameter.
		"""
		self.angle = angle
		
	def getRotationAngle(self):
		"""
		Returns the angle parameter.
		"""
		return self.angle	
		
	def setBaseLength(self,base_length):
		"""
		Sets the base base_length parameter.
		"""
		self.base_length = base_length
		
	def getBaseLength(self):
		"""
		Returns the base base_length parameter.
		"""
		return self.base_length		

class ErrorCntrlStraightRotationY(BaseErrorController):
	"""
	Subclass of BaseErrorController. It defines the error controller for 
	the rotation around y-axis Error node. 
	"""
	def __init__(self,name = "no_name"):
		BaseErrorController.__init__(self,name,"ErrorCntrlStraightRotationY")
		self.short_type_name = "RotY"
		self.entranceErrorAccNode = ErrorStraightRotationYNode()
		self.exitErrorAccNode = ErrorStraightRotationYNode()
		self.entranceErrorAccNode.setErrorControllerParamFunc(self.getEntanceNodeRotationParams)
		self.exitErrorAccNode.setErrorControllerParamFunc(self.getExitNodeRotationParams)
		self.angle = 0.
		self.base_length = 0.
		
	def getEntanceNodeRotationParams(self):
		"""
		Returns angle, base_length, and direction parameters for entrance Error node
		"""
		return (self.angle,self.base_length,+1)
		
	def getExitNodeRotationParams(self):
		"""
		Returns angle, base_length, and direction parameters for exit Error node
		"""
		return (self.angle,self.base_length,-1)

	def setRotationAngle(self,angle):
		"""
		Sets the angle parameter.
		"""
		self.angle = angle
		
	def getRotationAngle(self):
		"""
		Returns the angle parameter.
		"""
		return self.angle	
		
	def setBaseLength(self,base_length):
		"""
		Sets the base base_length parameter.
		"""
		self.base_length = base_length
		
	def getBaseLength(self):
		"""
		Returns the base base_length parameter.
		"""
		return self.base_length	

