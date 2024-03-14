"""
This module is a collection of the linac accelerator nodes which are the subclasses of 
the AccNode class. We cannot use here TEAPOT nodes from the TEAPOT package 
directly because we use the field as a parameter for quads and dipole correctors
instead of k1 = 1/(B*rho)*(dB/dr).
The abstract AbstractRF_Gap class is a parent class for all RF gap model classes.
"""

import os
import math

# import the finalization function 
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

# Import the linac specific tracking from linac_tracking. This module has
# the following functions duplicated the original TEAPOT functions
# drift - linac drift tracking
# quad1 - linac quad linear part of tracking
# quad2 - linac quad non-linear part of tracking
import linac_tracking

class BaseLinacNode(AccNodeBunchTracker):
	""" 
	The base abstract class of the linac accelerator elements hierarchy. 
	It cannot be tilted. The direct subclasses of this class will be markers, 
	user defined nodes etc.
	"""
	def __init__(self, name = "none"):
		"""
		Constructor. Creates the base linac element. This is a superclass for all linac elements.
		"""
		AccNodeBunchTracker.__init__(self,name)
		self.setType("baseLinacNode")
		self.setParam("pos",0.)
		self.__linacSeqence = None
		# by default we use the TEAPOT tracker module
		self.tracking_module = TPB
		
	def setLinacTracker(self, switch = True):
		"""
		This method will switch tracker module to the linac specific traker by default
		"""
		if(switch):
			self.tracking_module = linac_tracking
		else:
			self.tracking_module = TPB		
		
	def isRFGap(self):
		"""
		Returns False. The RF Gap node returns True.
		"""
		return False
		
	def setSequence(self, seq):
		"""
		Sets the seqence.
		"""
		self.__linacSeqence = seq
		
	def setPosition(self,pos):
		"""
		Sets the position of the node inside the sequence. 
		If the node has non-zero length, the position is usually 
		at the center.
		"""
		self.setParam("pos",pos)
		
	def getPosition(self):
		"""
		Returns the position of the node inside the sequence. 
		If the node has non-zero length, the position is usually 
		at the center.
		"""		
		return self.getParam("pos")
		
	def getSequence(self):
		"""
		Returns the seqence.
		"""
		return self.__linacSeqence

	def trackDesign(self, paramsDict):
		"""
		The RF First Gap nodes will reload this method to setup the design time of passage 
		of the bunch through this node
		"""
		self.track(paramsDict)
		
	def trackDesignBunch(self, bunch, paramsDict = {}, actionContainer = None):
		"""
		It tracks the bunch through the AccNodeBunchTracker instance.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Design Bunch Tracking")
		paramsDict["bunch"] = bunch
		
		def trackDesign(paramsDict):
			node = paramsDict["node"]
			node.trackDesign(paramsDict)
			
		actionContainer.addAction(trackDesign, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict)
		actionContainer.removeAction(trackDesign, AccActionsContainer.BODY)				
		

class MarkerLinacNode(BaseLinacNode):
	"""
	This is a marker. It does nothing. If the user wants to perform operations with bunch 
	he/shi should specify tracking and design tracking functions.
	"""
	def __init__(self, name = "none"):
		BaseLinacNode.__init__(self,name)
		self.setType("markerLinacNode")


class LinacNode(BaseLinacNode):
	"""
	The abstract class of the linac accelerator elements hierarchy that can be tilted.
	"""
	def __init__(self, name = "none"):
		BaseLinacNode.__init__(self,name)
		self.__tiltNodeIN  = TiltElement()
		self.__tiltNodeOUT = TiltElement()
		self.__tiltNodeIN.setName(name+"_tilt_in")
		self.__tiltNodeOUT.setName(name+"_tilt_out")
		self.addChildNode(self.__tiltNodeIN,AccNode.ENTRANCE)
		self.addChildNode(self.__tiltNodeOUT,AccNode.EXIT)
		self.addParam("tilt",self.__tiltNodeIN.getTiltAngle())
		self.setType("linacNode")

	def setTiltAngle(self, angle = 0.):
		"""
		Sets the tilt angle for the tilt operation.
		"""
		self.setParam("tilt", angle)
		self.__tiltNodeIN.setTiltAngle(angle)
		self.__tiltNodeOUT.setTiltAngle( (-1.0) * angle)

	def getTiltAngle(self):
		"""
		Returns the tilt angle for the tilt operation.
		"""
		return self.__tiltNodeIN.getTiltAngle()

	def getNodeTiltIN(self):
		"""
		Returns the Tilt Node instance before this node
		"""
		return self.__tiltNodeIN 
 
	def getNodeTiltOUT(self):
		"""
		Returns the  Tilt Node instance after this node
		"""
		return self.__tiltNodeOUT	


class LinacMagnetNode(LinacNode):
	"""
	The abstract class of the linac magnet.
	"""
	def __init__(self, name = "none"):
		LinacNode.__init__(self,name)
		self.__fringeFieldIN = FringeField(self)
		self.__fringeFieldOUT = FringeField(self)
		self.__fringeFieldIN.setName(name+"_fringe_in")
		self.__fringeFieldOUT	.setName(name+"_fringe_out")	
		self.addChildNode(self.__fringeFieldIN,AccNode.ENTRANCE)
		self.getChildNodes(AccNode.EXIT).insert(0,self.__fringeFieldOUT)
		self.setType("linacMagnet")		
		
	def getField(self):
		"""
		Returns the field value for the magnet.
		Abstract method. It should be implemented in the subclass.
		"""
		return None
		
	def setField(self,field):
		"""
		Sets the field value for the magnet.
		Abstract method. It should be implemented in the subclass.
		"""
		pass
		
	def getNodeFringeFieldIN(self):
		"""
		Returns the FringeField instance before this node
		"""
		return self.__fringeFieldIN 
 
	def getNodeFringeFieldOUT(self):
		"""
		Returns the FringeField instance after this node
		"""
		return self.__fringeFieldOUT 
		
	def setFringeFieldFunctionIN(self, trackFunction = None):
		"""
		Sets the fringe field function that will track the bunch
		through the fringe at the entrance of the node.
		"""
		self.__fringeFieldIN.setFringeFieldFunction(trackFunction)

	def setFringeFieldFunctionOUT(self, trackFunction = None):
		"""
		Sets the fringe field function that will track the bunch
		through the fringe at the exit of the element.
		"""
		self.__fringeFieldOUT.setFringeFieldFunction(trackFunction)

	def getFringeFieldFunctionIN(self, trackFunction = None):
		"""
		Returns the fringe field function that will track the bunch
		through the fringe at the entrance of the node.
		"""
		return self.__fringeFieldIN.getFringeFieldFunction()

	def getFringeFieldFunctionOUT(self, trackFunction = None):
		"""
		Returns the fringe field function that will track the bunch
		through the fringe at the exit of the element.
		"""
		return self.__fringeFieldOUT.getFringeFieldFunction()

	def setUsageFringeFieldIN(self,usage = True):
		"""
		Sets the property describing if the IN fringe
		field will be used in calculation.
		"""
		self.__fringeFieldIN.setUsage(usage)

	def getUsageFringeFieldIN(self):
		"""
		Returns the property describing if the IN fringe
		field will be used in calculation.
		"""
		return self.__fringeFieldIN.getUsage()

	def setUsageFringeFieldOUT(self,usage = True):
		"""
		Sets the property describing if the OUT fringe
		field will be used in calculation.
		"""
		self.__fringeFieldOUT.setUsage(usage)

	def getUsageFringeFieldOUT(self):
		"""
		Returns the property describing if the OUT fringe
		field will be used in calculation.
		"""
		return self.__fringeFieldOUT.getUsage()
		
#-----------------------------------------------------
#   The REAL LINAC NODES 
#-----------------------------------------------------

class Drift(BaseLinacNode):
	"""
	Drift element.
	"""
	def __init__(self, name = "drift"):
		"""
		Constructor. Creates the Drift element.
		"""
		BaseLinacNode.__init__(self,name)
		self.setType("drift")

	def track(self, paramsDict):
		"""
		The drift class implementation of the AccNode class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		self.tracking_module.drift(bunch, length)
		
class Quad(LinacMagnetNode):
	"""
	Quad Combined Function TEAPOT element.
	"""
	def __init__(self, name = "quad"):
		"""
		Constructor. Creates the Quad Combined Function element.
		Example: poleArr = [2,3] - for sextupoles and octupoles
		"""
		LinacMagnetNode.__init__(self,name)	
		self.addParam("dB/dr",0.)
		self.addParam("poles",[])
		self.addParam("kls",[])
		self.addParam("skews",[])
		self.setnParts(2)
		self.setType("linacQuad")

		def fringeIN(node,paramsDict):
			# B*rho = 3.335640952*momentum/charge [T*m] if momentum in GeV/c
			usageIN = node.getUsage()	
			if(not usageIN):
				return
			bunch = paramsDict["bunch"]
			charge = bunch.charge()
			momentum = bunch.getSyncParticle().momentum()
			#---- The charge sign will be accounted for inside tracking module functions.			
			kq = node.getParam("dB/dr")/bunch.B_Rho()
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			length = paramsDict["parentNode"].getLength()
			TPB.quadfringeIN(bunch,kq)
			if(length == 0.):
				return
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				k = klArr[i]/length
				skew = skewArr[i]
				TPB.multpfringeIN(bunch,pole,k,skew)

		def fringeOUT(node,paramsDict):
			# B*rho = 3.335640952*momentum/charge [T*m] if momentum in GeV/c
			usageOUT = node.getUsage()
			if(not usageOUT):
				return
			bunch = paramsDict["bunch"]
			charge = bunch.charge()
			momentum = bunch.getSyncParticle().momentum()
			#---- The charge sign will be accounted for inside tracking module functions
			kq = node.getParam("dB/dr")/bunch.B_Rho()	
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			length = paramsDict["parentNode"].getLength()
			bunch = paramsDict["bunch"]
			TPB.quadfringeOUT(bunch,kq)
			if(length == 0.):
				return
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				k = klArr[i]/length
				skew = skewArr[i]
				TPB.multpfringeOUT(bunch,pole,k,skew)

		self.setFringeFieldFunctionIN(fringeIN)
		self.setFringeFieldFunctionOUT(fringeOUT)
		self.getNodeTiltIN().setType("quad tilt in")
		self.getNodeTiltOUT().setType("quad tilt out")
		self.getNodeFringeFieldIN().setType("quad fringe in")
		self.getNodeFringeFieldOUT().setType("quad fringe out")
		
	def getField(self):
		"""
		Returns the field value for the magnet..
		"""
		return self.getParam("dB/dr")
		
	def setField(self,field):
		"""
		Sets the field value for the magnet.
		"""
		self.setParam("dB/dr",field)	
		
	def initialize(self):
		"""
		The  Quad Combined Function class implementation
		of the AccNode class initialize() method.
		The parts number should be even to have ability to put correctors
		at the center of the quad.
		"""
		nParts = self.getnParts()
		if(nParts < 2 and nParts%2 != 0):
			msg = "The Quad Combined Function class instance should have no less than 2 and even number of parts!"
			msg = msg + os.linesep
			msg = "Even number of parts are needed for the correctors at the center of the quad."
			msg = msg + os.linesep			
			msg = msg + "Method initialize():"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "nParts =" + str(nParts)
			orbitFinalize(msg)
		lengthStep = self.getLength()/nParts
		for i in range(nParts):
			self.setLength(lengthStep,i)
		"""
		#=============================================
		# This is an old TEAPOT-like implementation 
		# of the Quad slicing.
		#=============================================		
		lengthIN = (self.getLength()/(nParts - 1))/2.0
		lengthOUT = (self.getLength()/(nParts - 1))/2.0
		lengthStep = lengthIN + lengthOUT
		self.setLength(lengthIN,0)
		self.setLength(lengthOUT,nParts - 1)
		for i in range(nParts-2):
			self.setLength(lengthStep,i+1)
		"""	
		
	def track(self, paramsDict):
		"""
		The Quad Combined Function TEAPOT  class implementation
		of the AccNode class track(probe) method.
		"""
		bunch = paramsDict["bunch"]
		charge = bunch.charge()
		momentum = bunch.getSyncParticle().momentum()
		#---- The sign of dB/dr will be delivered to tracking module
		#---- functions as kq.
		#---- The charge sign will be accounted for inside tracking module
		#---- functions.
		kq = self.getParam("dB/dr")/bunch.B_Rho()
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
		poleArr = self.getParam("poles")
		klArr = self.getParam("kls")
		skewArr = self.getParam("skews")
		#print "debug name =",self.getName()," kq=",kq,"  L=",self.getLength(index)," index=",index
		#===============================================================
		# This is a 3-sub-parts implementation TEAPOT algorithm
		# Now the quad is divided into equally long parts, and
		# for each part we use the same tracking.
		# It will allow the backward tracking through the 
		# quad with Space Charge child nodes.
		#===============================================================
		step = length
		self.tracking_module.quad1(bunch, step/4, kq)
		self.tracking_module.quad2(bunch, step/4)
		for i in xrange(len(poleArr)):
			pole = poleArr[i]
			kl = klArr[i]/(2*nParts)
			skew = skewArr[i]
			TPB.multp(bunch,pole,kl,skew)
		self.tracking_module.quad2(bunch, step/4)
		self.tracking_module.quad1(bunch, step/2, kq)
		self.tracking_module.quad2(bunch, step/4)
		for i in xrange(len(poleArr)):
			pole = poleArr[i]
			kl = klArr[i]/(2*nParts)
			skew = skewArr[i]
			TPB.multp(bunch,pole,kl,skew)
		self.tracking_module.quad2(bunch, step/4)
		self.tracking_module.quad1(bunch, step/4, kq)
		"""
		#=============================================
		# This is an old TEAPOT-like implementation 
		# of the Quad tracking.
		#=============================================
		if(index == 0):
			self.tracking_module.quad1(bunch, length, kq)
			return
		if(index > 0 and index < (nParts-1)):
			self.tracking_module.quad2(bunch, length/2.0)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			self.tracking_module.quad2(bunch, length/2.0)
			self.tracking_module.quad1(bunch, length, kq)
			return
		if(index == (nParts-1)):
			self.tracking_module.quad2(bunch, length)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]*kq*length/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			self.tracking_module.quad2(bunch, length)
			self.tracking_module.quad1(bunch, length, kq)
		"""		
		return		

	def getTotalField(self,z):
		"""
		Returns the field of the quad.
		This function was added to make a uniform 
		interface with OverlappingQuadsNode.
		"""
		G = 0.
		if(abs(z) < self.getLength()/2.):
			G = self.getParam("dB/dr")
		return G


class Bend(LinacMagnetNode):
	"""
	Bend Combined Functions TEAPOT element.
	"""
	def __init__(self, name = "bend no name"):
		"""
		Constructor. Creates the Bend Combined Functions TEAPOT element .
		"""
		LinacMagnetNode.__init__(self,name)
		self.addParam("poles",[])
		self.addParam("kls",[])
		self.addParam("skews",[])

		self.addParam("ea1",0.)
		self.addParam("ea2",0.)
		self.addParam("rho",0.)
		self.addParam("theta",1.0e-36)
		
		self.setnParts(2)
		self.setType("bend linac")
		self.setUsageFringeFieldIN(True)
		self.setUsageFringeFieldOUT(True)
		
		def fringeIN(node,paramsDict):
			bunch = paramsDict["bunch"]
			length = paramsDict["parentNode"].getLength()
			usageIN = node.getUsage()
			e = node.getParam("ea1")
			rho = node.getParam("rho")
			poleArr = node.getParam("poles")
			klArr =  [-x*bunch.charge()*length for x in self.getParam("kls")]
			skewArr = node.getParam("skews")
			nParts = paramsDict["parentNode"].getnParts()
			if(e != 0.):
				inout = 0
				TPB.wedgedrift(bunch,e,inout)
				if(usageIN):
					frinout = 0
					TPB.wedgerotate(bunch, e, frinout)
					TPB.bendfringeIN(bunch, rho)
					if(length != 0.):
						for i in xrange(len(poleArr)):
							pole = poleArr[i]
							kl = klArr[i]/length
							skew = skewArr[i]
							TPB.multpfringeIN(bunch,pole,kl,skew)
					frinout = 1
					TPB.wedgerotate(bunch, e, frinout)                                
				TPB.wedgebendCF(bunch, e, inout, rho, len(poleArr), poleArr, klArr, skewArr, nParts - 1)
			else:
				if(usageIN):
					TPB.bendfringeIN(bunch, rho)
					if(length != 0.):
						for i in xrange(len(poleArr)):
							pole = poleArr[i]
							kl = klArr[i]/length
							skew = skewArr[i]
							TPB.multpfringeIN(bunch,pole,kl,skew)

		def fringeOUT(node,paramsDict):
			bunch = paramsDict["bunch"]
			length = paramsDict["parentNode"].getLength()
			usageOUT = node.getUsage()
			e = node.getParam("ea2")
			rho = node.getParam("rho")
			poleArr = node.getParam("poles")
			klArr =  [-x*bunch.charge()*length for x in self.getParam("kls")]
			skewArr = node.getParam("skews")
			nParts = paramsDict["parentNode"].getnParts()
			if(e != 0.):
				inout = 1
				TPB.wedgebendCF(bunch, e, inout, rho, len(poleArr), poleArr, klArr, skewArr, nParts - 1)
				if(usageOUT):
					frinout = 0
					TPB.wedgerotate(bunch, -e, frinout)
					TPB.bendfringeOUT(bunch, rho)
					if(length != 0.):
						for i in xrange(len(poleArr)):
							pole = poleArr[i]
							kl = klArr[i]/length
							skew = skewArr[i]
							TPB.multpfringeOUT(bunch,pole,kl,skew)
					frinout = 1
					TPB.wedgerotate(bunch, -e, frinout)
				TPB.wedgedrift(bunch,e,inout)
			else:
				if(usageOUT):
					TPB.bendfringeOUT(bunch, rho)
					if(length != 0.):
						for i in xrange(len(poleArr)):
							pole = poleArr[i]
							kl = klArr[i]/length
							skew = skewArr[i]
							TPB.multpfringeOUT(bunch,pole,kl,skew)

		self.setFringeFieldFunctionIN(fringeIN)
		self.setFringeFieldFunctionOUT(fringeOUT)

	def initialize(self):
		"""
		The  Bend Combined Functions TEAPOT class implementation of
		the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts < 2):
			msg = "The Bend Combined Functions TEAPOT class instance should have more than 2 parts!"
			msg = msg + os.linesep
			msg = msg + "Method initialize():"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "nParts =" + str(nParts)
			orbitFinalize(msg)
		rho = self.getLength()/self.getParam("theta")
		self.addParam("rho",rho)
		lengthIN = (self.getLength()/(nParts - 1))/2.0
		lengthOUT = (self.getLength()/(nParts - 1))/2.0
		lengthStep = lengthIN + lengthOUT
		self.setLength(lengthIN,0)
		self.setLength(lengthOUT,nParts - 1)
		for i in xrange(nParts-2):
			self.setLength(lengthStep,i+1)

	def track(self, paramsDict):
		"""
		The Bend Combined Functions TEAPOT  class implementation of
		the AccNodeBunchTracker class track(probe) method.
		"""
                bunch = paramsDict["bunch"]
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
		poleArr = self.getParam("poles")
		klArr =  [-x*bunch.charge()*self.getLength() for x in self.getParam("kls")]
		skewArr = self.getParam("skews")		
		theta = self.getParam("theta")/(nParts - 1)
		if(index == 0):
			TPB.bend1(bunch, length, theta/2.0)
			return
		if(index > 0 and index < (nParts-1)):
			TPB.bend2(bunch, length/2.0)
			TPB.bend3(bunch, theta/2.0)
			TPB.bend4(bunch,theta/2.0)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			TPB.bend4(bunch,theta/2.0)
			TPB.bend3(bunch, theta/2.0)
			TPB.bend2(bunch, length/2.0)
			TPB.bend1(bunch, length, theta)
			return
		if(index == (nParts-1)):
			TPB.bend2(bunch, length)
			TPB.bend3(bunch, theta/2.0)
			TPB.bend4(bunch, theta/2.0)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			TPB.bend4(bunch, theta/2.0)
			TPB.bend3(bunch, theta/2.0)
			TPB.bend2(bunch, length)
			TPB.bend1(bunch, length, theta/2.0)
		return

class DCorrectorH(LinacMagnetNode):
	"""
	The Horizontal Dipole Corrector.
	"""
	def __init__(self, name = "correctorh"):
		"""
		Constructor. Creates the Horizontal Dipole Corrector element .
		"""
		LinacMagnetNode.__init__(self,name)
		self.addParam("B",0.)
		self.addParam("effLength",0.)
		self.setType("dch")	
		self.setnParts(1)

	def getField(self):
		"""
		Returns the field value for the magnet..
		"""
		return self.getParam("B")
		
	def setField(self,field):
		"""
		Sets the field value for the magnet.
		"""
		self.setParam("B",field)

	def track(self, paramsDict):
		"""
		The Horizontal Dipole Corrector class implementation of
		the AccNode class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getParam("effLength")/nParts
		field = self.getParam("B")
		bunch = paramsDict["bunch"]
		charge = bunch.charge()
		syncPart = bunch.getSyncParticle()
		momentum = syncPart.momentum()
		# dp/p = Q*c*B*L/p p in GeV/c c = 2.99792*10^8/10^9
		kick = -field*charge*length*0.299792/momentum
		self.tracking_module.kick(bunch,kick,0.,0.)

class DCorrectorV(LinacMagnetNode):
	"""
	The Vertical Dipole Corrector.
	"""
	def __init__(self, name = "correctorv"):
		"""
		Constructor. Creates the Vertical Dipole Corrector element .
		"""
		LinacMagnetNode.__init__(self,name)
		self.addParam("B",0.)
		self.addParam("effLength",0.)
		self.setType("dcv")	
		self.setnParts(1)

	def getField(self):
		"""
		Returns the field value for the magnet..
		"""
		return self.getParam("B")
		
	def setField(self,field):
		"""
		Sets the field value for the magnet.
		"""
		self.setParam("B",field)

	def track(self, paramsDict):
		"""
		The Vertical Dipole Corrector class implementation of
		the AccNode class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getParam("effLength")/nParts
		field = self.getParam("B")
		bunch = paramsDict["bunch"]
		charge = bunch.charge()
		syncPart = bunch.getSyncParticle()
		momentum = syncPart.momentum()
		# dp/p = Q*c*B*L/p p in GeV/c, c = 2.99792*10^8/10^9
		kick = field*charge*length*0.299792/momentum
		self.tracking_module.kick(bunch,0,kick,0.)

class ThickKick(LinacMagnetNode):
	"""
	Thick kiker linac node.
	"""
	def __init__(self, name = "thick_kick"):
		"""
		Constructor. Creates the thick kiker linac node.
		"""
		LinacMagnetNode.__init__(self,name)	
		self.addParam("Bx",0.)
		self.addParam("By",[])
		self.setnParts(1)
		self.setType("thickKick")
		
	def getFieldBx(self):
		"""
		Returns the field Bx value for the magnet. Vertical kick.
		"""
		return self.getParam("Bx")
		
	def setFieldBx(self,field):
		"""
		Sets the field Bx value for the magnet. Vertical kick.
		"""
		self.setParam("Bx",field)	
		
	def getFieldBy(self):
		"""
		Returns the field By value for the magnet. Horizontal kick.
		"""
		return self.getParam("By")
		
	def setFieldBy(self,field):
		"""
		Sets the field By value for the magnet. Horizontal kick.
		"""
		self.setParam("By",field)		
		
	def initialize(self):
		"""
		The  Thick Kick class implementation of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		lengthStep = self.getLength()/nParts
		for i in range(nParts):
			self.setLength(lengthStep,i)
		
	def track(self, paramsDict):
		"""
		The Thick Kick  class implementation of the AccNode class track(probe) method.
		"""
		bunch = paramsDict["bunch"]
		charge = bunch.charge()
		momentum = bunch.getSyncParticle().momentum()
		Bx = self.getParam("Bx")
		By = self.getParam("By")	
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
		#print "debug name =",self.getName()," Bx=",Bx," By=",By,"  L=",self.getLength(index)," index=",index
		#==========================================
		# dp/p = Q*c*B*L/p p in GeV/c, c = 2.99792*10^8/10^9
		kickY = +Bx*charge*length*0.299792/momentum
		kickX = -By*charge*length*0.299792/momentum
		self.tracking_module.drift(bunch, length/2.0)
		self.tracking_module.kick(bunch,kickX,kickY,0.)
		self.tracking_module.drift(bunch, length/2.0)

class AbstractRF_Gap(BaseLinacNode):
	"""
	This is an abstarct class for all RF Gap classes.
	"""
	def __init__(self, name = "abstractfgap"):
		"""
		Constructor for the abstract RF gap.
		"""
		BaseLinacNode.__init__(self,name)	
		self.addParam("gap_phase",0.)
		self.addParam("rfCavity", None)
		self.setLength(0.)
		self.setType("abstractrfgap")	
		self.__isFirstGap = False
	
	def initialize(self):
		"""
		The AbstractRF_Gap class initialize() method.
		"""
		BaseLinacNode.initialize(self)
	
	def isRFGap(self):
		"""
		Returns True.
		"""
		return True

	def isFirstRFGap(self):
		"""
		Returns True if it is the first gap in RF cavity. 
		"""
		return self.__isFirstGap

	def setAsFirstRFGap(self, isFirst):
		"""
		Sets if it is the first gap in RF cavity. 
		"""
		self.__isFirstGap = isFirst
	
	def setRF_Cavity(self, rf_cav):
		"""
		Sets the parent RF Cavity.
		"""
		self.addParam("rfCavity",rf_cav)

	def getRF_Cavity(self):
		"""
		Returns the parent RF Cavity.
		"""
		return self.getParam("rfCavity")
	
	def setGapPhase(self, gap_phase):
		"""
		Sets the rf gap phase.
		"""
		self.setParam("gap_phase",gap_phase)

	def getGapPhase(self):
		"""
		Returns the rf gap phase.
		"""
		return self.getParam("gap_phase")
		
	def track(self, paramsDict):
		"""
		The subclasses AbstractRF_Gap class should implement
		the BaseLinacNode class track(probe) method.
		"""
		pass
	
	def trackDesign(self, paramsDict):
		"""
		The subclasses AbstractRF_Gap class should implement
		the BaseLinacNode class trackDesign(probe) method.
		"""
		pass

#------------------------------------------------------------
#     Auxilary classes 
#------------------------------------------------------------	
		
class TiltElement(BaseLinacNode):
	"""
	The class to do tilt at the entrance of an element.
	"""
	def __init__(self, name = "tilt no name", angle = 0.):
		"""
		Constructor. Creates the Tilt element.
		"""
		AccNode.__init__(self,name)
		self.__angle = angle
		self.setType("tilt teapot")

	def setTiltAngle(self, angle = 0.):
		"""
		Sets the tilt angle for the tilt operation.
		"""
		self.__angle = angle

	def getTiltAngle(self):
		"""
		Returns the tilt angle for the tilt operation.
		"""
		return self.__angle

	def track(self, paramsDict):
		"""
		It is tracking the dictionary with parameters through
		the titlt node.
		"""
		if(self.__angle != 0.):
			bunch = paramsDict["bunch"]
			TPB.rotatexy(bunch,self.__angle)


class FringeField(BaseLinacNode):
	"""
	The class is a base class for the fringe field classes for others elements.
	"""
	def __init__(self,  parentNode,  trackFunction = None , name = "fringe field no name"):
		"""
		Constructor. Creates the Fringe Field element.
		"""
		AccNode.__init__(self,name)
		self.setParamsDict(parentNode.getParamsDict())
		self.__trackFunc = trackFunction
		self.__usage = False
		self.setType("fringeField teapot")

	def track(self, paramsDict):
		"""
		It is tracking the dictionary with parameters through
		the fringe field node.
		"""
		if(self.__trackFunc != None and self.__usage == True):
			self.__trackFunc(self,paramsDict)

	def setFringeFieldFunction(self, trackFunction = None):
		"""
		Sets the fringe field function that will track the bunch through the fringe.
		"""
		self.__trackFunc = trackFunction

	def getFringeFieldFunction(self):
		"""
		Returns the fringe field function.
		"""
		return self.__trackFunc

	def setUsage(self,usage = True):
		"""
		Sets the boolean flag describing if the fringe
		field will be used in calculation.
		"""
		self.__usage = usage

	def getUsage(self):
		"""
		Returns the boolean flag describing if the fringe
		field will be used in calculation.
		"""
		return self.__usage


