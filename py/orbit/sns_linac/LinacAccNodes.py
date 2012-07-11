"""
This module is a collection of the linac accelerator nodes which are the subclasses of 
the AccNode class. We cannot use here TEAPOT nodes from the TEAPOT package 
directly because we use the field as a parameter for quads and dipole correctors
instead of k1 = 1/(B*rho)*(dB/dr). The RF Cavities are different from the ring RF.
"""

import os
import math

# import the function that creates multidimensional arrays
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

# from linac import the RF gap classes
from linac import BaseRfGap, MatrixRfGap

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
		self.__linacSeqence = None
		
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
		self.__params["tilt"] = angle
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
		

class TiltNode(BaseLinacNode):
	"""
	The class to do tilt at the entrance of an node.
	"""
	def __init__(self, name = "tilt", angle = 0.):
		"""
		Constructor. Creates the Tilt Node.
		"""
		AccNode.__init__(self,name)
		self.__angle = angle
		self.setType("tilt")

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

#-----------------------------------------------------
#    LINAC ABST NODES ELEMENTS
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
		TPB.drift(bunch, length)
		
class Quad(LinacMagnetNode):
	"""
	Quad Combined Function TEAPOT element.
	"""
	def __init__(self, name = "quad"):
		"""
		Constructor. Creates the Quad Combined Function element .
		"""
		LinacMagnetNode.__init__(self,name)	
		self.addParam("dB/dr",0.)
		self.addParam("poles",[])
		self.addParam("kls",[])
		self.addParam("skews",[])
		self.setnParts(2)
		self.setType("linacQuad")
		
		# B*rho = 3.335640952*momentum [T*m] if momentum in GeV/c
		def fringeIN(node,paramsDict):
			usageIN = node.getUsage()	
			if(not usageIN):
				return
			bunch = paramsDict["bunch"]	
			momentum = bunch.getSyncParticle().momentum()
			kq = node.getParam("dB/dr")/(3.335640952*momentum)
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			length = paramsDict["parentNode"].getLength()
			TPB.quadfringeIN(bunch,kq)
			if(length == 0.):
				return
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				k = klArr[i]*kq
				skew = skewArr[i]
				TPB.multpfringeIN(bunch,pole,k,skew)

		def fringeOUT(node,paramsDict):
			usageOUT = node.getUsage()
			if(not usageOUT):
				return
			bunch = paramsDict["bunch"]
			momentum = bunch.getSyncParticle().momentum()
			kq = node.getParam("dB/dr")/(3.335640952*momentum)	
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
				k = klArr[i]*kq
				skew = skewArr[i]
				TPB.multpfringeOUT(bunch,pole,k,skew)

		self.setFringeFieldFunctionIN(fringeIN)
		self.setFringeFieldFunctionOUT(fringeOUT)
		self.getNodeTiltIN().setType("quad tilt in")
		self.getNodeTiltOUT().setType("quad tilt out")
		self.getNodeFringeFieldIN().setType("quad fringe in")
		self.getNodeFringeFieldOUT().setType("quad fringe out")
		

	def initialize(self):
		"""
		The  Quad Combined Function class implementation
		of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts < 2 and nParts%2 != 0):
			msg = "The Quad Combined Function class instance should have no less than 2 and even number of parts!"
			msg = msg + os.linesep
			msg = msg + "Method initialize():"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "nParts =" + str(nParts)
			orbitFinalize(msg)
		lengthIN = (self.getLength()/(nParts - 1))/2.0
		lengthOUT = (self.getLength()/(nParts - 1))/2.0
		lengthStep = lengthIN + lengthOUT
		self.setLength(lengthIN,0)
		self.setLength(lengthOUT,nParts - 1)
		for i in xrange(nParts-2):
			self.setLength(lengthStep,i+1)


	def track(self, paramsDict):
		"""
		The Quad Combined Function TEAPOT  class implementation
		of the AccNode class track(probe) method.
		"""
		bunch = paramsDict["bunch"]		
		momentum = bunch.getSyncParticle().momentum()
		kq = self.getParam("dB/dr")/(3.335640952*momentum)	
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
		poleArr = self.getParam("poles")
		klArr = self.getParam("kls")
		skewArr = self.getParam("skews")
		#print "debug name =",self.getName()," kq=",kq,"  L=",self.getLength()		
		if(index == 0):
			TPB.quad1(bunch, length, kq)
			return
		if(index > 0 and index < (nParts-1)):
			TPB.quad2(bunch, length/2.0)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			TPB.quad2(bunch, length/2.0)
			TPB.quad1(bunch, length, kq)
			return
		if(index == (nParts-1)):
			#print "debug before xp",	bunch.xp(0)
			TPB.quad2(bunch, length)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]*kq*length/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			TPB.quad2(bunch, length)
			TPB.quad1(bunch, length, kq)
			#print "debug after xp",	bunch.xp(0)
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
		syncPart = bunch.getSyncParticle()
		momentum = syncPart.momentum()
		q = bunch.charge()
		# dp/p = Q*c*B*L/p p in GeV/c c = 2.99792*10^8/10^9		
		kick = q*field*length*0.299792/momentum
		TPB.kick(bunch,kick,0.,0.)

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
		syncPart = bunch.getSyncParticle()
		momentum = syncPart.momentum()
		q = bunch.charge()
		# dp/p = Q*c*B*L/p p in GeV/c, c = 2.99792*10^8/10^9		
		kick = q*field*length*0.299792/momentum
		TPB.kick(bunch,0,kick,0.)

class BaseRF_Gap(BaseLinacNode):
	"""
	The simplest RF gap representation. The only E*T*L defines all effects of the node.
	"""
	def __init__(self, name = "baserfgap"):
		"""
		Constructor for the simplest RF gap. E0TL parameter is in GeV. Phases are in radians.
		It has 3 parts with lengthes: 0.5 + 0. + 0.5 
		"""
		BaseLinacNode.__init__(self,name)
		self.addParam("E0TL",0.)
		self.addParam("modePhase",0.)		
		self.addParam("rfCavity", None)
		self.setType("baserfgap")	
		self.__isFirstGap = False
		self.setnParts(3)	
		self.cppGapModel = MatrixRfGap()

	def setCppGapModel(self, cppGapModel = MatrixRfGap):
		"""
		This method will set the fast c++ simple model for the RF Gap. 
		By default it is Matrix RF Gap model which is a linear transport matrix.
		"""
		self.cppGapModel = cppGapModel()
	
	def initialize(self):
		"""
		The Ring RF TEAPOT class implementation
		of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts != 3):
			msg = "The simple Rf gap should have 3 parts!"
			msg = msg + os.linesep
			msg = msg + "Method initialize():"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "nParts =" + str(nParts)
			msg = msg + os.linesep
			msg = msg + "lenght =" + str(self.getLength())
			orbitFinalize(msg)
		length = self.getLength()
		self.setLength(length/2.0,0)
		self.setLength(0.,1)
		self.setLength(length/2.0,2)
	
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
	
	def track(self, paramsDict):
		"""
		The simplest RF gap class implementation of
		the AccNode class track(probe) method.
		"""
		index = self.getActivePartIndex()
		length = self.getLength(index)
		bunch = paramsDict["bunch"]		
		syncPart = bunch.getSyncParticle()
		if(index == 0 or index == 2):
			gapOffset = 0.
			if(self.hasParam("gapOffset")): gapOffset = self.getParam("gapOffset")
			if(index == 2): gapOffset = -gapOffset
			TPB.drift(bunch, length + gapOffset)
			return
		E0TL = self.getParam("E0TL")		
		modePhase = self.getParam("modePhase")*math.pi
		rfCavity = self.getRF_Cavity()
		frequency = rfCavity.getFrequency()	
		rfPhase = rfCavity.getPhase() + modePhase
		phase = rfPhase
		arrival_time = syncPart.time()
		designArrivalTime = rfCavity.getDesignArrivalTime()
		if(self.__isFirstGap):
			if(rfCavity.isDesignSetUp()):
				phase = math.fmod(frequency*(arrival_time - designArrivalTime)*2.0*math.pi + rfPhase,2.0*math.pi)
			else:
				sequence = self.getSequence()
				accLattice = sequence.getLinacAccLattice()
				msg = "The BaseRF_Gap class. You have to run trackDesign on the LinacAccLattice first to initialize all RF Cavities' phases!"
				msg = msg + os.linesep
				msg = msg + "Lattice =" + accLattice.getName()				
				msg = msg + os.linesep
				msg = msg + "Sequence =" + sequence.getName()				
				msg = msg + os.linesep
				msg = msg + "RF Cavity =" + rfCavity.getName()				
				msg = msg + os.linesep
				msg = msg + "Name of element=" + self.getName()
				msg = msg + os.linesep
				msg = msg + "Type of element=" + self.getType()
				msg = msg + os.linesep
				orbitFinalize(msg)				
		else:
			phase = math.fmod(frequency*(arrival_time - designArrivalTime)*2.0*math.pi+rfPhase,2.0*math.pi)	
		#------------------------------------------------------
		#call rf gap with E0TL phase phase of the gap and a longitudinal shift parameter	
		self.cppGapModel.trackBunch(bunch,frequency,E0TL,phase)
		#print "debug RF E0TL=",E0TL," phase=",phase*180./math.pi," eKin[MeV]=",bunch.getSyncParticle().kinEnergy()*1.0e+3		
		
	def trackDesign(self, paramsDict):
		"""
		The RF First Gap node setups the design time of passage 
		of the bunch through this node.
		"""
		index = self.getActivePartIndex()
		length = self.getLength(index)
		bunch = paramsDict["bunch"]
		if(index == 0 or index == 2):
			gapOffset = 0.
			if(self.hasParam("gapOffset")): gapOffset = self.getParam("gapOffset")
			if(index == 2): gapOffset = -gapOffset
			TPB.drift(bunch, length + gapOffset)
			return		
		E0TL = self.getParam("E0TL")			
		rfCavity = self.getRF_Cavity()
		modePhase = self.getParam("modePhase")*math.pi	
		arrival_time = bunch.getSyncParticle().time()
		frequency = rfCavity.getFrequency()	
		rfPhase = rfCavity.getDesignPhase() + modePhase
		phase = rfPhase
		if(self.__isFirstGap):
			rfCavity.setDesignArrivalTime(arrival_time)
			rfCavity.setDesignSetUp(True)			
		else:
			first_gap_arr_time = rfCavity.getDesignArrivalTime()
			#print "debug name=",self.getName()," delta_phase=",frequency*(arrival_time - first_gap_arr_time)*360.0," rfPhase=",rfPhase*180/math.pi
			phase = math.fmod(frequency*(arrival_time - first_gap_arr_time)*2.0*math.pi+rfPhase,2.0*math.pi)		
		#print "debug name=",self.getName()," arr_time=",arrival_time," phase=",phase*180./math.pi," E0TL=",E0TL*1.0e+3," freq=",frequency
		#------------------------------------------------------
		#call rf gap with E0TL phase phase of the gap and a longitudinal shift parameter	
		self.cppGapModel.trackBunch(bunch,frequency,E0TL,phase)
		#syncPart = bunch.getSyncParticle()
		#eKin = syncPart.kinEnergy()
		#print "debug RF E0TL=",E0TL," phase=",phase*180./math.pi," eKin[MeV]=",bunch.getSyncParticle().kinEnergy()*1.0e+3


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
		if(self.__trackFunc != None):
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

#----------------------------------------------------------------
# Classes that are specific for the linac model
#----------------------------------------------------------------

class RF_Cavity(NamedObject,ParamsDictObject):
	"""
	This is the class to keep refernces to the RF Gaps which are BaseLinacNode
	subclasses. This class does not belong to the AccNodes.
	"""
	def __init__(self, name = "none"):
		NamedObject.__init__(self, name)
		ParamsDictObject.__init__(self)
		self.__rfGaps = []
		self.addParam("frequency",0.)
		self.addParam("phase",0.)
		self.addParam("amp",0.)		
		self.addParam("designPhase",0.)
		self.addParam("designAmp",0.)		
		self.addParam("designArrivalTime",0.)
		self.addParam("isDesignSetUp",False)
		
	def setDesignSetUp(self,designOnOf):
		""" Sets the design set up information (yes,no). """
		self.setParam("isDesignSetUp",designOnOf)	

	def isDesignSetUp(self):
		""" Returns the design set up information (yes,no). """
		return self.getParam("isDesignSetUp")	
		
	def setDesignArrivalTime(self,time):
		""" Sets the design arrival time for the first RF gap. """
		self.setParam("designArrivalTime",time)
		
	def getDesignArrivalTime(self):
		""" Returns the design arrival time for the first RF gap. """
		return self.getParam("designArrivalTime")
		
	def setDesignPhase(self,phase):
		""" Sets the design phase for the first RF gap. """
		self.setParam("designPhase",phase)
		
	def getDesignPhase(self):
		""" Returns the design phase for the first RF gap. """
		return self.getParam("designPhase")

	def setDesignAmp(self,Amp):
		""" Sets the design Amp for the RF cavity. """
		self.setParam("designAmp",Amp)
		
	def getDesignAmp(self):
		""" Returns the design Amp for the RF cavity. """
		return self.getParam("designAmp")

	def setPhase(self,phase):
		""" Sets the phase for the first RF gap. """
		self.setParam("phase",phase)
		
	def getPhase(self):
		""" Returns the phase for the first RF gap. """
		return self.getParam("phase")

	def setAmp(self,Amp):
		""" Sets the Amp for RF cavity. """
		self.setParam("Amp",Amp)
		
	def getAmp(self):
		""" Returns the Amp for RF cavity. """
		return self.getParam("Amp")
		
	def setFrequency(self,freq):
		""" Sets the frequency in Hz. """
		self.setParam("frequency",freq)
		
	def getFrequency(self):
		""" Returns the frequency in Hz. """
		return self.getParam("frequency")
		
	def addRF_GapNode(self,rfGap):
		""" Adds the rf gap to the cavity."""
		self.__rfGaps.append(rfGap)
		rfGap.setRF_Cavity(self)
		if(len(self.__rfGaps) == 1):
			rfGap.setAsFirstRFGap(True)
		else:
			rfGap.setAsFirstRFGap(False)
		
	def getRF_GapNodes(self):
		""" Returns the array with rf gaps. """
		return self.__rfGaps[:]
	
class Sequence(NamedObject,ParamsDictObject):
	"""
	This is the class to keep refernces to AccNodes that constitute the accelerator sequence.
	"""
	def __init__(self, name = "none"):
		NamedObject.__init__(self, name)
		ParamsDictObject.__init__(self)
		self.__linacNodes = []
		self.addParam("position",0.)	
		self.addParam("length",0.)
		self.addParam("linacAccLattice",None)
		
	def setLinacAccLattice(self,	lattice):
		self.addParam("linacAccLattice",lattice)	
		
	def getLinacAccLattice(self):
		return self.getParam("linacAccLattice")
		
	def addNode(self,node, index = -1):
		""" Adds the Linac Node to the sequence. """
		node.setSequence(self)		
		if(index < 0):
			self.__linacNodes.append(node)
		else:
			self.__linacNodes.insert(index,node)
		
	def getNodes(self):
		""" Returns the array with Linac Nodes. """
		return self.__linacNodes
		
	def setPosition(self, pos):
		""" Sets the position of the sequence. """
		return self.setParam("position",pos)		
		
	def getPosition(self):
		""" Returns the position of the sequence. """
		return self.getParam("position")

	def getLength(self):
		"""
		Returns the total length of the sequence [m].
		"""		
		return self.getParam("length")
		
	def setLength(self, length):
		"""
		Sets the total length of the sequence [m].
		"""		
		return self.setParam("length",length)	
