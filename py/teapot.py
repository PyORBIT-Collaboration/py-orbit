"""
The module includes classes for all TEAPOT elements. The approach is based
on the original ORBIT approach developed by J. Holmes.
"""

import sys
import os
import math

#import teapot base functions from wrapper around c++ functions
import teapot_base as TPB

#import the function that creates multidimensional arrays
from pyORBIT_utils import orbitFinalize

#import the general accelerator elements and lattice
from lattice import AccLattice, AccElement, AccLine

#import the MAD parser. It will be used only for TEAPOT class
from mad_parser import MADparser, LattElement

"""
Drift
Multipole
Quad
Bend
Solenoid
ringRF
Kicker
"""

class TEAPOT:
	"""
	The shell class for the TEAPOT tracker class collection.
	It calls a MAD parser and the TEAPOT accelerator element factory.
	"""
	def __init__(self):
		pass

	def getLattice(self, mad_file_name, lineName):
		"""
		Method. It returns the teapot lattice as an
		instance of AccLattice class from lattice.py package.
		"""
		lattice = AccLattice()
		parser = MADparser()
		parser.parse(mad_file_name)
		accLines = parser.getMADLines()
		if(not accLines.has_key(lineName)):
			print "==============================="
			print "MAD file:",mad_file_name
			print "Can not find accelerator line:",lineName
			print "STOP."
			sys.exit(1)
		#accelerator lines and elements from mad_parser package
		accMAD_Line = accLines[lineName]
		lattice.setName(lineName)
		accMADElements = accMAD_Line.getElements()
		#make TEAPOT lattice elements by using TEAPOT
		# element factory
		for madElm in accMADElements:
			elm = _teapotFactory.getElement(madElm)
			lattice.addChildNode(elm)
		lattice.initialize()
		return lattice

class _teapotFactory:
	"""
	The factory class to produce TEAPOT accelerator elements.
	"""
	def __init__(self):
		pass

	def getElement(self, madElm):
		"""
		Method. It produces the TEAPOT accelerator elements.
		"""
		#madElm = LattElement(" "," ")
		params_ini = madElm.getParameters()
		params = {}
		for par_key in params_ini.iterkeys():
			params[par_key.lower()] = params_ini[par_key]
		#Length parameter
		length = 0.
		if(params.has_key("l")):
			length = params["l"]
		#TILT parameter - if it is there tilt = True,
		#and if it has value tiltAngle = value
		tilt = False
		tiltAngle = None
		if(params.has_key("tilt")):
			tilt = True
			if(params["tilt"] != None):
				tiltAngle = params["tilt"]
		#============DRIFT Element ==================================
		elm = None
		if(madElm.getType().lower() == "drift"):
			elm = DriftTEAPOT(madElm.getName())
		#============Dipole Element SBEND or RBEND ===================
		if(madElm.getType().lower()  == "sbend" or \
			 madElm.getType().lower()  == "rbend"):
			elm = BendTEAPOT(madElm.getName())

			theta = 0.
			if(params.has_key("angle")):
				theta = params["angle"]

			ea1 = 0.
			if(params.has_key("e1")):
				ea1 = params["e1"]

			ea2 = 0.
			if(params.has_key("e2")):
				ea2 = params["e2"]

			if(madElm.getType().lower() == "rbend" ):
				length = length*(theta/2.0)/math.sin(theta/2.0)
				ea1 = ea1 + theta/2.0
				ea2 = ea2 + theta/2.0

			elm.addParam("theta",theta)
			elm.addParam("ea1",ea1)
			elm.addParam("ea2",ea2)
			if(tilt):
				if(tiltAngle == None):
					tiltAngle = math.math.pi/2.0
					if(theta < 0.):
						tiltAngle = - tiltAngle

			if(params.has_key("k1")):
				k1 = params["k1"]
				params["k1l"] = k1*length

			if(params.has_key("k2")):
				k2 = params["k2"]
				params["k2l"] = k2*length

			if(params.has_key("k3")):
				k3 = params["k3"]
				params["k3l"] = k3*length
		#===========QUAD quadrupole element =====================
		if(madElm.getType().lower()  == "quad" or \
			 madElm.getType().lower()  == "quadrupole"):
			elm = QuadTEAPOT(madElm.getName())
			kq = 0.
			if(params.has_key("k1")):
				kq = params["k1"]
			elm.addParam("kq",kq)
			if(tilt):
				if(tiltAngle == None):
					tiltAngle = math.math.pi/4.0
		#===========Sextupole element =====================
		if(madElm.getType().lower()  == "sextupole"):
			elm = MultipoleTEAPOT(madElm.getName())
			k2 = 0.
			if(params.has_key("k2")):
				k2 = params["k2"]
			params["k2l"] = length*k2
			if(tilt):
				if(tiltAngle == None):
					tiltAngle = math.math.pi/6.0
		#===========Octupole element =====================
		if(madElm.getType().lower()  == "octupole"):
			elm = MultipoleTEAPOT(madElm.getName())
			k3 = 0.
			if(params.has_key("k3")):
				k3 = params["k3"]
			params["k2l"] = length*k3
			if(tilt):
				if(tiltAngle == None):
					tiltAngle = math.math.pi/8.0
		#===========Multipole element =====================
		if(madElm.getType().lower()  == "multipole"):
			elm = MultipoleTEAPOT(madElm.getName())
		#===========Solenoid element ======================
		if(madElm.getType().lower()  == "solenoid"):
			elm = SolenoidTEAPOT(madElm.getName())
			ks = 0.
			if(params.has_key("ks")):
				ks = params["ks"]
			elm.addParam("B",ks)
		#===========Kicker element ======================
		if(madElm.getType().lower()  == "kicker"):
			elm = KickTEAPOT(madElm.getName())
			hkick = 0.
			if(params.has_key("hkick")):
				hkick = params["hkick"]
			vkick = 0.
			if(params.has_key("vkick")):
				hkick = params["vkick"]
			elm.addParam("kx",hkick)
			elm.addParam("ky",vkick)
		#===========HKicker element ======================
		if(madElm.getType().lower()  == "hkicker" or \
			 madElm.getType() .lower() == "hkick"):
			elm = KickTEAPOT(madElm.getName())
			hkick = 0.
			if(params.has_key("hkick")):
				hkick = params["hkick"]
			elm.addParam("kx",hkick)
		#===========VKicker element ======================
		if(madElm.getType().lower()  == "vkicker" or \
			 madElm.getType().lower()  == "vkick"):
			elm = KickTEAPOT(madElm.getName())
			vkick = 0.
			if(params.has_key("vkick")):
				hkick = params["vkick"]
			elm.addParam("ky",vkick)
		#===========RF Cavity element ======================
		if(madElm.getType().lower() == "rfcavity"):
			elm = RingRFTEAPOT(madElm.getName())
			#the MAD RF element has a L parameter,
			#but it does not contribute in the ring length!
			"""
			if(length != 0.):
				drft_1 = DriftTEAPOT(madElm.getName()+"_drift")
				drft_2 = DriftTEAPOT(madElm.getName()+"_drift")
				drft_1.setLength(length/2.0)
				drft_2.setLength(length/2.0)
				elm.addChildNodeAtStart(drft_1)
				elm.addChildNodeAtFinish(drft_2)
			"""
			volt = 0.
			if(params.has_key("volt")):
				volt = params["volt"]
			harmon = 0.
			if(params.has_key("harmon")):
				harmon = params["harmon"]
			phase_s = 0.
			if(params.has_key("lag")):
				phase_s = (-1.0) * params["lag"]*2.0*math.pi
			elm.addRF(harmon,volt,phase_s)
		#==========Others elements such as markers,monitor,rcollimator
		if(madElm.getType().lower() == "marker" or \
			 madElm.getType().lower() == "monitor" or \
			 madElm.getType().lower() == "rcolimator"):
			elm = NodeTEAPOT(madElm.getName())
		#------------------------------------------------
		#ready to finish
		#------------------------------------------------
		if(elm == None):
			print "======== Can not create elementwith type:",madElm.getType()
			print "You have to fix the _teapotFactory class."
			print "Stop."
			sys.exit(1)
		#set length
		elm.setLength(length)
		#set tilt angle
		if(tilt == True and tiltAngle != None):
			elm.setTiltAngle(tiltAngle)
		#set K1L,K2L,K3L,... and  T1L,T2L,T3L,...
		poles = []
		kls = []
		skews = []
		for i in xrange(20):
			pole = i
			kl_param = None
			skew = 0
			if(params.has_key("k"+str(pole)+"l")):
				kl_param = params["k"+str(pole)+"l"]
			if(params.has_key("t"+str(pole)+"l")):
				skew = 1
			if(kl_param != None):
				poles.append(pole)
				kls.append(kl_param)
				skews.append(skew)
		if(len(poles) > 0):
			elm.addParams({"poles":poles,"kls":kls,"skews":skews})
		return elm

	getElement = classmethod(getElement)


class BaseTEAPOT(AccElement):
	""" The base class of the TEAPOT accelerator elements hierarchy. """
	def __init__(self, name = "no name"):
		"""
		Constructor of the base TEAPOT element. This is a superclass for all TEAPOT elements.
		"""
		AccElement.__init__(self,name)
		self.setType("base teapot")

	def _initialize(self, paramsDict):
		"""
		The drift class implementation of the AccElement class _initialize(probe) method.
		"""
		nParts = self.getnParts()
		lengthStep = self.getLength()/nParts
		for i in xrange(nParts):
			self.setPartLength(i,lengthStep)

class NodeTEAPOT(BaseTEAPOT):
	def __init__(self, name = "no name"):
		"""
		Constructor of the real TEAPOT element. This is a superclass for all real TEAPOT elements.
		The term real means that element is included in TEAPOT list of elements.
		For instance TILT is not included.
		"""
		BaseTEAPOT.__init__(self,name)
		self.__tiltNodeIN  = TiltTEAPOT()
		self.__tiltNodeOUT = TiltTEAPOT()
		self.__fringeFieldIN = FringeFieldTEAPOT(self)
		self.__fringeFieldOUT = FringeFieldTEAPOT(self)
		self.addChildNodeAtStart(self.__tiltNodeIN)
		self.addChildNodeAtStart(self.__fringeFieldIN)
		self.addChildNodeAtFinish(self.__fringeFieldOUT)
		self.addChildNodeAtFinish(self.__tiltNodeOUT)
		self.addParam("tilt",self.__tiltNodeIN.getTiltAngle())
		self.setType("node teapot")

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


class DriftTEAPOT(NodeTEAPOT):
	"""
	Drift TEAPOT element.
	"""
	def __init__(self, name = "drift no name"):
		"""
		Constructor of the Drift TEAPOT element.
		"""
		NodeTEAPOT.__init__(self,name)
		self.setType("drift teapot")

	def track(self, paramsDict):
		"""
		The drift class implementation of the AccElement class track(probe) method.
		"""
		length = self.getPartLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		TPB.drift(bunch, length)

class SolenoidTEAPOT(NodeTEAPOT):
	"""
	Solenoid TEAPOT element.
	"""
	def __init__(self, name = "solenoid no name"):
		"""
		Constructor of the Solenoid TEAPOT element .
		"""
		NodeTEAPOT.__init__(self,name)
		self.setnParts(1)

		def fringeIN(node,paramsDict):
			B = node.getParam("B")
			bunch = paramsDict["bunch"]
			TPB.solnfringeIN(bunch,B)

		def fringeOUT(node,paramsDict):
			B = node.getParam("B")
			bunch = paramsDict["bunch"]
			TPB.solnfringeOUT(bunch,B)

		self.setFringeFieldFunctionIN(fringeIN)
		self.setFringeFieldFunctionOUT(fringeOUT)

		self.addParam("B",0.)
		self.setType("solenoid teapot")

	def track(self, paramsDict):
		"""
		The Solenoid TEAPOT  class implementation of the AccElement class track(probe) method.
		"""
		index = self.getActivePartIndex()
		length = self.getPartLength(index)
		bunch = paramsDict["bunch"]
		B = node.getParam("B")
		TPB.soln(bunch, length, B)

class MultipoleTEAPOT(NodeTEAPOT):
	"""
	Multipole Combined Function TEAPOT element.
	"""
	def __init__(self, name = "multipole no name"):
		"""
		Constructor of the Multipole Combined Function TEAPOT element.
		"""
		NodeTEAPOT.__init__(self,name)
		self.setnParts(2)
		self.addParam("poles",[])
		self.addParam("kls",[])
		self.addParam("skews",[])

		def fringeIN(node,paramsDict):
			usageIN = node.getUsage()
			if(not usageIN):
				return
			length = paramsDict["parentNode"].getLength()
			if(length == 0.):
				return
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			bunch = paramsDict["bunch"]
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]
				skew = skewArr[i]
				TPB.multpfringeIN(bunch,pole,kl/length,skew)

		def fringeOUT(node,paramsDict):
			usageOUT = node.getUsage()
			if(not usageOUT):
				return
			length = paramsDict["parentNode"].getLength()
			if(length == 0.):
				return
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			bunch = paramsDict["bunch"]
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]
				skew = skewArr[i]
				TPB.multpfringeOUT(bunch,pole,kl/length,skew)

		self.setFringeFieldFunctionIN(fringeIN)
		self.setFringeFieldFunctionOUT(fringeOUT)

		self.setType("multipole teapot")

	def _initialize(self, paramsDict):
		"""
		The  Multipole Combined Function TEAPOT class
		implementation of the AccElement class _initialize(probe) method.
		"""
		nParts = self.getnParts()
		if(nParts < 2):
			msg = "The Multipole TEAPOT class instance should have more than 2 parts!"
			msg = msg + os.linesep
			msg = msg + "Method _initialize(self, paramsDict):"
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
		self.setPartLength(0,lengthIN)
		self.setPartLength(nParts - 1,lengthOUT)
		for i in xrange(nParts-2):
			self.setPartLength(i+1,lengthStep)


	def track(self, paramsDict):
		"""
		The Multipole Combined Function TEAPOT  class
		implementation of the AccElement class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getPartLength(index)
		bunch = paramsDict["bunch"]
		poleArr = self.getParam("poles")
		klArr = self.getParam("kls")
		skewArr = self.getParam("skews")
		if(index == 0):
			TPB.drift(bunch, length)
			return
		if(index > 0 and index < (nParts-1)):
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			TPB.drift(bunch, length)
			return
		if(index == (nParts-1)):
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			TPB.drift(bunch, length)
		return

class QuadTEAPOT(NodeTEAPOT):
	"""
	Quad Combined Function TEAPOT element.
	"""
	def __init__(self, name = "quad no name"):
		"""
		Constructor of the Quad Combined Function TEAPOT element .
		"""
		NodeTEAPOT.__init__(self,name)
		self.setnParts(2)
		self.addParam("kq",0.)
		self.addParam("poles",[])
		self.addParam("kls",[])
		self.addParam("skews",[])

		def fringeIN(node,paramsDict):
			usageIN = node.getUsage()
			if(not usageIN):
				return
			kq = node.getParam("kq")
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			length = paramsDict["parentNode"].getLength()
			bunch = paramsDict["bunch"]
			TPB.quadfringeIN(bunch,kq)
			if(length == 0.):
				return
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]
				skew = skewArr[i]
				TPB.multpfringeIN(bunch,pole,kl/length,skew)

		def fringeOUT(node,paramsDict):
			usageOUT = node.getUsage()
			if(not usageOUT):
				return
			kq = node.getParam("kq")
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
				kl = klArr[i]
				skew = skewArr[i]
				TPB.multpfringeOUT(bunch,pole,kl/length,skew)

		self.setFringeFieldFunctionIN(fringeIN)
		self.setFringeFieldFunctionOUT(fringeOUT)

		self.setType("quad teapot")

	def _initialize(self, paramsDict):
		"""
		The  Quad Combined Function TEAPOT class implementation
		of the AccElement class _initialize(probe) method.
		"""
		nParts = self.getnParts()
		if(nParts < 2):
			msg = "The Quad Combined Function TEAPOT class instance should have more than 2 parts!"
			msg = msg + os.linesep
			msg = msg + "Method _initialize(self, paramsDict):"
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
		self.setPartLength(0,lengthIN)
		self.setPartLength(nParts - 1,lengthOUT)
		for i in xrange(nParts-2):
			self.setPartLength(i+1,lengthStep)

	def track(self, paramsDict):
		"""
		The Quad Combined Function TEAPOT  class implementation
		of the AccElement class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getPartLength(index)
		kq = self.getParam("kq")
		poleArr = self.getParam("poles")
		klArr = self.getParam("kls")
		skewArr = self.getParam("skews")
		bunch = paramsDict["bunch"]
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
			TPB.quad2(bunch, length)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			TPB.quad2(bunch, length)
			TPB.quad1(bunch, length, kq)
		return

class BendTEAPOT(NodeTEAPOT):
	"""
	Bend Combined Functions TEAPOT element.
	"""
	def __init__(self, name = "bend no name"):
		"""
		Constructor of the Bend Combined Functions TEAPOT element .
		"""
		NodeTEAPOT.__init__(self,name)
		self.setnParts(2)

		self.addParam("poles",[])
		self.addParam("kls",[])
		self.addParam("skews",[])

		self.addParam("ea1",0.)
		self.addParam("ea2",0.)
		self.addParam("rho",0.)
		self.addParam("theta",1.0e-36)

		def fringeIN(node,paramsDict):
			usageIN = node.getUsage()
			e = node.getParam("ea1")
			rho = node.getParam("rho")
			poleArr = node.getParam("poles")[:]
			klArr = node.getParam("kls")[:]
			skewArr = node.getParam("skews")[:]
			length = paramsDict["parentNode"].getLength()
			bunch = paramsDict["bunch"]
			nParts = paramsDict["parentNode"].getnParts()
			if(e != 0.):
				bunch = paramsDict["bunch"]
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
			usageOUT = node.getUsage()
			e = node.getParam("ea2")
			rho = node.getParam("rho")
			poleArr = node.getParam("poles")[:]
			klArr = node.getParam("kls")[:]
			skewArr = node.getParam("skews")[:]
			length = paramsDict["parentNode"].getLength()
			bunch = paramsDict["bunch"]
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

		self.setType("bend teapot")

	def _initialize(self, paramsDict):
		"""
		The  Bend Combined Functions TEAPOT class implementation of
		the AccElement class _initialize(probe) method.
		"""
		nParts = self.getnParts()
		if(nParts < 2):
			msg = "The Bend Combined Functions TEAPOT class instance should have more than 2 parts!"
			msg = msg + os.linesep
			msg = msg + "Method _initialize(self, paramsDict):"
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
		self.setPartLength(0,lengthIN)
		self.setPartLength(nParts - 1,lengthOUT)
		for i in xrange(nParts-2):
			self.setPartLength(i+1,lengthStep)

	def track(self, paramsDict):
		"""
		The Bend Combined Functions TEAPOT  class implementation of
		the AccElement class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getPartLength(index)
		poleArr = self.getParam("poles")
		klArr = self.getParam("kls")
		skewArr = self.getParam("skews")
		bunch = paramsDict["bunch"]
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

class RingRFTEAPOT(NodeTEAPOT):
	"""
	Ring RF TEAPOT element.
	"""
	def __init__(self, name = "ringRF no name"):
		"""
		Constructor of the Ring RF TEAPOT element.
		Harmonics numbers are 1,2,3, ...
		Voltages are in Volts.
		Phases are in radians.
		"""
		NodeTEAPOT.__init__(self,name)
		self.setnParts(1)
		self.addParam("harmonics",[])
		self.addParam("voltages",[])
		self.addParam("phases",[])
		self.setType("ringRF teapot")

	def addRF(self, harmonic, voltage, phase):
		"""
		Method. It adds the RF component with sertain
		harmonic, voltage, and phase.
		"""
		self.getParam("harmonics").append(harmonic)
		self.getParam("voltages").append(voltage)
		self.getParam("phases").append(phase)

	def _initialize(self, paramsDict):
		"""
		The Ring RF TEAPOT class implementation
		of the AccElement class _initialize(probe) method.
		"""
		nParts = self.getnParts()
		if(nParts > 1):
			msg = "The Ring RF TEAPOT class instance should not have more than 1 part!"
			msg = msg + os.linesep
			msg = msg + "Method _initialize(self, paramsDict):"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "nParts =" + str(nParts)
			msg = msg + os.linesep
			msg = msg + "lenght =" + str(self.getLength())
			orbitFinalize(msg)
		self.setPartLength(0,self.getLength())

	def track(self, paramsDict):
		"""
		The Ring RF TEAPOT class implementation
		of the AccElement class track(probe) method.
		"""
		harmArr = self.getParam("harmonics")
		voltArr = self.getParam("voltages")
		phaseArr = self.getParam("phases")
		bunch = paramsDict["bunch"]
		for i in xrange(len(harmArr)):
			TPB.ringRF(bunch,harmArr[i],voltArr[i],phaseArr[i])

class KickTEAPOT(NodeTEAPOT):
	"""
	Kick TEAPOT element.
	"""
	def __init__(self, name = "kick no name"):
		"""
		Constructor of the Kick TEAPOT element .
		"""
		NodeTEAPOT.__init__(self,name)
		self.setnParts(2)
		self.addParam("kx",0.)
		self.addParam("ky",0.)
		self.addParam("dE",0.)
		self.setType("kick teapot")

	def _initialize(self, paramsDict):
		"""
		The  Kicker TEAPOT class implementation of
		the AccElement class _initialize(probe) method.
		"""
		nParts = self.getnParts()
		if(nParts < 2):
			msg = "The Kick TEAPOT class instance should have more than 2 parts!"
			msg = msg + os.linesep
			msg = msg + "Method _initialize(self, paramsDict):"
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
		self.setPartLength(0,lengthIN)
		self.setPartLength(nParts - 1,lengthOUT)
		for i in xrange(nParts-2):
			self.setPartLength(i+1,lengthStep)

	def track(self, paramsDict):
		"""
		The Kick TEAPOT  class implementation of
		the AccElement class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getPartLength(index)
		kx = self.getParam("kx")/(nParts-1)
		ky = self.getParam("ky")/(nParts-1)
		dE = self.getParam("dE")/(nParts-1)
		bunch = paramsDict["bunch"]
		if(index == 0):
			TPB.drift(bunch, length)
			TPB.kick(bunch,kx,ky,dE)
			return
		if(index > 0 and index < (nParts-1)):
			TPB.drift(bunch, length)
			TPB.kick(bunch,kx,ky,dE)
			return
		if(index == (nParts-1)):
			TPB.drift(bunch, length)
		return

class TiltTEAPOT(BaseTEAPOT):
	"""
	The class to do tilt at the entrance of an TEAPOT element.
	"""
	def __init__(self, name = "tilt no name", angle = 0.):
		"""
		Constructor of the Tilt TEAPOT element.
		"""
		AccElement.__init__(self,name)
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

class FringeFieldTEAPOT(BaseTEAPOT):
	"""
	The class is a base class for the fringe field classes for others TEAPOT elements.
	"""
	def __init__(self,  parentNode,  trackFunction = None , name = "fringe field no name"):
		"""
		Constructor of the Fringe Field TEAPOT element.
		"""
		AccElement.__init__(self,name)
		self.setParams(parentNode.getParams())
		self.__trackFunc = trackFunction
		self.__usage = True
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