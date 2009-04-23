"""
Module. Includes classes for all TEAPOT elements. The approach is based
on the original ORBIT approach developed by J. Holmes.
"""

import sys
import os
import math

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

# import the function that creates multidimensional arrays
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer

# import the MAD parser to construct lattices of TEAPOT elements.
from orbit.parsers.mad_parser import MAD_Parser, MAD_LattElement

"""
Drift
Bend
Quad
Multipole
Solenoid
Kicker
RingRF
"""

class TEAPOT_Lattice(AccLattice):
	"""
	The subclass of the AccLattice class. Shell class for the TEAPOT nodes collection.
	TEAPOT has the ability to read MAD files.
	"""
	def __init__(self, name = "no name"):
		AccLattice.__init__(self,name)

	def readMAD(self, mad_file_name, lineName):
		"""
		It creates the teapot lattice from MAD file.
		"""
		parser = MAD_Parser()
		parser.parse(mad_file_name)
		accLines = parser.getMAD_LinesDict()
		if(not accLines.has_key(lineName)):
			print "==============================="
			print "MAD file: ", mad_file_name
			print "Can not find accelerator line: ", lineName
			print "STOP."
			sys.exit(1)
		# accelerator lines and elements from mad_parser package
		accMAD_Line = accLines[lineName]
		self.setName(lineName)
		accMADElements = accMAD_Line.getElements()
		# make TEAPOT lattice elements by using TEAPOT
		# element factory
		for madElem in accMADElements:
			elem = _teapotFactory.getElement(madElem)
			self.addNode(elem)
		self.initialize()

	def initialize(self):
		AccLattice.initialize(self)
		#set up ring length for RF nodes
		ringRF_Node = RingRFTEAPOT()
		for node in self.getNodes():
			if(node.getType() == ringRF_Node.getType()):
				node.getParamsDict()["ring_lenght"] = self.getLength()		
	

	def getSubLattice(self, index_start = -1, index_stop = -1,):
		"""
		It returns the new TEAPOT_Lattice with children with indexes 
		between index_start and index_stop inclusive
		"""
		return self._getSubLattice(TEAPOT_Lattice(),index_start,index_stop)

	def trackBunch(self, bunch, paramsDict = {}, actionContainer = None):
		"""
		It tracks the bunch through the lattice.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Bunch Tracking")
		paramsDict["bunch"] = bunch
		
		def track(paramsDict):
			node = paramsDict["node"]
			node.track(paramsDict)
			
		actionContainer.addAction(track, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict)
		actionContainer.removeAction(track, AccActionsContainer.BODY)
		

class _teapotFactory:
	"""
	Class. Factory to produce TEAPOT accelerator elements.
	"""
	def __init__(self):
		pass

	def getElement(self, madElem):
		"""
		Method. Produces TEAPOT accelerator elements.
		"""
		# madElem = MAD_LattElement(" "," ")
		params_ini = madElem.getParameters()
		params = {}
		for par_key in params_ini.iterkeys():
			params[par_key.lower()] = params_ini[par_key]
		# Length parameter
		length = 0.
		if(params.has_key("l")): length = params["l"]
		# TILT parameter - if it is there tilt = True,
		# and if it has value tiltAngle = value
		tilt = False
		tiltAngle = None
		if(params.has_key("tilt")):
			tilt = True
			if(params["tilt"] != None):
				tiltAngle = params["tilt"]
		# ============DRIFT Element ==================================
		elem = None
		if(madElem.getType().lower() == "drift"):
			elem = DriftTEAPOT(madElem.getName())
		# ============Dipole Element SBEND or RBEND ===================
		if(madElem.getType().lower()  == "sbend" or \
			 madElem.getType().lower()  == "rbend"):
			elem = BendTEAPOT(madElem.getName())

			theta = 0.
			if(params.has_key("angle")):
				theta = params["angle"]

			ea1 = 0.
			if(params.has_key("e1")):
				ea1 = params["e1"]

			ea2 = 0.
			if(params.has_key("e2")):
				ea2 = params["e2"]

			if(madElem.getType().lower() == "rbend" ):
				length = length*(theta/2.0)/math.sin(theta/2.0)
				ea1 = ea1 + theta/2.0
				ea2 = ea2 + theta/2.0

			elem.addParam("theta",theta)
			elem.addParam("ea1",ea1)
			elem.addParam("ea2",ea2)
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
		# ===========QUAD quadrupole element =====================
		if(madElem.getType().lower()  == "quad" or \
			 madElem.getType().lower()  == "quadrupole"):
			elem = QuadTEAPOT(madElem.getName())
			kq = 0.
			if(params.has_key("k1")):
				kq = params["k1"]
			elem.addParam("kq",kq)
			if(tilt):
				if(tiltAngle == None):
					tiltAngle = math.math.pi/4.0
		# ===========Sextupole element =====================
		if(madElem.getType().lower()  == "sextupole"):
			elem = MultipoleTEAPOT(madElem.getName())
			k2 = 0.
			if(params.has_key("k2")):
				k2 = params["k2"]
			params["k2l"] = length*k2
			if(tilt):
				if(tiltAngle == None):
					tiltAngle = math.math.pi/6.0
		# ===========Octupole element =====================
		if(madElem.getType().lower()  == "octupole"):
			elem = MultipoleTEAPOT(madElem.getName())
			k3 = 0.
			if(params.has_key("k3")):
				k3 = params["k3"]
			params["k2l"] = length*k3
			if(tilt):
				if(tiltAngle == None):
					tiltAngle = math.math.pi/8.0
		# ===========Multipole element =====================
		if(madElem.getType().lower()  == "multipole"):
			elem = MultipoleTEAPOT(madElem.getName())
		# ===========Solenoid element ======================
		if(madElem.getType().lower()  == "solenoid"):
			elem = SolenoidTEAPOT(madElem.getName())
			ks = 0.
			if(params.has_key("ks")):
				ks = params["ks"]
			elem.addParam("B",ks)
		# ===========Kicker element ======================
		if(madElem.getType().lower()  == "kicker"):
			elem = KickTEAPOT(madElem.getName())
			hkick = 0.
			if(params.has_key("hkick")):
				hkick = params["hkick"]
			vkick = 0.
			if(params.has_key("vkick")):
				hkick = params["vkick"]
			elem.addParam("kx",hkick)
			elem.addParam("ky",vkick)
		# ===========HKicker element ======================
		if(madElem.getType().lower()  == "hkicker" or \
			 madElem.getType() .lower() == "hkick"):
			elem = KickTEAPOT(madElem.getName())
			hkick = 0.
			if(params.has_key("hkick")):
				hkick = params["hkick"]
			elem.addParam("kx",hkick)
		# ===========VKicker element ======================
		if(madElem.getType().lower()  == "vkicker" or \
			 madElem.getType().lower()  == "vkick"):
			elem = KickTEAPOT(madElem.getName())
			vkick = 0.
			if(params.has_key("vkick")):
				hkick = params["vkick"]
			elem.addParam("ky",vkick)
		# ===========RF Cavity element ======================
		if(madElem.getType().lower() == "rfcavity"):
			elem = RingRFTEAPOT(madElem.getName())
			# the MAD RF element has a L parameter,
			# but it does not contribute in the ring length!
			"""
			if(length != 0.):
				drft_1 = DriftTEAPOT(madElem.getName()+"_drift")
				drft_2 = DriftTEAPOT(madElem.getName()+"_drift")
				drft_1.setLength(length/2.0)
				drft_2.setLength(length/2.0)
				elem.addChildNode(drft_1,AccNode.ENTRANCE)
				elem.addChildNode(drft_2,AccNode.EXIT)
			"""
			volt = 0.
			if(params.has_key("volt")):
				volt = params["volt"]
			harmon = 1
			if(params.has_key("harmon")):
				harmon = params["harmon"]
			phase_s = 0.
			if(params.has_key("lag")):
				phase_s = (-1.0) * params["lag"]*2.0*math.pi
			elem.addRF(harmon,volt,phase_s)
		# ==========Others elements such as markers,monitor,rcollimator
		if(madElem.getType().lower() == "marker" or \
			 madElem.getType().lower() == "monitor" or \
			 madElem.getType().lower() == "rcolimator"):
			elem = NodeTEAPOT(madElem.getName())
		# ------------------------------------------------
		# ready to finish
		# ------------------------------------------------
		if(elem == None):
			print "======== Can not create elementwith type:",madElem.getType()
			print "You have to fix the _teapotFactory class."
			print "Stop."
			sys.exit(1)
		# set length
		elem.setLength(length)
		# set tilt angle
		if(tilt == True and tiltAngle != None):
			elem.setTiltAngle(tiltAngle)
		# set K1L,K2L,K3L,... and  T1L,T2L,T3L,...
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
			elem.addParam("poles",poles)
			elem.addParam("kls",kls)
			elem.addParam("skews",skews)
		return elem

	getElement = classmethod(getElement)


class BaseTEAPOT(AccNode):
	""" The base abstract class of the TEAPOT accelerator elements hierarchy. """
	def __init__(self, name = "no name"):
		"""
		Constructor. Creates the base TEAPOT element. This is a superclass for all TEAPOT elements.
		"""
		AccNode.__init__(self,name)
		self.setType("base teapot")
		
	def trackBunch(self, bunch, paramsDict = {}, actionContainer = None):
		"""
		It tracks the bunch through the BaseTEAPOT instance.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Bunch Tracking")
		paramsDict["bunch"] = bunch
		
		def track(paramsDict):
			node = paramsDict["node"]
			node.track(paramsDict)
			
		actionContainer.addAction(track, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict)
		actionContainer.removeAction(track, AccActionsContainer.BODY)		
		
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the element. Each element 
		should implement this method.
		"""
		pass	

class NodeTEAPOT(BaseTEAPOT):
	def __init__(self, name = "no name"):
		"""
		Constructor. Creates the real TEAPOT element. This is a superclass for all real TEAPOT elements.
		The term real means that element is included in TEAPOT list of elements.
		For instance TILT is not included.
		"""
		BaseTEAPOT.__init__(self,name)
		self.__tiltNodeIN  = TiltTEAPOT()
		self.__tiltNodeOUT = TiltTEAPOT()
		self.__fringeFieldIN = FringeFieldTEAPOT(self)
		self.__fringeFieldOUT = FringeFieldTEAPOT(self)
		self.addChildNode(self.__tiltNodeIN,AccNode.ENTRANCE)
		self.addChildNode(self.__fringeFieldIN,AccNode.ENTRANCE)
		self.addChildNode(self.__fringeFieldOUT,AccNode.EXIT)
		self.addChildNode(self.__tiltNodeOUT,AccNode.EXIT)
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

	def getNodeFringeFieldIN(self):
		"""
		Returns the FringeFieldTEAPOT instance before this TEAPOT node
		"""
		return self.__fringeFieldIN 
 
	def getNodeFringeFieldOUT(self):
		"""
		Returns the FringeFieldTEAPOT instance after this TEAPOT node
		"""
		return self.__fringeFieldOUT 
		
	def getNodeTiltIN(self):
		"""
		Returns the TiltTEAPOT instance before this TEAPOT node
		"""
		return self.__tiltNodeIN 
 
	def getNodeTiltOUT(self):
		"""
		Returns the  TiltTEAPOT instance after this TEAPOT node
		"""
		return self.__tiltNodeOUT		
		
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
		Constructor. Creates the Drift TEAPOT element.
		"""
		NodeTEAPOT.__init__(self,name)
		self.setType("drift teapot")

	def track(self, paramsDict):
		"""
		The drift class implementation of the AccNode class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		TPB.drift(bunch, length)

class SolenoidTEAPOT(NodeTEAPOT):
	"""
	Solenoid TEAPOT element.
	"""
	def __init__(self, name = "solenoid no name"):
		"""
		Constructor. Creates the Solenoid TEAPOT element .
		"""
		NodeTEAPOT.__init__(self,name)
		self.setType("solenoid teapot")

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

	def track(self, paramsDict):
		"""
		The Solenoid TEAPOT  class implementation of the AccNode class track(probe) method.
		"""
		index = self.getActivePartIndex()
		length = self.getLength(index)
		bunch = paramsDict["bunch"]
		B = node.getParam("B")
		TPB.soln(bunch, length, B)

class MultipoleTEAPOT(NodeTEAPOT):
	"""
	Multipole Combined Function TEAPOT element.
	"""
	def __init__(self, name = "multipole no name"):
		"""
		Constructor. Creates the Multipole Combined Function TEAPOT element.
		"""
		NodeTEAPOT.__init__(self,name)
		self.addParam("poles",[])
		self.addParam("kls",[])
		self.addParam("skews",[])
		
		self.setnParts(2)
		
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

	def initialize(self):
		"""
		The  Multipole Combined Function TEAPOT class
		implementation of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts < 2):
			msg = "The Multipole TEAPOT class instance should have more than 2 parts!"
			msg = msg + os.linesep
			msg = msg + "Method initialize():"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "nParts =" + str(self.getnParts())
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
		The Multipole Combined Function TEAPOT  class
		implementation of the AccNode class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
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
		Constructor. Creates the Quad Combined Function TEAPOT element .
		"""
		NodeTEAPOT.__init__(self,name)

		self.addParam("kq",0.)
		self.addParam("poles",[])
		self.addParam("kls",[])
		self.addParam("skews",[])
		self.setnParts(2)

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
		self.getNodeTiltIN().setType("quad tilt in")
		self.getNodeTiltOUT().setType("quad tilt out")
		self.getNodeFringeFieldIN().setType("quad fringe in")
		self.getNodeFringeFieldOUT().setType("quad fringe out")

		self.setType("quad teapot")

	def initialize(self):
		"""
		The  Quad Combined Function TEAPOT class implementation
		of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts < 2):
			msg = "The Quad Combined Function TEAPOT class instance should have more than 2 parts!"
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
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
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
		Constructor. Creates the Bend Combined Functions TEAPOT element .
		"""
		NodeTEAPOT.__init__(self,name)
		self.addParam("poles",[])
		self.addParam("kls",[])
		self.addParam("skews",[])

		self.addParam("ea1",0.)
		self.addParam("ea2",0.)
		self.addParam("rho",0.)
		self.addParam("theta",1.0e-36)
		
		self.setnParts(2)
		
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
		the AccNode class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
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
	def __init__(self, name = "RingRF no name"):
		"""
		Constructor. Creates the Ring RF TEAPOT element.
		Harmonics numbers are 1,2,3, ...
		Voltages are in Volts.
		Phases are in radians.
		"""
		NodeTEAPOT.__init__(self,name)
		self.addParam("harmonics",[])
		self.addParam("voltages",[])
		self.addParam("phases",[])
		self.addParam("ring_length",0.)
		self.setType("RingRF teapot")
		
		self.setnParts(1)
		

	def addRF(self, harmonic, voltage, phase):
		"""
		Method. It adds the RF component with sertain
		harmonic, voltage, and phase.
		"""
		self.getParam("harmonics").append(harmonic)
		self.getParam("voltages").append(voltage)
		self.getParam("phases").append(phase)

	def initialize(self):
		"""
		The Ring RF TEAPOT class implementation
		of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts > 1):
			msg = "The Ring RF TEAPOT class instance should not have more than 1 part!"
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

	def track(self, paramsDict):
		"""
		The Ring RF TEAPOT class implementation
		of the AccNode class track(probe) method.
		"""
		harmArr = self.getParam("harmonics")
		voltArr = self.getParam("voltages")
		phaseArr = self.getParam("phases")
		ring_length = self.getParam("ring_length")
		bunch = paramsDict["bunch"]
		for i in range(len(harmArr)):
			TPB.RingRF(bunch,ring_length,harmArr[i],voltArr[i],phaseArr[i])

class KickTEAPOT(NodeTEAPOT):
	"""
	Kick TEAPOT element.
	"""
	def __init__(self, name = "kick no name"):
		"""
		Constructor. Creates the Kick TEAPOT element .
		"""
		NodeTEAPOT.__init__(self,name)
		self.addParam("kx",0.)
		self.addParam("ky",0.)
		self.addParam("dE",0.)
		self.setType("kick teapot")
		
		self.setnParts(2)
		
	def initialize(self):
		"""
		The  Kicker TEAPOT class implementation of
		the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts < 2):
			msg = "The Kick TEAPOT class instance should have more than 2 parts!"
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
		The Kick TEAPOT  class implementation of
		the AccNode class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
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
		Constructor. Creates the Tilt TEAPOT element.
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

class FringeFieldTEAPOT(BaseTEAPOT):
	"""
	The class is a base class for the fringe field classes for others TEAPOT elements.
	"""
	def __init__(self,  parentNode,  trackFunction = None , name = "fringe field no name"):
		"""
		Constructor. Creates the Fringe Field TEAPOT element.
		"""
		AccNode.__init__(self,name)
		self.setParamsDict(parentNode.getParamsDict())
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
