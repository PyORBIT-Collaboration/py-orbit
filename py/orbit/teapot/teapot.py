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
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import the MAD parser to construct lattices of TEAPOT elements.
from orbit.parsers.mad_parser import MAD_Parser, MAD_LattElement
# import the MADX parser to construct lattices of TEAPOT elements.
from orbit.parsers.madx_parser import MADX_Parser, MADX_LattElement

# import aperture
from aperture import Aperture

# monitor
from bunch import BunchTwissAnalysis


"""
Drift
Bend
Quad
Multipole
Solenoid
Kicker
RingRF
monitor
"""


class TEAPOT_Lattice(AccLattice):
	"""
	The subclass of the AccLattice class. Shell class for the TEAPOT nodes collection.
	TEAPOT has the ability to read MAD files.
	"""
	def __init__(self, name = "no name"):
		AccLattice.__init__(self,name)
		self.useCharge = 1

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
			elems = _teapotFactory.getElements(madElem)
			for elem in elems:
				self.addNode(elem)
		self._addChildren()
		self.initialize()

	def readMADX(self, madx_file_name, seqName):
		"""
			It creates the teapot lattice from MAD file.
		"""
		parser = MADX_Parser()
		parser.parse(madx_file_name)
		
		if(not seqName == parser.getSequenceName()):
			print "==============================="
			print "MADX file: ", madx_file_name
			print "Can not find accelerator sequence: ", seqName
			print "STOP."
			sys.exit(1)
		
		self.setName(parser.getSequenceName())
		accMADElements = parser.getSequenceList()
		# make TEAPOT lattice elements by using TEAPOT
		# element factory
		for madElem in accMADElements:
			#print madElem.getParameters()
			elems = _teapotFactory.getElements(madElem)
			for elem in elems:
				self.addNode(elem)
		self._addChildren()
		self.initialize()

	def _addChildren(self):
		"""
		It does nothing here. It is a place holder for TEAPOT_Ring class.
		"""
		pass

	def initialize(self):
		AccLattice.initialize(self)
		#set up ring length for RF nodes
		ringRF_Node = RingRFTEAPOT()
		for node in self.getNodes():
			if(node.getType() == ringRF_Node.getType()):
				node.getParamsDict()["ring_length"] = self.getLength()

	def getSubLattice(self, index_start = -1, index_stop = -1,):
		"""
		It returns the new TEAPOT_Lattice with children with indexes 
		between index_start and index_stop inclusive
		"""
		new_teapot_lattice = self._getSubLattice(TEAPOT_Lattice(),index_start,index_stop)
		new_teapot_lattice.setUseRealCharge(self.getUseRealCharge())
		return new_teapot_lattice

	def trackBunch(self, bunch, paramsDict = {}, actionContainer = None, index_start = -1, index_stop = -1):
		"""
		It tracks the bunch through the lattice. Indexes index_start and index_stop are inclusive.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Bunch Tracking")
		paramsDict["bunch"] = bunch
		paramsDict["useCharge"] = self.useCharge
		
		def track(paramsDict):
			node = paramsDict["node"]
			node.track(paramsDict)
			
		actionContainer.addAction(track, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict,index_start,index_stop)
		actionContainer.removeAction(track, AccActionsContainer.BODY)

	def setUseRealCharge(self, useCharge = 1):
		""" If useCharge != 1 the trackBunch(...) method will assume the charge = +1 """ 
		self.useCharge = useCharge
	
	def getUseRealCharge(self):
		""" If useCharge != 1 the trackBunch(...) method will assume the charge = +1 """ 
		return self.useCharge


class TEAPOT_Ring(TEAPOT_Lattice):
	"""
		The subclass of the TEAPOT_Lattice class. Shell class for the TEAPOT nodes 
		collection for rings. TEAPOT has the ability to read MAD files.
		"""
	def __init__(self, name = "no name"):
		TEAPOT_Lattice.__init__(self,name)	

	def _addChildren(self):
		"""
		Adds Bunch wrapping nodes to all lattice nodes of 1st level.
		These wrapping nodes will move particles from head ( tail ) to tail ( head)
		if the longitudinal positions too big / small(negative). 
		"""
		TEAPOT_Lattice.initialize(self)
		for node in self.getNodes():
			length = node.getLength()
			if(length > 0.):
				bunchwrapper = BunchWrapTEAPOT(node.getName()+":Bunch_Wrap:Exit")
				bunchwrapper.getParamsDict()["ring_length"] = self.getLength()
				node.addChildNode(bunchwrapper, AccNode.EXIT)				
		#---- adding turn counter node at the end of lattice
		turn_counter = TurnCounterTEAPOT()
		self.getNodes().append(turn_counter)

class _teapotFactory:
	"""
	Class. Factory to produce TEAPOT accelerator elements.
	"""
	def __init__(self):
		pass

	def getElements(self, madElem):
		"""
		Method produces TEAPOT accelerator elements. It returns the array of TEAPOT nodes.
		Usually there is only one such element, but sometimes (RF with a non zero length) 
		there could be more than one element in the returned array 
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
		# ============Aperture Element ==================================
		if(madElem.getType().lower() == "aperture"):
			elem = ApertureTEAPOT(madElem.getName())
			if(params.has_key("apertype")):
				elem.addParam("aperture", params["aperture"])
				elem.addParam("apertype", params["apertype"])
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
			params["k3l"] = length*k3
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
				vkick = params["vkick"]
			elem.addParam("kx",hkick)
			elem.addParam("ky",vkick)
		# ===========HKicker element ======================
		if(madElem.getType().lower() == "hkicker" or \
                   madElem.getType().lower() == "hkick"):
			elem = KickTEAPOT(madElem.getName())
			hkick = 0.
			if(params.has_key("hkick")):
				hkick = params["hkick"]
			elem.addParam("kx",hkick)
		# ===========VKicker element ======================
		if(madElem.getType().lower() == "vkicker" or \
                   madElem.getType().lower() == "vkick"):
			elem = KickTEAPOT(madElem.getName())
			vkick = 0.
			if(params.has_key("vkick")):
				vkick = params["vkick"]
			elem.addParam("ky",vkick)
		# ===========RF Cavity element ======================
		if(madElem.getType().lower() == "rfcavity"):
			elem = RingRFTEAPOT(madElem.getName())
			drft_1 = DriftTEAPOT(madElem.getName()+"_drift1")
			drft_2 = DriftTEAPOT(madElem.getName()+"_drift2")
			drft_1.setLength(length/2.0)
			drft_2.setLength(length/2.0)
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
			return [drft_1,elem,drft_2]
		# ==========Others elements such as markers,monitor,rcollimator
		if(madElem.getType().lower() == "marker" or \
			 #madElem.getType().lower() == "monitor" or \
			 #madElem.getType().lower() == "hmonitor" or \
			 #madElem.getType().lower() == "vmonitor" or \
			 madElem.getType().lower() == "rcolimator"):
			elem = NodeTEAPOT(madElem.getName())
		if(madElem.getType().lower() =="monitor"):
			elem = 	MonitorTEAPOT(madElem.getName())
			drft_1 = DriftTEAPOT(madElem.getName()+"_drift1")
			drft_2 = DriftTEAPOT(madElem.getName()+"_drift2")
			drft_1.setLength(length/2.0)
			drft_2.setLength(length/2.0)
			xAvg = 0.0
			yAvg = 0.0
			elem.addParam("xAvg",xAvg)
			elem.addParam("yAvg",yAvg)
			if length > 0:
				return [drft_1, elem, drft_2]
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
		if(params.has_key("kls")) and (params.has_key("poles") and params.has_key("skews")):
			elem.addParam("kls",params["kls"])
			elem.addParam("poles",params["poles"])
			elem.addParam("skews",params["skews"])
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
		return [elem]

	getElements = classmethod(getElements)

class BaseTEAPOT(AccNodeBunchTracker):
	""" The base abstract class of the TEAPOT accelerator elements hierarchy. """
	def __init__(self, name = "no name"):
		"""
		Constructor. Creates the base TEAPOT element. This is a superclass for all TEAPOT elements.
		"""
		AccNodeBunchTracker.__init__(self,name)
		self.setType("base teapot")

class TurnCounterTEAPOT(BaseTEAPOT):
	def __init__(self, name = "TurnCounter"):
		"""
		Constructor. Creates the TEAPOT for turn count in the Ring lattice.
		"""
		BaseTEAPOT.__init__(self,name)
		self.setType("turn counter")
		
	def track(self, paramsDict):
		"""
		The Turn Counter class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		bunch = paramsDict["bunch"]
		if(bunch.hasBunchAttrInt("TurnNumber") != 0):
			turn = bunch.bunchAttrInt("TurnNumber")
			bunch.bunchAttrInt("TurnNumber",turn + 1)

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
		self.__tiltNodeIN.setName(name+"_tilt_in")
		self.__tiltNodeOUT.setName(name+"_tilt_out")
		self.__fringeFieldIN.setName(name+"_fringe_in")
		self.__fringeFieldOUT	.setName(name+"_fringe_out")	
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
		self.setParam("tilt", angle)
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
		The drift class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		TPB.drift(bunch, length)


class ApertureTEAPOT(NodeTEAPOT):
	"""
	Aperture TEAPOT element.
	"""
	def __init__(self, name = "aperture no name"):
		"""
		Constructor. Creates the aperutre element.
		"""
		NodeTEAPOT.__init__(self,name)
		self.setType("aperture")
		self.addParam("aperture", [])
		self.addParam("apertype", 0.0)
		
	def initialize(self):
	
		shape = self.getParam("apertype")
		dim = self.getParam("aperture")
		if len(dim) > 0:
			if shape == 1:
				self.aperture = Aperture(shape, dim[0], 0.0, 0.0, 0.0)
			if shape == 2:
				self.aperture = Aperture(shape, dim[0], dim[1], 0.0, 0.0)
			if shape == 3:
				self.aperture = Aperture(shape, dim[0], dim[1], 0.0, 0.0)

	def track(self, paramsDict):
		"""
		The aperture class implementation of the ApertueNode.
		"""
		bunch = paramsDict["bunch"]
		lostbunch = paramsDict["lostbunch"]
		self.aperture.checkBunch(bunch,lostbunch)

class MonitorTEAPOT(NodeTEAPOT):
	"""
	Monitor TEAPOT element.
	"""
	def __init__(self, name = "Monitor no name"):
		"""
		Constructor. Creates the aperutre element.
		"""
		NodeTEAPOT.__init__(self,name)
		self.setType("monitor teapot")
		self.twiss =  BunchTwissAnalysis()
			
	def track(self, paramsDict):
		"""
		The bunchtuneanalysis-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		self.twiss.analyzeBunch(bunch)
		self.addParam("xAvg",self.twiss.getAverage(0))
		self.addParam("xpAvg",self.twiss.getAverage(1))
		self.addParam("yAvg",self.twiss.getAverage(2))
		self.addParam("ypAvg",self.twiss.getAverage(3))


class BunchWrapTEAPOT(NodeTEAPOT):
	"""
		Drift TEAPOT element.
		"""
	def __init__(self, name = "drift no name"):
		"""
			Constructor. Creates the Bunch wrapper TEAPOT element used in Ring lattices.
		"""
		NodeTEAPOT.__init__(self,name)
		self.setType("bunch_wrap_teapot")
		self.addParam("ring_length",0.)
			   
	def track(self, paramsDict):
		"""
			The bunch wrap class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		length = self.getParam("ring_length")
		TPB.wrapbunch(bunch, length)


class SolenoidTEAPOT(NodeTEAPOT):
	"""
	Solenoid TEAPOT element.
	"""
	def __init__(self, name = "solenoid no name"):
		"""
		Constructor. Creates the Solenoid TEAPOT element.
		"""
		NodeTEAPOT.__init__(self, name)
		self.setType("solenoid teapot")
		self.addParam("B", 0.)
		self.waveform = None

	def track(self, paramsDict):
		"""
		The Solenoid TEAPOT class implementation of the
                AccNodeBunchTracker class track(probe) method.
		"""
		index = self.getActivePartIndex()
		length = self.getLength(index)
		strength = 1.0
		if(self.waveform):
			strength = self.waveform.getStrength()
		B = strength * self.getParam("B")
		bunch = paramsDict["bunch"]
		useCharge = 1
		if(paramsDict.has_key("useCharge")): \
                  useCharge = paramsDict["useCharge"]
		TPB.soln(bunch, length, B, useCharge)

	def setWaveform(self, waveform):
		"""
		Sets the time dependent waveform function
		"""
		self.waveform = waveform


class MultipoleTEAPOT(NodeTEAPOT):
	"""
	Multipole Combined Function TEAPOT element.
	"""
	def __init__(self, name = "multipole no name"):
		"""
		Constructor. Creates the Multipole
                Combined Function TEAPOT element.
		"""
		NodeTEAPOT.__init__(self,name)
		self.addParam("poles", [])
		self.addParam("kls", [])
		self.addParam("skews", [])
		self.setnParts(2)
		self.waveform = None

                def fringeIN(node, paramsDict):
			usageIN = node.getUsage()
			if(not usageIN):
				return
			length = paramsDict["parentNode"].getLength()
			if(length == 0.):
				return
                        strength = 1.0
                        if(self.waveform):
                                strength = self.waveform.getStrength()
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			bunch = paramsDict["bunch"]
			useCharge = 1
			if(paramsDict.has_key("useCharge")): \
                          useCharge = paramsDict["useCharge"]
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = strength * klArr[i]
				skew = skewArr[i]
				TPB.multpfringeIN(bunch, pole, kl/length, \
                                                  skew, useCharge)

		def fringeOUT(node,paramsDict):
			usageOUT = node.getUsage()
			if(not usageOUT):
				return
			length = paramsDict["parentNode"].getLength()
			if(length == 0.):
				return
                        strength = 1.0
                        if(self.waveform):
                                strength = self.waveform.getStrength()
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			bunch = paramsDict["bunch"]
			useCharge = 1
			if(paramsDict.has_key("useCharge")): \
                          useCharge = paramsDict["useCharge"]
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = strength * klArr[i]
				skew = skewArr[i]
				TPB.multpfringeOUT(bunch, pole, kl/length, \
                                                   skew, useCharge)

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
		lengthIN = (self.getLength() / (nParts - 1)) / 2.0
		lengthOUT = (self.getLength() / (nParts - 1)) / 2.0
		lengthStep = lengthIN + lengthOUT
		self.setLength(lengthIN, 0)
		self.setLength(lengthOUT, nParts - 1)
		for i in xrange(nParts - 2):
			self.setLength(lengthStep, i + 1)

	def track(self, paramsDict):
		"""
		The Multipole Combined Function TEAPOT  class
		implementation of the AccNodeBunchTracker class
                track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
		strength = 1.0
		if(self.waveform):
			strength = self.waveform.getStrength()
		poleArr = self.getParam("poles")
		klArr = self.getParam("kls")
		skewArr = self.getParam("skews")
		bunch = paramsDict["bunch"]
		useCharge = 1
		if(paramsDict.has_key("useCharge")): \
                  useCharge = paramsDict["useCharge"]
		if(index == 0):
			TPB.drift(bunch, length)
			return
		if(index > 0 and index < (nParts - 1)):
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = strength * klArr[i] / (nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch, pole, kl, skew, useCharge)
			TPB.drift(bunch, length)
			return
		if(index == (nParts-1)):
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = strength * klArr[i] / (nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch, pole, kl, skew, useCharge)
			TPB.drift(bunch, length)
		return

	def setWaveform(self, waveform):
		"""
		Sets the time dependent waveform function
		"""
		self.waveform = waveform


class QuadTEAPOT(NodeTEAPOT):
	"""
	Quad Combined Function TEAPOT element.
	"""
	def __init__(self, name = "quad no name"):
		"""
		Constructor. Creates the Quad
                Combined Function TEAPOT element.
		"""
		NodeTEAPOT.__init__(self,name)

		self.addParam("kq", 0.)
		self.addParam("poles", [])
		self.addParam("kls", [])
		self.addParam("skews", [])
		self.setnParts(2)
		self.waveform = None

		def fringeIN(node,paramsDict):
			usageIN = node.getUsage()
			if(not usageIN):
				return
                        strength = 1.0
                        if(self.waveform):
                                strength = self.waveform.getStrength()
			kq = strength * node.getParam("kq")
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			length = paramsDict["parentNode"].getLength()
			bunch = paramsDict["bunch"]
			useCharge = 1
			if(paramsDict.has_key("useCharge")): \
                          useCharge = paramsDict["useCharge"]
			TPB.quadfringeIN(bunch, kq, useCharge)
			if(length == 0.):
				return
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = strength * klArr[i]
				skew = skewArr[i]
				TPB.multpfringeIN(bunch, pole, kl/length, \
                                                  skew, useCharge)

		def fringeOUT(node,paramsDict):
			usageOUT = node.getUsage()
			if(not usageOUT):
				return
                        strength = 1.0
                        if(self.waveform):
                                strength = self.waveform.getStrength()
			kq = strength * node.getParam("kq")
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			length = paramsDict["parentNode"].getLength()
			bunch = paramsDict["bunch"]
			useCharge = 1
			if(paramsDict.has_key("useCharge")): \
                          useCharge = paramsDict["useCharge"]
			TPB.quadfringeOUT(bunch, kq, useCharge)
			if(length == 0.):
				return
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = strength * klArr[i]
				skew = skewArr[i]
				TPB.multpfringeOUT(bunch, pole, kl/length, \
                                                   skew,useCharge)

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
		lengthIN = (self.getLength() / (nParts - 1)) / 2.0
		lengthOUT = (self.getLength() / (nParts - 1)) / 2.0
		lengthStep = lengthIN + lengthOUT
		self.setLength(lengthIN, 0)
		self.setLength(lengthOUT, nParts - 1)
		for i in xrange(nParts - 2):
			self.setLength(lengthStep, i + 1)

	def track(self, paramsDict):
		"""
		The Quad Combined Function TEAPOT  class implementation
		of the AccNodeBunchTracker class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
                strength = 1.0
                if(self.waveform):
                        strength = self.waveform.getStrength()
		kq = strength * self.getParam("kq")
		poleArr = self.getParam("poles")
		klArr = self.getParam("kls")
		skewArr = self.getParam("skews")
		bunch = paramsDict["bunch"]
		useCharge = 1
		if(paramsDict.has_key("useCharge")): \
                  useCharge = paramsDict["useCharge"]
		if(index == 0):
			TPB.quad1(bunch, length, kq, useCharge)
			return
		if(index > 0 and index < (nParts - 1)):
			TPB.quad2(bunch, length/2.0)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = strength * klArr[i] / (nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch, pole, kl, skew, useCharge)
			TPB.quad2(bunch, length / 2.0)
			TPB.quad1(bunch, length, kq, useCharge)
			return
		if(index == (nParts - 1)):
			TPB.quad2(bunch, length)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = strength * klArr[i] / (nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch, pole, kl, skew, useCharge)
			TPB.quad2(bunch, length)
			TPB.quad1(bunch, length, kq, useCharge)
		return

	def setWaveform(self, waveform):
		"""
		Sets the time dependent waveform function
		"""
		self.waveform = waveform

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
			useCharge = 1
			if(paramsDict.has_key("useCharge")): useCharge = paramsDict["useCharge"]
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
							TPB.multpfringeIN(bunch,pole,kl,skew,useCharge)
					frinout = 1
					TPB.wedgerotate(bunch, e, frinout)
				TPB.wedgebendCF(bunch, e, inout, rho, len(poleArr), poleArr, klArr, skewArr, nParts - 1, useCharge)
			else:
				if(usageIN):
					TPB.bendfringeIN(bunch, rho)
					if(length != 0.):
						for i in xrange(len(poleArr)):
							pole = poleArr[i]
							kl = klArr[i]/length
							skew = skewArr[i]
							TPB.multpfringeIN(bunch,pole,kl,skew,useCharge)

		def fringeOUT(node,paramsDict):
			usageOUT = node.getUsage()
			e = node.getParam("ea2")
			rho = node.getParam("rho")
			poleArr = node.getParam("poles")[:]
			klArr = node.getParam("kls")[:]
			skewArr = node.getParam("skews")[:]
			length = paramsDict["parentNode"].getLength()
			bunch = paramsDict["bunch"]
			useCharge = 1
			if(paramsDict.has_key("useCharge")): useCharge = paramsDict["useCharge"]
			nParts = paramsDict["parentNode"].getnParts()
			if(e != 0.):
				inout = 1
				TPB.wedgebendCF(bunch, e, inout, rho, len(poleArr), poleArr, klArr, skewArr, nParts - 1, useCharge)
				if(usageOUT):
					frinout = 0
					TPB.wedgerotate(bunch, -e, frinout)
					TPB.bendfringeOUT(bunch, rho)
					if(length != 0.):
						for i in xrange(len(poleArr)):
							pole = poleArr[i]
							kl = klArr[i]/length
							skew = skewArr[i]
							TPB.multpfringeOUT(bunch,pole,kl,skew,useCharge)
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
							TPB.multpfringeOUT(bunch,pole,kl,skew,useCharge)

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
		the AccNodeBunchTracker class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
		poleArr = self.getParam("poles")
		klArr = self.getParam("kls")
		skewArr = self.getParam("skews")
		bunch = paramsDict["bunch"]
		useCharge = 1
		if(paramsDict.has_key("useCharge")): useCharge = paramsDict["useCharge"]
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
				TPB.multp(bunch,pole,kl,skew,useCharge)
			TPB.bend4(bunch,theta/2.0)
			TPB.bend3(bunch,theta/2.0)
			TPB.bend2(bunch,length/2.0)
			TPB.bend1(bunch,length,theta)
			return
		if(index == (nParts-1)):
			TPB.bend2(bunch, length)
			TPB.bend3(bunch, theta/2.0)
			TPB.bend4(bunch, theta/2.0)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew,useCharge)
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
			msg = msg + "length =" + str(self.getLength())
			orbitFinalize(msg)

	def track(self, paramsDict):
		"""
		The Ring RF TEAPOT class implementation
		of the AccNodeBunchTracker class track(probe) method.
		"""
		harmArr = self.getParam("harmonics")
		voltArr = self.getParam("voltages")
		phaseArr = self.getParam("phases")
		ring_length = self.getParam("ring_length")
		bunch = paramsDict["bunch"]
		useCharge = 1
		if(paramsDict.has_key("useCharge")): useCharge = paramsDict["useCharge"]	
		length = self.getLength(self.getActivePartIndex())
		for i in range(len(harmArr)):
			#print "debug rl=",ring_length," harm=",harmArr[i]," v=",voltArr[i],
			#print " ph0=",phaseArr[i]," L=",self.getLength()
			TPB.RingRF(bunch,ring_length,harmArr[i],voltArr[i],phaseArr[i],useCharge)
		
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
		self.waveform = None
	
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
		The Kick TEAPOT class implementation of
		the AccNodeBunchTracker class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
		strength = 1.0
		if(self.waveform):
			strength = self.waveform.getStrength()
		kx = strength * self.getParam("kx")/(nParts-1)
		ky = strength * self.getParam("ky")/(nParts-1)
		dE = self.getParam("dE")/(nParts-1)
		bunch = paramsDict["bunch"]
		useCharge = 1
		if(paramsDict.has_key("useCharge")): useCharge = paramsDict["useCharge"]
		if(index == 0):
			TPB.drift(bunch, length)
			TPB.kick(bunch,kx,ky,dE,useCharge)
			return
		if(index > 0 and index < (nParts-1)):
			TPB.drift(bunch, length)
			TPB.kick(bunch,kx,ky,dE,useCharge)
			return
		if(index == (nParts-1)):
			TPB.drift(bunch, length)
		return

	def setWaveform(self, waveform):
		"""
		Sets the time dependent waveform function
		"""
		self.waveform = waveform


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
