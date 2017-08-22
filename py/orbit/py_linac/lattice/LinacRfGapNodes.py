"""
This package is a collection of the RF gap node implementations.
The RF Cavities and gaps in them are different from the ring RF.
"""

import os
import math
import sys

#---- MPI module function and classes
import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

# import from orbit Python utilities
from orbit.utils import orbitFinalize
from orbit.utils import phaseNearTargetPhase, phaseNearTargetPhaseDeg
from orbit.utils import speed_of_light

# import from orbit c++ utilities
from orbit_utils import Polynomial
from orbit_utils import Function

# from LinacAccLattice import Sequence
from LinacAccLatticeLib import Sequence
from LinacAccNodes import Drift, BaseLinacNode

# from linac import the RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF, RfGapThreePointTTF
from linac import BaseRfGap_slow, RfGapTTF_slow, RfGapThreePointTTF_slow

# The abstract RF gap import
from LinacAccNodes import AbstractRF_Gap

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

# Import the linac specific tracking from linac_tracking. This module has
# the following functions duplicated the original TEAPOT functions
# drift - linac drift tracking
# quad1 - linac quad linear part of tracking
# quad2 - linac quad non-linear part of tracking
import linac_tracking

from bunch import Bunch


class BaseRF_Gap(AbstractRF_Gap):
	"""
	The simplest RF gap representation. The only E0*T*L or E0*L and TTFs define 
	all effects of the node. By default the Matrix RF Gap model is used. 
	This model can be replaced later with a more complex RF gap model by using 
	the setCppGapModel(...) method. User should provide the necessary parameters
	for each type of RF gap model.
	MatrixRfGap - E0TL, mode
	BaseRfGap - E0TL, mode
	RfGapTTF - E0L, T,S,Tp,Sp, beta_min, beta_max
	The phase of the first gap in the cavity is defined by the parent cavity instance.
	The relative amplitude is also defined by the parent cavity instance.
	The 'mode' parameter is a shift of the phase between two gaps in PI units.
	"""
	def __init__(self, name = "baserfgap"):
		"""
		Constructor for the simplest RF gap. 
		E0L and E0TL parameters are in GeV. Phases are in radians.
		"""
		AbstractRF_Gap.__init__(self,name)		
		self.addParam("E0TL",0.)
		self.addParam("mode",0.)		
		self.addParam("gap_phase",0.)
		self.addParam("rfCavity", None)
		#----- TTF model params ------
		self.addParam("E0L",0.)	
		self.polyT = Polynomial(0)
		self.polyT.coefficient(0,1.0)
		self.polyS = Polynomial(0)
		self.polyS.coefficient(0,0.0)
		self.polyTp = Polynomial(0)
		self.polyTp.coefficient(0,0.0)
		self.polySp = Polynomial(0)
		self.polySp.coefficient(0,0.0)
		self.addParam("beta_min",0.)
		self.addParam("beta_max",1.)
		#-----------------------------
		self.addParam("rfCavity",None)
		self.addParam("EzFile","no_file")
		self.setType("baserfgap")	
		self.__isFirstGap = False
		#---- by default we use the TTF model 
		#---- which is a Transit-Time-Factor model from Parmila
		#self.cppGapModel = MatrixRfGap()
		#self.cppGapModel = BaseRfGap()
		self.cppGapModel = RfGapTTF()

	def setLinacTracker(self, switch = True):
		"""
		This method will switch RF gap model to slower one where transformations 
		coefficients are calculated for each particle in the bunch.
		"""
		AbstractRF_Gap.setLinacTracker(self,switch)
		if(switch):
			if(isinstance(self.cppGapModel,BaseRfGap)):
				self.cppGapModel = BaseRfGap_slow()
			if(isinstance(self.cppGapModel,RfGapTTF)):
				self.cppGapModel = RfGapTTF_slow()				
		else:
			if(isinstance(self.cppGapModel,BaseRfGap)):
				self.cppGapModel = BaseRfGap()
			if(isinstance(self.cppGapModel,RfGapTTF)):
				self.cppGapModel = RfGapTTF()				

	def setnParts(self, n = 1):
		"""
		Method. Sets the number of body parts of the node. 
		For the RF gap with zero length it will be only 1.
		"""
		BaseLinacNode.setnParts(self,1)


	def setCppGapModel(self, cppGapModel = MatrixRfGap()):
		"""
		This method will set the fast c++ simple model for the RF Gap. 
		By default it is Matrix RF Gap model which is a linear transport matrix.
		"""
		self.cppGapModel = cppGapModel
	
	def initialize(self):
		"""
		The BaseRF_Gap  class implementation
		of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts != 1):
			msg = "The BaseRF_Gap RF gap should have 1 parts!"
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
		self.setLength(0.,0)

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
		
	def getTTF_Polynimials(self):
		"""
		Returns the T,S,Tp,Sp, polynomials in the TTF model.
		"""
		return (self.polyT,self.polyS,self.polyTp,self.polySp)	
		
	def getBetaMinMax(self):
		"""
		Returns beta min and max for TTF model polynomials.
		"""	
		return (self.getParam("beta_min"),self.getParam("beta_max")	)
		
	def setBetaMinMax(self,beta_min,beta_max):
		"""
		Sets beta min and max for TTF model polynomials.
		"""	
		self.setParam("beta_min", beta_min)
		self.setParam("beta_max", beta_max)
		
	def track(self, paramsDict):
		"""
		The simplest RF gap class implementation of
		the AccNode class track(probe) method.
		"""
		bunch = paramsDict["bunch"]		
		syncPart = bunch.getSyncParticle()
		E0TL = self.getParam("E0TL")		
		E0L = self.getParam("E0L")
		modePhase = self.getParam("mode")*math.pi
		rfCavity = self.getRF_Cavity()
		frequency = rfCavity.getFrequency()	
		phase = rfCavity.getPhase() + modePhase
		rf_ampl = rfCavity.getAmp()
		arrival_time = syncPart.time()
		designArrivalTime = rfCavity.getDesignArrivalTime()
		if(self.__isFirstGap):
			if(rfCavity.isDesignSetUp()):
				#print "debug RF =",self.getName(),"  phase=",(phase*180./math.pi - 180.)
				phase = math.fmod(frequency*(arrival_time - designArrivalTime)*2.0*math.pi + phase,2.0*math.pi)
				#print "debug RF =",self.getName(),"  phase=",(phase*180./math.pi - 180.)
			else:
				sequence = self.getSequence()
				accLattice = sequence.getLinacAccLattice()
				msg = "The BaseRF_Gap class. You have to run trackDesign on the LinacAccLattice first to initialize all RF Cavities' phases!"
				msg = msg + os.linesep
				if(accLattice != None):				
					msg = msg + "Lattice =" + accLattice.getName()				
					msg = msg + os.linesep
				if(sequence != None):				
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
			phase = math.fmod(frequency*(arrival_time - designArrivalTime)*2.0*math.pi+phase,2.0*math.pi)	
		#---- rf gap input phase -----
		self.setGapPhase(phase)
		#call rf gap model to track the bunch	
		if(isinstance(self.cppGapModel,MatrixRfGap) or isinstance(self.cppGapModel,BaseRfGap) or isinstance(self.cppGapModel,BaseRfGap_slow)):
			self.cppGapModel.trackBunch(bunch,frequency,E0TL*rf_ampl,phase)
		else:
			self.ttf_track_bunch__(bunch,frequency,E0L*rf_ampl,phase)
		#print "debug delta_time in deg=",frequency*(arrival_time - designArrivalTime)*380.
		#print "debug RF =",self.getName()," E0TL=",E0TL," phase=",(phase*180./math.pi - 180.)," eKin[MeV]=",bunch.getSyncParticle().kinEnergy()*1.0e+3		
		
	def trackDesign(self, paramsDict):
		"""
		The RF First Gap node setups the design time of passage 
		of the bunch through this node.
		"""
		bunch = paramsDict["bunch"]
		eKin_in = bunch.getSyncParticle().kinEnergy()
		E0TL = self.getParam("E0TL")			
		E0L = self.getParam("E0L")			
		rfCavity = self.getRF_Cavity()
		modePhase = self.getParam("mode")*math.pi	
		arrival_time = bunch.getSyncParticle().time()
		frequency = rfCavity.getFrequency()
		phase = rfCavity.getPhase() + modePhase
		if(self.__isFirstGap):
			rfCavity.setDesignArrivalTime(arrival_time)
			rfCavity.setDesignSetUp(True)		
			rfCavity._setDesignPhase(rfCavity.getPhase())
			rfCavity._setDesignAmp(rfCavity.getAmp())
		else:
			first_gap_arr_time = rfCavity.getDesignArrivalTime()
			#print "debug name=",self.getName()," delta_phase=",frequency*(arrival_time - first_gap_arr_time)*360.0," phase=",phase*180/math.pi
			phase = math.fmod(frequency*(arrival_time - first_gap_arr_time)*2.0*math.pi+phase,2.0*math.pi)		
		#print "debug design name=",self.getName()," arr_time=",arrival_time," phase=",phase*180./math.pi," E0TL=",E0TL*1.0e+3," freq=",frequency
		#---- rf gap input phase -----	
		self.setGapPhase(phase)		
		#call rf gap model to track the bunch
		rf_ampl = rfCavity.getDesignAmp()	
		if(isinstance(self.cppGapModel,MatrixRfGap) or isinstance(self.cppGapModel,BaseRfGap) or isinstance(self.cppGapModel,BaseRfGap_slow)):
			self.cppGapModel.trackBunch(bunch,frequency,E0TL*rf_ampl,phase)
		else:
			self.ttf_track_bunch__(bunch,frequency,E0L*rf_ampl,phase)
		eKin_out = bunch.getSyncParticle().kinEnergy()
		#print "debug name=",self.getName()," phase=",(phase*180./math.pi-180.)," Ein=",eKin_in*1000.,"  Eout=",eKin_out*1000.," dE=",(eKin_out-eKin_in)*1000.

	def ttf_track_bunch__(self,bunch,frequency,E0L,phase):
		"""
		Tracks the bunch through the TTF thin gap model. This private method was 
		introduced to to check the beta TTF limits in the polynomial representation
		of T,T',S,and S' functions of the relativistic beta. 
		"""
		beta = bunch.getSyncParticle().beta()
		beta_min = self.getParam("beta_min")	
		beta_max = self.getParam("beta_max")	
		if(beta < beta_min or  beta > beta_max):
			sequence = self.getSequence()
			accLattice = sequence.getLinacAccLattice()
			rfCavity = self.getRF_Cavity()			
			msg = "The Python BaseRF_Gap class. The beta for SyncPart is not in the range [min:max]!"
			msg = msg + os.linesep
			if(accLattice != None):				
				msg = msg + "Lattice =" + accLattice.getName()				
				msg = msg + os.linesep
			if(sequence != None):				
				msg = msg + "Sequence =" + sequence.getName()				
				msg = msg + os.linesep
			msg = msg + "RF Cavity =" + rfCavity.getName()				
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "beta=" + str(beta)
			msg = msg + os.linesep				
			msg = msg + "beta min=" + str(beta_min)
			msg = msg + os.linesep				
			msg = msg + "beta max=" + str(beta_max)
			msg = msg + os.linesep				
			orbitFinalize(msg)				
		self.cppGapModel.trackBunch(bunch,frequency,E0L,phase,self.polyT,self.polyS,self.polyTp,self.polySp)

#-----------------------------------------------------------------------
#
# This part of the package is for classes related to the axis RF fields
#
#------------------------------------------------------------------------

class RF_AxisFieldsStore:
	
	"""
	The dictionary with the axis field Functions 
	with the input file names as keys.
	This is a collection of the static methods.
	"""
	
	#---- static_axis_field_dict[file_name] = Function
	static_axis_field_dict = {}
	
	def __init__(self):
		pass
	
	@classmethod
	def addAxisFieldsForAccSeq(cls,accLattice,accSeqNamesList,dir_location = ""):
		"""
		This method add to the store the axis RF fields of all RF gap nodes 
		(BaseRF_Gap class instance with "EzFile" parameter) from the set of accSeqences.
		The dir_location string variable will be added to the rf_gap.getParam("EzFile") to get
		the file names.
		"""		
		for accNamesSeq in accSeqList:
			accSeq = accLattice.getSequence(accNamesSeq)
			cavs = accSeq.getRF_Cavities()
			for cav in cavs:
				rf_gaps = cav.getRF_GapNodes()
				for rf_gap in rf_gaps:
					cls.addAxisField(rf_gap.getParam("EzFile"),dir_location)
	
	@classmethod
	def addAxisField(cls,fl_name,dir_location = ""):
		"""
		This method add to the store the axis RF field for the RF gap node. 
		The dir_location string variable will be added to the fl_name to get
		the file name.
		Returns the axis RF field function.
		"""
		if(cls.static_axis_field_dict.has_key(fl_name)): 
			return cls.static_axis_field_dict[fl_name]
		comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
		data_type = mpi_datatype.MPI_DOUBLE
		rank = orbit_mpi.MPI_Comm_rank(comm)
		main_rank = 0
		x_arr = []
		y_arr = []
		if(rank == 0):
			fl_in = open(dir_location + fl_name,"r")
			lns = fl_in.readlines()
			fl_in.close()
			for ln in lns:
				res_arr = ln.split()
				if(len(res_arr) == 2):
					x = float(res_arr[0])
					y = float(res_arr[1])
					x_arr.append(x)		
					y_arr.append(y)	
		x_arr = orbit_mpi.MPI_Bcast(x_arr,data_type,main_rank,comm)
		y_arr = orbit_mpi.MPI_Bcast(y_arr,data_type,main_rank,comm)
		function = Function()
		for ind in range(len(x_arr)):
			function.add(x_arr[ind],y_arr[ind])
		#---- setting the const step (if function will allow it) 
		#---- will speed up function calculation later
		function.setConstStep(1)
		cls.static_axis_field_dict[fl_name] = function
		return function
		
	@classmethod
	def getAxisFieldFunction(cls,fl_name):
		"""
		This method returns the Function with the RF axis field for a particular
		name of the file. If store does not have this function it will return
		the None object.
		"""
		if(cls.static_axis_field_dict.has_key(fl_name)):
			return cls.static_axis_field_dict[fl_name]
		else:
			return None

	@classmethod
	def getSize(cls):
		"""
		This method returns the number of Functions with the RF axis fields in this
		store.
		"""
		return len(cls.static_axis_field_dict.keys())
		
		
class AxisFieldRF_Gap(AbstractRF_Gap):
	"""
	The RF gap representation that uses the RF axis field. User have to provide the 
	input file with this field. This function should be normalized to the integral of 1.
	The absolute value of the field will be calculated as cavAmp*E0L*Field(z).
	The three point tracker RfGapThreePointTTF will be used to track the Bunch instance.
	The  longitudinal step during the tracking z_step should be defined externally. The 
	default value is 1 cm. The minimal and maximal longitudinal coordinates z_min
	and z_max could be used directly from the axis field file or can be corrected 
	externally to avoid overlapping of electric fields from neighboring gaps.
	The instance of this class has the reference to the BaseRF_Gap instance and uses 
	it as a source of information.
	"""
	
	#---- static test bunch for the design phase calculation 
	static_test_bunch = Bunch()
	
	def __init__(self, baserf_gap):
		"""
		Constructor for the axis field RF gap. 
		E0L parameter is in GeV. Phases are in radians.
		"""
		AbstractRF_Gap.__init__(self,baserf_gap.getName())
		self.setAsFirstRFGap(baserf_gap.isFirstRFGap())
		self.baserf_gap = baserf_gap
		self.setType("axis_field_rfgap")			
		self.addParam("E0TL",self.baserf_gap.getParam("E0TL"))
		self.addParam("mode",self.baserf_gap.getParam("mode"))
		self.addParam("gap_phase",self.baserf_gap.getParam("gap_phase"))
		self.addParam("rfCavity",self.baserf_gap.getParam("rfCavity"))
		self.addParam("E0L",self.baserf_gap.getParam("E0L"))
		self.addParam("EzFile",self.baserf_gap.getParam("EzFile"))
		self.setPosition(self.baserf_gap.getPosition())
		#---- aperture parameters
		if(baserf_gap.hasParam("aperture") and baserf_gap.hasParam("aprt_type")):
			self.addParam("aperture",baserf_gap.getParam("aperture"))
			self.addParam("aprt_type",baserf_gap.getParam("aprt_type"))
		#---- axis field related parameters
		self.axis_field_func = None
		self.z_step = 0.01
		self.z_min = 0.
		self.z_max = 0.
		self.z_tolerance = 0.000001    # in meters
		self.phase_tolerance = 0.001   # in degrees
		#---- gap_phase_vs_z_arr keeps [pos,phase] pairs after the tracking
		self.gap_phase_vs_z_arr = []
		#---- The position of the particle during the run. 
		#---- It is used for the path length accounting.
		self.part_pos = 0.
		#---- The RF gap model - three points model
		self.cppGapModel = RfGapThreePointTTF()

	def setLinacTracker(self, switch = True):
		"""
		This method will switch RF gap model to slower one where transformations 
		coefficients are calculated for each particle in the bunch.
		"""
		AbstractRF_Gap.setLinacTracker(self,switch)
		if(switch):
			self.cppGapModel = RfGapThreePointTTF_slow()			
		else:
			self.cppGapModel = RfGapThreePointTTF()

	def readAxisFieldFile(self,dir_location = "", file_name = "", z_step = 0.01):
		"""
		Method. Reads the axis field from the file. User have to call this method.
		There is no other source of information about the axis field.
		"""
		if(file_name == ""):
			self.axis_field_func = RF_AxisFieldsStore.addAxisField(self.baserf_gap.getParam("EzFile"),dir_location)
		else:
			self.axis_field_func = RF_AxisFieldsStore.addAxisField(file_name,dir_location)
		z_min = self.axis_field_func.getMinX()
		z_max = self.axis_field_func.getMaxX()
		self.z_step = z_step
		self.setZ_Min_Max(z_min,z_max)
		
	def getAxisFieldFunction(self):
		"""
		It returns the axis field function.
		"""
		return self.axis_field_func
		
	def setAxisFieldFunction(self,axis_field_func):
		"""
		It sets the axis field function.
		"""
		self.axis_field_func = axis_field_func	
		
	def getZ_Step(self):
		"""
		Returns the longitudinal step during the tracking.
		"""
		return self.z_step
		
	def setZ_Step(self,z_step):
		if(self.axis_field_func == None):
			msg = "Class AxisFieldRF_Gap: You have to get the axis field from a file first!"
			msg = msg + os.linesep
			msg = "Call readAxisFieldFile(dir_location,file_name) method first!"
			orbitFinalize(msg)				
		length = self.getLength()
		nParts = int(length*1.0000001/z_step)
		if(nParts < 1): nParts = 1
		self.z_step = length/nParts
		#---- this will set the even distribution of the lengths between parts
		self.setnParts(nParts)

	def getZ_Min_Max(self):
		"""
		Returns the tuple (z_min,z_max) with the limits of the axis field.
		These parameters define the length of the node. The center of the node
		is at 0.
		"""
		return (self.z_min,self.z_max)

	def setZ_Min_Max(self,z_min,z_max):
		"""
		Sets the actual longitudinal sizes of the node. It is used for small correction
		of the length to avoid fields overlapping from neighbouring gaps.
		"""
		self.z_min = z_min
		self.z_max = z_max
		length = self.z_max - self.z_min
		self.setLength(length)
		self.setZ_Step(self.z_step)

	def getEzFiled(self,z):
		"""
		Returns the Ez field on the axis of the RF gap in V/m. 
		"""
		rfCavity = self.getRF_Cavity()
		E0L = 1.0e+9*self.getParam("E0L")
		rf_ampl = rfCavity.getAmp()
		Ez = E0L*rf_ampl*self.axis_field_func.getY(z)	
		return Ez

	def getRF_Cavity(self):
		"""
		Returns the parent RF Cavity.
		"""
		return self.getParam("rfCavity")
				
	def track(self, paramsDict):
		"""
		The AxisFieldRF_Gap class implementation of
		the AccNode class track(probe) method.
		User have to track the design bunch first to setup all gaps arrival time. 
		"""
		rfCavity = self.getRF_Cavity()
		if(not rfCavity.isDesignSetUp()):
			sequence = self.getSequence()
			accLattice = sequence.getLinacAccLattice()
			msg  = "The AxisFieldRF_Gap class. "
			msg += "You have to run trackDesign on the LinacAccLattice"
			msg += "first to initialize all RF Cavities' phases!"
			msg += os.linesep
			if(accLattice != None):				
				msg = msg + "Lattice =" + accLattice.getName()				
				msg = msg + os.linesep
			if(sequence != None):				
				msg = msg + "Sequence =" + sequence.getName()				
				msg = msg + os.linesep
			msg = msg + "RF Cavity =" + rfCavity.getName()				
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			orbitFinalize(msg)
		#-----------------------------------------
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		part_length = self.getLength(index)		
		bunch = paramsDict["bunch"]		
		syncPart = bunch.getSyncParticle()	
		eKin_in = syncPart.kinEnergy()
		E0L = 1.0e+9*self.getParam("E0L")
		modePhase = self.getParam("mode")*math.pi
		frequency = rfCavity.getFrequency()	
		rf_ampl = rfCavity.getAmp()
		arrival_time = syncPart.time()
		designArrivalTime = rfCavity.getDesignArrivalTime()
		phase_shift = rfCavity.getPhase() - rfCavity.getDesignPhase()
		phase = rfCavity.getFirstGapEtnrancePhase() + phase_shift
		#----------------------------------------
		phase = math.fmod(frequency*(arrival_time - designArrivalTime)*2.0*math.pi + phase,2.0*math.pi)		
		if(index == 0):
			self.part_pos = self.z_min 
			self.gap_phase_vs_z_arr = [[self.part_pos,phase],]
		zm = self.part_pos
		z0 = zm + part_length/2
		zp = z0 + part_length/2
		Em = E0L*rf_ampl*self.axis_field_func.getY(zm)
		E0 = E0L*rf_ampl*self.axis_field_func.getY(z0)
		Ep = E0L*rf_ampl*self.axis_field_func.getY(zp)			
		#---- advance the particle position
		self.tracking_module.drift(bunch,part_length/2)
		self.part_pos += part_length/2	
		#call rf gap model to track the bunch
		time_middle_gap = syncPart.time() - arrival_time
		delta_phase = math.fmod(2*math.pi*time_middle_gap*frequency,2.0*math.pi)
		self.gap_phase_vs_z_arr.append([self.part_pos,phase+delta_phase])
		#---- this part is the debugging ---START---
		#eKin_out = syncPart.kinEnergy()
		#s  = "debug pos[mm]= %7.2f "%(self.part_pos*1000.)
		#s += " ekin= %9.6f"%(syncPart.kinEnergy()*1000.)
		#s += " phase = %9.2f "%(phaseNearTargetPhaseDeg((phase+delta_phase)*180./math.pi,0.))
		#s += " dE= %9.6f "%((eKin_out-eKin_in)*1000.)		
		#print s
		#---- this part is the debugging ---STOP---
		self.cppGapModel.trackBunch(bunch,part_length/2,Em,E0,Ep,frequency,phase+delta_phase+modePhase)
		self.tracking_module.drift(bunch,part_length/2)
		#---- advance the particle position
		self.part_pos += part_length/2
		time_middle_gap = syncPart.time() - arrival_time
		delta_phase = math.fmod(2*math.pi*time_middle_gap*frequency,2.0*math.pi)
		self.gap_phase_vs_z_arr.append([self.part_pos,phase+delta_phase])
		#---- this part is the debugging ---START---
		#eKin_out = syncPart.kinEnergy()
		#s  = "debug pos[mm]= %7.2f "%(self.part_pos*1000.)
		#s += " ekin= %9.6f"%(syncPart.kinEnergy()*1000.)
		#s += " phase = %9.2f "%(phaseNearTargetPhaseDeg((phase+delta_phase)*180./math.pi,0.))
		#s += " dE= %9.6f "%((eKin_out-eKin_in)*1000.)
		#print s
		#---- this part is the debugging ---STOP---
		#---- Calculate the phase at the center
		if(index == (nParts - 1)):
			pos_old = self.gap_phase_vs_z_arr[0][0]			
			phase_gap = self.gap_phase_vs_z_arr[0][1]
			ind_min = -1
			for ind in range(1,len(self.gap_phase_vs_z_arr)):
				[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
				if(math.fabs(pos) >= math.fabs(pos_old)):
					ind_min = ind -1
					phase_gap = self.gap_phase_vs_z_arr[ind_min][1]
					phase_gap = phaseNearTargetPhase(phase_gap,0.)
					self.gap_phase_vs_z_arr[ind_min][1] = phase_gap
					break
				pos_old = pos
			self.setGapPhase(phase_gap)
			#---- wrap all gap part's phases around the central one
			if(ind_min > 0):
				for ind in range(ind_min-1,-1,-1):
					[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
					[pos,phase_gap1] = self.gap_phase_vs_z_arr[ind+1]
					self.gap_phase_vs_z_arr[ind][1] = phaseNearTargetPhase(phase_gap,phase_gap1)
				for ind in range(ind_min+1,len(self.gap_phase_vs_z_arr)):
					[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
					[pos,phase_gap1] = self.gap_phase_vs_z_arr[ind-1]				
					self.gap_phase_vs_z_arr[ind][1] = phaseNearTargetPhase(phase_gap,phase_gap1)

	def trackDesign(self, paramsDict):
		"""
		The method is tracking the design synchronous particle through the RF Gap.
		If the gap is a first gap in the cavity we put the arrival time as 
		a cavity parameter. The pair of the cavity design phase and this arrival time 
		at the first gap are used during the real bunch tracking.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		part_length = self.getLength(index)
		bunch = paramsDict["bunch"]
		syncPart = bunch.getSyncParticle()
		eKin_in = syncPart.kinEnergy()
		#---- parameter E0L is in GeV, but cppGapModel = RfGapThreePointTTF() uses fields in V/m
		E0L = 1.0e+9*self.getParam("E0L")
		modePhase = self.getParam("mode")*math.pi
		rfCavity = self.getRF_Cavity()
		rf_ampl = rfCavity.getDesignAmp()
		arrival_time = syncPart.time()
		frequency = rfCavity.getFrequency()
		phase = rfCavity.getFirstGapEtnrancePhase()
		#---- calculate the entance phase
		if(self.isFirstRFGap() and index == 0):
			rfCavity.setDesignArrivalTime(arrival_time)
			phase = self.__calculate_first_part_phase(bunch)
			rfCavity.setFirstGapEtnrancePhase(phase)
			rfCavity.setFirstGapEtnranceDesignPhase(phase)
			rfCavity.setDesignSetUp(True)		
			rfCavity._setDesignPhase(rfCavity.getPhase())
			rfCavity._setDesignAmp(rfCavity.getAmp())
			#print "debug firs gap first part phase=",phase*180./math.pi," arr time=",arrival_time
		else:
			first_gap_arr_time = rfCavity.getDesignArrivalTime()
			#print "debug name=",self.getName()," delta_phase=",frequency*(arrival_time - first_gap_arr_time)*360.0," phase=",phase*180/math.pi
			phase = math.fmod(frequency*(arrival_time - first_gap_arr_time)*2.0*math.pi+phase,2.0*math.pi)		
		if(index == 0):
			self.part_pos = self.z_min 
			self.gap_phase_vs_z_arr = [[self.part_pos,phase],]
		#print "debug design name=",self.getName()," index=",index," pos=",self.part_pos," arr_time=",arrival_time," phase=",phase*180./math.pi," freq=",frequency
		zm = self.part_pos
		z0 = zm + part_length/2
		zp = z0 + part_length/2
		Em = E0L*rf_ampl*self.axis_field_func.getY(zm)
		E0 = E0L*rf_ampl*self.axis_field_func.getY(z0)
		Ep = E0L*rf_ampl*self.axis_field_func.getY(zp)
		#---- advance the particle position
		self.tracking_module.drift(bunch,part_length/2)
		self.part_pos += part_length/2	
		#call rf gap model to track the bunch
		time_middle_gap = syncPart.time() - arrival_time
		delta_phase = math.fmod(2*math.pi*time_middle_gap*frequency,2.0*math.pi)
		self.gap_phase_vs_z_arr.append([self.part_pos,phase+delta_phase])
		#---- this part is the debugging ---START---
		#eKin_out = syncPart.kinEnergy()
		#s  = "debug pos[mm]= %7.2f "%(self.part_pos*1000.)
		#s += " ekin= %9.6f"%(syncPart.kinEnergy()*1000.)
		#s += " phase = %9.2f "%(phaseNearTargetPhaseDeg((phase+delta_phase)*180./math.pi,0.))
		#s += " dE= %9.6f "%((eKin_out-eKin_in)*1000.)		
		#print s
		#---- this part is the debugging ---STOP---
		self.cppGapModel.trackBunch(bunch,part_length/2,Em,E0,Ep,frequency,phase+delta_phase+modePhase)
		self.tracking_module.drift(bunch,part_length/2)
		#---- advance the particle position
		self.part_pos += part_length/2
		time_middle_gap = syncPart.time() - arrival_time
		delta_phase = math.fmod(2*math.pi*time_middle_gap*frequency,2.0*math.pi)
		self.gap_phase_vs_z_arr.append([self.part_pos,phase+delta_phase])
		#---- this part is the debugging ---START---
		#eKin_out = syncPart.kinEnergy()
		#s  = "debug pos[mm]= %7.2f "%(self.part_pos*1000.)
		#s += " ekin= %9.6f"%(syncPart.kinEnergy()*1000.)
		#s += " phase = %9.2f "%(phaseNearTargetPhaseDeg((phase+delta_phase)*180./math.pi,0.))
		#s += " dE= %9.6f "%((eKin_out-eKin_in)*1000.)
		#print s
		#---- this part is the debugging ---STOP---
		#---- Calculate the phase at the center
		if(index == (nParts - 1)):
			pos_old = self.gap_phase_vs_z_arr[0][0]			
			phase_gap = self.gap_phase_vs_z_arr[0][1]
			ind_min = -1
			for ind in range(1,len(self.gap_phase_vs_z_arr)):
				[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
				if(math.fabs(pos) >= math.fabs(pos_old)):
					ind_min = ind -1
					phase_gap = self.gap_phase_vs_z_arr[ind_min][1]
					phase_gap = phaseNearTargetPhase(phase_gap,0.)
					self.gap_phase_vs_z_arr[ind_min][1] = phase_gap
					break
				pos_old = pos
			self.setGapPhase(phase_gap)
			#---- wrap all gap part's phases around the central one
			if(ind_min > 0):
				for ind in range(ind_min-1,-1,-1):
					[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
					[pos,phase_gap1] = self.gap_phase_vs_z_arr[ind+1]
					self.gap_phase_vs_z_arr[ind][1] = phaseNearTargetPhase(phase_gap,phase_gap1)
				for ind in range(ind_min+1,len(self.gap_phase_vs_z_arr)):
					[pos,phase_gap] = self.gap_phase_vs_z_arr[ind]
					[pos,phase_gap1] = self.gap_phase_vs_z_arr[ind-1]				
					self.gap_phase_vs_z_arr[ind][1] = phaseNearTargetPhase(phase_gap,phase_gap1)

	def calculate_first_part_phase(self,bunch_in):
		"""
		The privat method should be exposed to the AxisField_and_Quad_RF_Gap class
		"""
		phase_start = self.__calculate_first_part_phase(bunch_in)
		return phase_start

	def __calculate_first_part_phase(self,bunch_in):
		rfCavity = self.getRF_Cavity()
		#---- the design phase at the center of the RF gap 
		#---- (this is from a thin gap approach)
		frequency = rfCavity.getFrequency()
		modePhase = self.getParam("mode")*math.pi
		phase_cavity = rfCavity.getPhase()
		#---- parameter E0L is in GeV, but cppGapModel = RfGapThreePointTTF() uses fields in V/m
		E0L_local = 1.0e+9*rfCavity.getAmp()*self.getParam("E0L")		
		cav_ampl = rfCavity.getAmp()
		#---- we have to find the phase_start 
		#---- which is the phase at the distance z_min before the gap center
		#---- z_min by defenition is negative
		bunch = AxisFieldRF_Gap.static_test_bunch
		bunch_in.copyEmptyBunchTo(bunch)
		syncPart = bunch.getSyncParticle()
		syncPart.time(0.)
		eKin_init = syncPart.kinEnergy()		
		#print "debug eKin[MeV]= %9.5f"%(syncPart.kinEnergy()*1000.)
		beta = syncPart.beta()
		phase_adv = 2.0*math.pi*frequency*math.fabs(self.z_min)/(beta*speed_of_light)
		#print "debug phase diff at start=",phase_adv*180./math.pi
		phase_start = phaseNearTargetPhase(phase_cavity - phase_adv,0.)
		#print "debug phase at start=",phase_start*180./math.pi
		phase_cavity_new = phase_cavity + 10*self.phase_tolerance
		while(math.fabs(phase_cavity_new-phase_cavity) > self.phase_tolerance*math.pi/180.):
			bunch_in.copyEmptyBunchTo(bunch)
			syncPart.time(0.)
			syncPart.kinEnergy(eKin_init)
			z_old = self.z_min
			z = self.z_min + self.z_step
			while(z < 0.):
				if((z+ self.z_step) > 0.): 
					z = 0.
					if(math.fabs(z - z_old) < self.z_tolerance):
						break
				half_step = (z - z_old)/2
				zm = z_old
				z0 = zm + half_step
				zp = z0 + half_step
				self.tracking_module.drift(bunch,half_step)
				time_gap = syncPart.time()
				delta_phase = 2*math.pi*time_gap*frequency 
				Em = E0L_local*self.axis_field_func.getY(zm)
				E0 = E0L_local*self.axis_field_func.getY(z0)
				Ep = E0L_local*self.axis_field_func.getY(zp)
				#s  = "debug z[mm]= %7.2f "%(z0*1000.)
				#s += " ekin= %9.5f"%(syncPart.kinEnergy()*1000.)
				#s += " phase = %9.2f "%(phaseNearTargetPhaseDeg((phase_start+delta_phase+modePhase)*180./math.pi,0.))
				#print s	
				self.cppGapModel.trackBunch(bunch,half_step,Em,E0,Ep,frequency,phase_start+delta_phase+modePhase)
				self.tracking_module.drift(bunch,half_step)
				#time_gap = syncPart.time()
				#delta_phase = 2*math.pi*time_gap*frequency 				
				#s  = "debug z[mm]= %7.2f "%(zp*1000.)
				#s += " ekin= %9.5f"%(syncPart.kinEnergy()*1000.)
				#s += " phase = %9.2f "%(phaseNearTargetPhaseDeg((phase_start+delta_phase+modePhase)*180./math.pi,0.))
				#print s			
				z_old = z
				z =  z_old + self.z_step
			time_gap = syncPart.time()
			delta_phase =2*math.pi*time_gap*frequency 
			phase_cavity_new = phaseNearTargetPhase(phase_start+delta_phase,0.)
			#s  = " phase_diff = %8.4f "%(delta_phase*180./math.pi)
			#s += " phase_cavity = %8.4f "%(phase_cavity*180./math.pi)
			#s += " new = %8.4f "%(phase_cavity_new *180./math.pi)
			#s += " phase_start = %8.4f "%(phase_start*180./math.pi)
			#s += " eKin[MeV]= %9.5f "%(syncPart.kinEnergy()*1000.)
			#s += " dE[MeV]= %9.6f "%(syncPart.kinEnergy()*1000. - 2.5)
			#print "debug "+s
			phase_start -= 0.8*(phase_cavity_new - phase_cavity)
		#---- undo the last change in the while loop
		phase_start += 0.8*(phase_cavity_new - phase_cavity)
		#print "debug phase_start=",phase_start*180./math.pi
		return phase_start

