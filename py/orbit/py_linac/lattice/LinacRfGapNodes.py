"""
This package is a collection of the RF gap node implementations.
The RF Cavities and gaps in them are different from the ring RF.
"""

import os
import math

# import from orbit utilities
from orbit.utils import orbitFinalize
from orbit_utils import Polynomial

# from LinacAccLattice import Sequence
from LinacAccLatticeLib import Sequence

# from linac import the RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF

# The abstract RF gap import
from LinacAccNodes import AbstractRF_Gap

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
		self.cppGapModel = MatrixRfGap()
		self.cppGapModel = BaseRfGap()
		self.cppGapModel = RfGapTTF()

	def setnParts(self, n = 1):
		"""
		Method. Sets the number of body parts of the node. For the RF gap it will be only 1.
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
		The Ring RF TEAPOT class implementation
		of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts != 1):
			msg = "The simple Rf gap should have 1 parts!"
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
		if(isinstance(self.cppGapModel,MatrixRfGap) or isinstance(self.cppGapModel,BaseRfGap)):
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
		if(isinstance(self.cppGapModel,MatrixRfGap) or isinstance(self.cppGapModel,BaseRfGap)):
			self.cppGapModel.trackBunch(bunch,frequency,E0TL*rf_ampl,phase)
		else:
			self.ttf_track_bunch__(bunch,frequency,E0L*rf_ampl,phase)
		eKin_out = bunch.getSyncParticle().kinEnergy()
		#print "debug name=",self.getName()," phase=",(phase*180./math.pi-180.)," Ein=",eKin_in*1000.,"  Eout=",eKin_out*1000.," dE=",(eKin_out-eKin_in)*1000.

	def ttf_track_bunch__(self,bunch,frequency,E0L,phase):
		"""
		Tracks the bunch through the TTF thin gap model. 
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

