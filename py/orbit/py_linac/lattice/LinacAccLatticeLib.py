"""
The general Linac Accelerator lattice. It is a subclass of the Acclattice. 
It tracks the bunch through the linac accelerator nodes. In addition to the 
usual accelerator lattice it has the sequences and RF cavities. The sequences and 
RF cavities are containers for accelerator nodes (seqencies) and RF gaps. 
The Sequence class is used to keep iformation about positions of elements that are inside. 
The Cavity class keeps the refernce to RF gaps and a value of the cavity amplitude.
"""

import os
import math

# import the utilities
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject, phaseNearTargetPhase

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer

# import Sequence and RF_Cavity
from LinacAccNodes import Quad, AbstractRF_Gap, MarkerLinacNode

# import orbit Bunch
from bunch import Bunch

class LinacAccLattice(AccLattice):
	"""
	The subclass of the AccLattice class. In the beginning the lattcie is empty.
	"""
	def __init__(self, name = None):
		AccLattice.__init__(self,name)
		self.__rfCavities = []
		self.__sequences = []
		
	def initialize(self):
		"""
		Method. Initializes the linac lattice, child node structures, and calculates 
		the one turn matrix.
		"""
		AccLattice.initialize(self)
		#-----add Sequences that are referenced from the AccNodes 
		self.__sequences = []
		for node in self.getNodes():
			seq = node.getSequence()
			if(not seq in self.__sequences):
				self.__sequences.append(seq)	
		#------define sequences' position and length (position of the beginning of the sequence)
		seqs = self.getSequences()
		for seq in seqs:
			nodes = seq.getNodes()
			if(len(nodes) == 0): continue
			node_start = seq.getNodes()[0]
			node_last = seq.getNodes()[len(nodes)-1]
			(pos_start,pos_tmp) = self.getNodePositionsDict()[node_start]
			(pos_tmp,pos_stop) = self.getNodePositionsDict()[node_last]
			seq.setPosition(pos_start)
			seq.setLength(pos_stop-pos_start)
		#------add RF cavities that are in referenced from the AccNodes (RF gaps)
		self.__rfCavities = []
		for node in self.getNodes():
			if(isinstance(node,AbstractRF_Gap)):
				rf_cavity = node.getRF_Cavity()
				if(rf_cavity == None):
					msg = "LinacAccLattice.initialize() - Problem with RF gap!"
					msg = msg + os.linesep
					msg = msg + "RF gap="+node.getName()+" does not have the parent RF cavity!"
					msg = msg + os.linesep
					msg = msg + "Stop."
					msg = msg + os.linesep
					orbitFinalize(msg)					
				if(not rf_cavity in self.__rfCavities):
					self.__rfCavities.append(rf_cavity)
					
	def getSubLattice(self, index_start = -1, index_stop = -1):
		"""
		It returns the new LinacAccLattice with children with indexes 
		between index_start and index_stop inclusive. 
		UNFORTUNATELY: At this moment it is not possible because the lattice has 
		RF gap with unique references to RF cavities and Sequences which are
		bidirectional.
		"""
		msg = "The getSubLattice method of LinacAccLattice!"
		msg = msg + os.linesep
		msg = msg + "That is not possible to get a sub-lattice for LinacAccLattice."
		msg = msg + os.linesep
		msg = msg + "Use the trackBunch method with start and end indexes."
		msg = msg + os.linesep
		msg = msg + "Stop."
		msg = msg + os.linesep
		orbitFinalize(msg)			
		return self._getSubLattice(LinacAccLattice(),index_start,index_stop)

	def trackBunch(self, bunch, paramsDict = None, actionContainer = None, index_start = -1, index_stop = -1):
		"""
		It tracks the bunch through the lattice.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Bunch Tracking")
		if(paramsDict == None): paramsDict = {}			
		paramsDict["bunch"] = bunch
		
		def track(paramsDict):
			node = paramsDict["node"]
			node.track(paramsDict)
			
		actionContainer.addAction(track, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict,index_start,index_stop)
		actionContainer.removeAction(track, AccActionsContainer.BODY)

	def trackDesignBunch(self, bunch_in, paramsDict = None, actionContainer = None, index_start = -1, index_stop = -1):
		"""
		This will track the design bunch through the linac and set up RF Cavities times of
		arrivals.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Design Bunch Tracking")
		if(paramsDict == None): paramsDict = {}		
		bunch = Bunch()
		bunch_in.copyEmptyBunchTo(bunch)
		bunch.getSyncParticle().time(0.)	
		paramsDict["bunch"] = bunch
		
		def trackDesign(localParamsDict):
			node = localParamsDict["node"]
			node.trackDesign(localParamsDict)
			
		actionContainer.addAction(trackDesign, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict,index_start,index_stop)
		actionContainer.removeAction(trackDesign, AccActionsContainer.BODY)
		return bunch

	def getRF_Cavity(self,name):
		""" Returns the cavity instance according to the name """
		for cav in self.__rfCavities:
			if(name == cav.getName()):
				return cav
		return None

	def getRF_Cavities(self):
		""" Returns the array with RF cavities. """
		return self.__rfCavities		

	def addSequence(self,seq):
		if(isinstance(seq, Sequence) == True):
			self.__sequences.append(seq)
		else:
			msg = "The LinacAccLattice, method addSequence(seq)!"
			msg = msg + os.linesep
			msg = msg + "seq is not a subclass of Sequence."
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)	

	def getSequence(self,name):
		""" Returns the sequence instance according to the name """
		for seq in self.__sequences:
			if(name == seq.getName()):
				return seq
		return None
		
	def getSequences(self):
		""" Returns the array with sequences. """
		return self.__sequences
		
	def getQuads(self, seq = None):
		""" Returns the list of all quads or just quads belong to a particular sequence. """ 
		quads = []
		for node in self.getNodes():
			if(isinstance(node,Quad)):
				if(seq == None):
					quads.append(node)
				else:
					if(node.getSequence() == seq):
						quads.append(node)
		return quads
		
	def getNodesOfClass(self, Class, seq = None):
		""" 
		Returns the list of all nodes which are instances of Class
		or these nodes which are also belong to a particular sequence. 
		""" 
		nodes = []
		for node in self.getNodes():
			if(isinstance(node,Class)):
				if(seq == None):
					nodes.append(node)
				else:
					if(node.getSequence() == seq):
						nodes.append(node)
		return nodes
		
	def getNodesOfClasses(self, Classes, seq = None):
		""" 
		Returns the list of all nodes which are instances of any Class in Classes array
		or these nodes which are also belong to a particular sequence. 
		""" 
		nodes = []
		for node in self.getNodes():
			info = False
			for Class in Classes:
				if(isinstance(node,Class)):
					info = True
			if(info):
				if(seq == None):
					nodes.append(node)
				else:
					if(node.getSequence() == seq):
						nodes.append(node)
		return nodes	

	def getRF_Gaps(self, rf_cav = None):
		""" Returns the list of all RF gaps or just gaps belong to a particular RF cavity. """
		gaps = []
		for node in self.getNodes():
			if(isinstance(node,AbstractRF_Gap)):
				if(rf_cav == None):
					gaps.append(node)
				else:
					if(node.getRF_Cavity() == rf_cav):
						gaps.append(node)
		return gaps
	
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
		#---- firstGapEntrancePhase is used by gaps instead of phase in the case
		#---- of non-zero length gaps
		self.__firstGapEntrancePhase = 0.
		self.__firstGapEntranceDesignPhase = 0.
		self.addParam("frequency",0.)
		self.addParam("phase",0.)
		self.addParam("amp",1.)		
		self.addParam("designPhase",0.)
		self.addParam("designAmp",1.)		
		self.addParam("designArrivalTime",0.)
		self.addParam("isDesignSetUp",False)
		self.addParam("pos",0.)
		
	def setDesignSetUp(self,designOnOf):
		""" Sets the design set up information (True,False). """
		self.setParam("isDesignSetUp",designOnOf)	

	def isDesignSetUp(self):
		""" Returns the design set up information (True,False). """
		return self.getParam("isDesignSetUp")	
		
	def setDesignArrivalTime(self,time):
		""" Sets the design arrival time for the first RF gap. """
		self.setParam("designArrivalTime",time)
		
	def getDesignArrivalTime(self):
		""" Returns the design arrival time for the first RF gap. """
		return self.getParam("designArrivalTime")
		
	def _setDesignPhase(self,phase):
		""" Sets the design phase for the first RF gap. This method is called from the design tracking. """
		self.setParam("designPhase",phase)
		
	def getDesignPhase(self):
		""" Returns the design phase for the first RF gap. """
		return self.getParam("designPhase")

	def _setDesignAmp(self,Amp):
		""" Sets the design Amp for the RF cavity. This method is called from the design tracking."""
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

	def setFirstGapEtnrancePhase(self,phase):
		""" Sets the phase at the first gap entrance if Length_of_gap > 0. """
		self.__firstGapEntrancePhase = phase
		
	def getFirstGapEtnrancePhase(self):
		""" Returns the phase at the first gap entrance if Length_of_gap > 0. """
		return self.__firstGapEntrancePhase

	def setFirstGapEtnranceDesignPhase(self,phase):
		""" Sets the design phase at the first gap entrance if Length_of_gap > 0. """
		self.__firstGapEntranceDesignPhase = phase
		
	def getFirstGapEtnranceDesignPhase(self):
		""" Returns the design phase at the first gap entrance if Length_of_gap > 0. """
		return self.__firstGapEntranceDesignPhase

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
		
	def setPosition(self,pos):
		""" Sets the position of the RF cavity in the sequence. """
		self.setParam("pos",pos)
		
	def getPosition(self):
		""" Returns the position of the RF cavity in the sequence. """
		return self.getParam("pos")		
		
	def addRF_GapNode(self,rfGap):
		""" Adds the rf gap to the cavity."""
		self.__rfGaps.append(rfGap)
		rfGap.setRF_Cavity(self)
		if(len(self.__rfGaps) == 1):
			rfGap.setAsFirstRFGap(True)
		else:
			rfGap.setAsFirstRFGap(False)
		
	def getAvgGapPhase(self):
		""" Returns average phase for all RF gaps in the cavity """ 
		avg_phase = 0.
		phase = 0.
		phase_arr = []
		if(len(self.__rfGaps) > 0): 
				phase = phaseNearTargetPhase(self.__rfGaps[0].getGapPhase(),0.)				
		for rfGap in self.__rfGaps:
			phase_new = phaseNearTargetPhase(rfGap.getGapPhase(),phase)
			phase_arr.append(phase_new*180./math.pi)
			avg_phase += phase_new
			phase = phase_new
		if(len(self.__rfGaps) > 0): avg_phase /= len(self.__rfGaps)
		return avg_phase
		
	def getAvgGapPhaseDeg(self):
		""" Returns average phase in degrees for all RF gaps in the cavity """
		return self.getAvgGapPhase()*180./math.pi	
		
	def removeAllGapNodes(self):
		""" Remove all rf gaps from this cavity. """
		self.__rfGaps = []
		
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
		self.__rfCavities = []		
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
		
	def setNodes(self, linacNodes):
		""" Set a new set of Linac Nodes. """
		self.__linacNodes = linacNodes
		for node in self.__linacNodes:
			node.setSequence(self)
			
	def addRF_Cavity(self,cav):
		""" Adds the RF cavity to the list inside this sequence. """ 
		if(isinstance(cav, RF_Cavity) == True):
			self.__rfCavities.append(cav)
		else:
			msg = "The Sequence class , method addRF_Cavity(cav)!"
			msg = msg + os.linesep
			msg = msg + "cav is not a subclass of RF_Cavity."
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)				
			
	def getRF_Cavity(self,name):
		""" Returns the cavity instance according to the name """
		for cav in self.__rfCavities:
			if(name == cav.getName()):
				return cav
		return None
		
	def getRF_Cavities(self):
		""" Returns the array with RF cavities. """
		return self.__rfCavities
			
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
