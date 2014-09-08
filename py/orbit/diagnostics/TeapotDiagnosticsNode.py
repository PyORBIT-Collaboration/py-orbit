"""
This module is a collimator node class for TEAPOT lattice
"""

import os
import math

# import the auxiliary classes
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

# import Diagnostics classes
from diagnostics import StatLats, StatLatsSetMember
from diagnostics import Moments, MomentsSetMember

#import Bunch diagnostics
from bunch import BunchTuneAnalysis

class TeapotStatLatsNode(DriftTEAPOT):
	""" 
	The statlats node class for TEAPOT lattice
	"""
	def __init__(self, filename , name = "statlats no name"):
		"""
		Constructor. Creates the StatLats TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,name)
		self.statlats = StatLats(filename)
		self.setType("statlats teapot")
		self.setLength(0.0)
		self.position = 0.0
		self.lattlength = 0.0
		self.file_out = open(filename,"w")

	def track(self, paramsDict):
		"""
		The statlats-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		self.statlats.writeStatLats(self.position,bunch,self.lattlength)
			
	def setPosition(self,pos):
		self.position = pos

	def closeStatLats(self):
		self.file_out.close()

	def setLatticeLength(self, lattlength):
		self.lattlength = lattlength

class TeapotStatLatsNodeSetMember(DriftTEAPOT):
	"""
		The statlats node class for TEAPOT lattice
	"""
	def __init__(self, file, name = "statlats no name"):
		"""
			Constructor. Creates the StatLats TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,name)
		self.statlats = StatLatsSetMember(file)
		self.setType("statlats teapot")
		self.setLength(0.0)
		self.position = 0.0
		self.lattlength = 0.0
		self.active = True
		self.file = file
	
	def track(self, paramsDict):
		"""
			The statlats-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		if(self.active):
			length = self.getLength(self.getActivePartIndex())
			bunch = paramsDict["bunch"]
			self.statlats.writeStatLats(self.position,bunch,self.lattlength)
	
	def setPosition(self,pos):
		self.position = pos

	def setLatticeLength(self, lattlength):
		self.lattlength = lattlength

	def activate(self):
		self.active = True

	def deactivate(self):
		self.active = False

	def resetFile(self, file):
		self.file = file
		self.statlats.resetFile(self.file)


class TeapotMomentsNode(DriftTEAPOT):
	"""
	The moments node class for TEAPOT lattice
	"""
	def __init__(self, filename, order, nodispersion = "true", name = "moments no name"):
		"""
		Constructor. Creates the StatLats TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,name)
		self.moments = Moments(filename, order, name)
		self.setType("moments teapot")
		self.setLength(0.0)
		self.position = 0.0
		self.lattlength = 0.0
		self.file_out = open(filename,"w")
	
	def track(self, paramsDict):
		"""
		The moments-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		self.moments.writeMoments(self.position,bunch,self.lattlength)
	
	def setPosition(self,pos):
		self.position = pos
	
	def closeMoments(self):
		self.file_out.close()
	
	def setLatticeLength(self, lattlength):
		self.lattlength = lattlength
	


class TeapotMomentsNodeSetMember(DriftTEAPOT):
	"""
	The moments node class for TEAPOT lattice
	"""
	def __init__(self, file, order, nodispersion = "true", name = "moments no name"):
		"""
		Constructor. Creates the Moments TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,str(name))
		
		self.file = file
		self.moments = MomentsSetMember(self.file, order, nodispersion="true")
		self.setType("moments teapot")
		self.setLength(0.0)
		self.position = 0.0
		self.lattlength = 0.0
		self.active = True
		
	
	def track(self, paramsDict):
		"""
		The moments-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		if(self.active):
			length = self.getLength(self.getActivePartIndex())
			bunch = paramsDict["bunch"]
			self.moments.writeMoments(self.position, bunch, self.lattlength)
	
	def setPosition(self,pos):
		self.position = pos
	
	def setLatticeLength(self, lattlength):
		self.lattlength = lattlength

	def activate(self):
		self.active = True

	def deactivate(self):
		self.active = False

	def resetFile(self, file):
		self.file = file
		self.moments.resetFile(self.file)
		

class TeapotTuneAnalysisNode(DriftTEAPOT):
	
			
	def __init__(self, name = "tuneanalysis no name"):
		"""
		Constructor. Creates the StatLats TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,name)
		self.bunchtune = BunchTuneAnalysis()
		self.setType("tune calculator teapot")
		self.lattlength = 0.0
		self.setLength(0.0)
		self.position = 0.0
	
	def track(self, paramsDict):
		"""
		The bunchtuneanalysis-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		self.bunchtune.analyzeBunch(bunch)
	
	def setPosition(self,pos):
		self.position = pos
	
	def setLatticeLength(self, lattlength):
		self.lattlength = lattlength

	def assignTwiss(self, betax, alphax, etax, etapx, betay, alphay):
		self.bunchtune.assignTwiss(betax, alphax, etax, etapx, betay, alphay)
