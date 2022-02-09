"""
This module contains bump classes for TEAPOT lattice
"""

import os
import math

# import the auxiliary classes
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

# import bump class
from orbit.bumps import simpleBump, TDsimpleBump


class TeapotSimpleBumpNode(DriftTEAPOT):
	""" 
	Kicker node class for TEAPOT lattice
	"""
	def __init__(self, bunch, xbump, xpbump, ybump, ypbump, \
                     name = "bump"):
		"""
		Constructor. Creates a Bump TEAPOT element.
		"""
		DriftTEAPOT.__init__(self, name)
		self.simplebump = simpleBump(bunch, xbump, xpbump, \
                                             ybump, ypbump);
		self.setType("Bump")
		self.setLength(0.0)

	def track(self, paramsDict):
		"""
		Simplebump-teapot class implementation of the
                AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		self.simplebump.bump()


class TDTeapotSimpleBumpNode(DriftTEAPOT):
	""" 
	Kicker node class for TEAPOT lattice
	"""
	def __init__(self, bunch, xbump, xpbump, ybump, ypbump, \
                     waveform, name = "bump"):
		"""
		Constructor. Creates a TDBump TEAPOT element.
		"""
		DriftTEAPOT.__init__(self, name)
		self.TDsimplebump = TDsimpleBump(bunch, xbump, xpbump, \
                                                 ybump, ypbump, waveform);
		self.setType("Bump")
		self.setLength(0.0)

	def track(self, paramsDict):
		"""
		TDSimplebump-teapot class implementation of the
                AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		self.TDsimplebump.bump()
