"""
Module. Includes classes for 1D longidutinal space charge accelerator nodes.
"""

import sys
import os
import math


# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

#import longitudinal space charge package
from spacecharge import LSpaceChargeCalc

class SC1D_AccNode(DriftTEAPOT):

	def __init__(self, b_a, length, nMacrosMin, useSpaceCharge, nBins, name = "long sc node"):
		"""
			Constructor. Creates the Foil TEAPOT element.
		"""
		DriftTEAPOT.__init__(self,name)
		self.lspacecharge = LSpaceChargeCalc(b_a, length, nMacrosMin, useSpaceCharge, nBins)
		self.setType("long sc node")
		self.setLength(0.0)
		
	
	def trackBunch(self, bunch):
		"""
			The foil-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		self.lspacecharge.trackBunch(bunch)		#put the track method here
		
	def track(self, paramsDict):
		"""
			The foil-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		self.lspacecharge.trackBunch(bunch)		#put the track method here

	def assignImpedance(self, py_cmplx_arr):
		#size = size(py_cmplx_arr);
		#Py_complex cmplx;
		#PyObject* py_cmplx;	
		self.lspacecharge.assignImpedance(py_cmplx_arr)
		print "Assigning the impadance array"