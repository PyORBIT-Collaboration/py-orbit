"""
This package is a collection of the lattice functions.
Some of these functions could be lattice class methods,
but I think they are too specific. So they will be functions.
"""

import os
import math
import sys


from orbit.py_linac.lattice import BaseLinacNode, Drift, Quad
from orbit.py_linac.lattice import AxisFieldRF_Gap
from orbit.py_linac.lattice import AxisField_and_Quad_RF_Gap

# import acc. nodes
from LinacAccNodes import Quad
from LinacRfGapNodes import AxisFieldRF_Gap

from LinacFieldOverlappingNodes import AxisField_and_Quad_RF_Gap
from LinacFieldOverlappingNodes import OverlappingQuadsNode

def GetGlobalQuadGradient(accLattice,z):
	"""
	The service function for the overlapping fields package.
	Returns the quad field for certain position in the lattice that
	has usual and overlapping quads.
	"""
	node_pos_dict = accLattice.getNodePositionsDict()
	nodes = accLattice.getNodes()
	G = 0.
	(node,index,posBefore,posAfter) = accLattice.getNodeForPosition(z)
	if(isinstance(node,Quad)):
		return node.getParam("dB/dr")
	if(isinstance(node,OverlappingQuadsNode)):
		G = node.getTotalField(z - (posBefore+posAfter)/2)			
		return G
	if(isinstance(node,AxisField_and_Quad_RF_Gap)):
		(z_min,z_max) = node.getZ_Min_Max()
		G = node.getTotalField((z - posBefore)+z_min)			
		return G		
	return G

def GetGlobalRF_AxisField(accLattice,z):
	"""
	The service function for the overlapping RF fields package.
	Returns the RF field on the axis of the RF cavities 
	for certain position in the lattice. If we have 
	the BaseRF_Gap instance we will get 0, because it is 
	an element with zero length.
	"""	
	node_pos_dict = accLattice.getNodePositionsDict()
	nodes =accLattice.getNodes()
	Ez = 0.
	(node,index,posBefore,posAfter) = accLattice.getNodeForPosition(z)
	if(isinstance(node,AxisField_and_Quad_RF_Gap) or isinstance(node,AxisFieldRF_Gap)):
		(z_min,z_max) = node.getZ_Min_Max()
		Ez = node.getEzFiled(z - posBefore + z_min)
		modePhase = node.getParam("mode")*math.pi
		Ez = Ez*math.cos(modePhase)
	return Ez