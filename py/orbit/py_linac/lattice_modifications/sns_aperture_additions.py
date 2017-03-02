#!/usr/bin/env python

#--------------------------------------------------------------
# Function to add the aperture class instances to the SNS linac lattice.
# These apertures are not belong to the particular accelerator elements,
# so we created them as markers: MEBT:ChpPlt:Entr and MEBT:ChpPlt:Exit
#--------------------------------------------------------------

import math
import sys
import os

from orbit.py_linac.lattice import LinacApertureNode
from orbit.py_linac.lattice import Quad

def AddMEBTChopperPlatesAperturesToSNS_Lattice(accLattice,aprtNodes):
	"""
	Function will add two Aperture nodes at the entrance and exit of
	MEBT chopper plates. It returns the list of Aperture nodes.
	"""
	x_size = 0.060
	y_size = 0.018
	shape = 3
	node_pos_dict = accLattice.getNodePositionsDict()
	node1 = accLattice.getNodesForName("MEBT:ChpPlt:Entr")[0]
	node2 = accLattice.getNodesForName("MEBT:ChpPlt:Exit")[0]
	for node in [node1,node2]:
		node_name = node.getName()
		(posBefore, posAfter) = node_pos_dict[node]
		apertureNode = LinacApertureNode(shape,x_size/2.0,y_size/2.0,posBefore)
		apertureNode.setName(node_name+":Aprt")
		apertureNode.setSequence(node.getSequence())
		node.addChildNode(apertureNode,node.ENTRANCE)
		aprtNodes.append(apertureNode)
	aprtNodes = sorted(aprtNodes, key = lambda x: x.getPosition(), reverse = False)
	return aprtNodes
		
