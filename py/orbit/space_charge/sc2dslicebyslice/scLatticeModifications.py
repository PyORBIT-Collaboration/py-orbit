"""
Module. Includes functions that will modify the accelerator lattice by inserting the SC accelerator nodes.
"""

# import SC acc. nodes
from orbit.space_charge.sc2dslicebyslice import SC2DSliceBySlice_AccNode

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import the boindary from c++ py module
from spacecharge import Boundary2D

#import the general SC lattice modification function 
from orbit.space_charge.scLatticeModifications import setSC_General_AccNodes

def setSC2DSliceBySliceAccNodes(lattice, sc_path_length_min, space_charge_calculator, boundary = None):
	"""
	It will put a set of a space charge SC2DSliceBySlice_AccNode into the lattice as child nodes of the first level accelerator nodes.
	The SC nodes will be inserted at the beginning of a particular part of the first level AccNode element.
	The distance between SC nodes should be more than sc_path_length_min, and the boundary is optional.
	The function will return the array of SC nodes as a convenience for the user.
	"""
	scNodes_arr = setSC_General_AccNodes(lattice, sc_path_length_min, space_charge_calculator, SC2DSliceBySlice_AccNode)
	for scNode in scNodes_arr:
		scNode.setName(scNode.getName()+"SC2DSliceBySlice")
		scNode.setBoundary(boundary)
	# initialize the lattice
	lattice.initialize()
	return scNodes_arr
		
