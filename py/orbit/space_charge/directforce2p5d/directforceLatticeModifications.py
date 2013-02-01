"""
Module. Includes functions that will modify the accelerator lattice by inserting the SC accelerator nodes.
"""

# import Direct Force acc. nodes
from orbit.space_charge.directforce2p5d import DirectForce2p5D_AccNode

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

#import the general SC lattice modification function 
from orbit.space_charge.scLatticeModifications import setSC_General_AccNodes

def setDirectForce2p5DAccNodes(lattice, sc_path_length_min, space_charge_calculator):
	"""
	It will put a set of a space charge SC2p5D_AccNode into the lattice as child nodes of the first level accelerator nodes.
	The SC nodes will be inserted at the beginning of a particular part of the first level AccNode element.
	The distance between SC nodes should be more than sc_path_length_min	The function will return the array of SC nodes as a convenience for the user.
	"""
	scNodes_arr = setSC_General_AccNodes(lattice, sc_path_length_min, space_charge_calculator, DirectForce2p5D_AccNode)
	for scNode in scNodes_arr:
		scNode.setName(scNode.getName()+"DirectForce2p5D")
	
	# initialize the lattice
	lattice.initialize()
	return scNodes_arr
