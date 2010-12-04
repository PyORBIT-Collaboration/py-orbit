"""
Module. Includes functions that will modify the accelerator lattice by inserting the SC accelerator nodes.
"""

# import SC acc. nodes
from orbit.space_charge.sc3d import SC3D_AccNode, SC_UniformEllipses_AccNode

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

#import the general SC lattice modification function 
from orbit.space_charge.scLatticeModifications import setSC_General_AccNodes

def setSC3DAccNodes(lattice, sc_path_length_min, space_charge_calculator):
	"""
	It will put a set of a space charge SC2p5D_AccNode into the lattice as child nodes of the first level accelerator nodes.
	The SC nodes will be inserted at the beginning of a particular part of the first level AccNode element.
	The distance between SC nodes should be more than sc_path_length_min, and the boundary is optional.
	The function will return the array of SC nodes as a convenience for the user.
	"""
	scNodes_arr = setSC_General_AccNodes(lattice, sc_path_length_min, space_charge_calculator, SC3D_AccNode)
	for scNode in scNodes_arr:
		scNode.setName(scNode.getName()+"SC3D")
	# initialize the lattice
	lattice.initialize()
	return scNodes_arr
		
def setUniformEllipsesSCAccNodes(lattice, sc_path_length_min, space_charge_calculator):
	"""
	It will put a set of a space charge SC_UniformEllipses_AccNode into the lattice as child nodes of the first level accelerator nodes.
	The SC nodes will be inserted at the beginning of a particular part of the first level AccNode element.
	The distance between SC nodes should be more than sc_path_length_min.
	The function will return the array of SC nodes as a convenience for the user.
	"""
	scNodes_arr = setSC_General_AccNodes(lattice, sc_path_length_min, space_charge_calculator, SC_UniformEllipses_AccNode)
	for scNode in scNodes_arr:
		scNode.setName(scNode.getName()+"UnifEllsSC")
	# initialize the lattice
	lattice.initialize()
	return scNodes_arr	

