"""
Module. Includes functions that will modify the accelerator lattice
by inserting one teapot node accelerator node.
"""

# import the auxiliary classes
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode,\
AccActionsContainer, AccNodeBunchTracker

# import RF accelerator nodes
from orbit.rf_cavities import RFNode

# import teapot drift class
from orbit.teapot import DriftTEAPOT

def addRFNode(lattice, position, rf_node):
	"""
	This will put one rf cavity node in the lattice
	"""
	length_tolerance = 0.0001
	lattice.initialize()
	position_start = position
	position_stop = position + rf_node.getLength()
	(node_start_ind, node_stop_ind, z, ind) = (-1, -1, 0., 0)
	for node in lattice.getNodes():
		if(position_start >= z and\
		position_start <= z + node.getLength()):
			node_start_ind = ind
		if(position_stop  >= z and\
		position_stop  <= z + node.getLength()):
			node_stop_ind  = ind
		ind += 1
		z += node.getLength()
	#-------Now we check that between start and end
	#-------we have only non-modified drift elements
	#-------If the rf was added first - that is a problem.
	#-------The collimation should be added first.
	for node in lattice.getNodes()[node_start_ind:node_stop_ind + 1]:
		#print "debug node = ", node.getName(),\
		#" type = ", node.getType()," L = ", node.getLength()
		if(not isinstance(node, DriftTEAPOT)):
			print "Non-drift node = ", node.getName(),\
			" type = ", node.getType()," L = ", node.getLength()
			orbitFinalize("We have a non-drift element at the\
			location of the rf node! Stop!")
			# Make array of nodes from rf node in the center and possibly
	# two drifts if their length is more than length_tollerance [m]
	nodes_new_arr = [rf_node,]
	drift_node_start = lattice.getNodes()[node_start_ind]
	drift_node_stop  = lattice.getNodes()[node_stop_ind]
	#-------Now we create two drift nodes surrounding the rf node.
	#-------If the length of one of these additional drifts is
	#-------less than length_tollerance [m] we skip this drift
	if(position_start > lattice.getNodePositionsDict()\
	[drift_node_start][0] +  length_tolerance):
		drift_node_start_new = DriftTEAPOT(drift_node_start.getName())
		drift_node_start_new.setLength(position_start\
		- lattice.getNodePositionsDict()[drift_node_start][0])
		nodes_new_arr.insert(0, drift_node_start_new)
	if(position_stop < lattice.getNodePositionsDict()\
	[drift_node_stop][1] - length_tolerance):
		drift_node_stop_new = DriftTEAPOT(drift_node_stop.getName())
		drift_node_stop_new.setLength(lattice.getNodePositionsDict\
		()[drift_node_stop][1] - position_stop)
		nodes_new_arr.append(drift_node_stop_new)
	#-------Now modify the lattice by replacing the chosen part
	#-------by the new nodes
	lattice.getNodes()[node_start_ind:node_stop_ind + 1] = nodes_new_arr
	# initialize the lattice
	lattice.initialize()
