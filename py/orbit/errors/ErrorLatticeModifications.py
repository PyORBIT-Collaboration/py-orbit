"""
Module. Includes functions that will modify the accelerator lattice
by inserting an error node into a teapot accelerator node.
"""

# import the auxiliary classes
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode,\
AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

def addErrorNode(lattice, position, Error_Node):
	"""
	This will put one error node into the lattice
	"""
	length_tolerance = 0.0001
	lattice.initialize()
	position_start = position
	position_stop  = position + Error_Node.getLength()
	(node_start_ind, node_stop_ind, z, ind) = (-1, -1, 0.0, 0)
	for node in lattice.getNodes():
		if(position_start >= z and position_start <= z + node.getLength()):
			node_start_ind = ind
		if(position_stop  >= z and position_stop  <= z + node.getLength()):
			node_stop_ind  = ind
		ind += 1
		z += node.getLength()

	"""
	Check that between start and end there are only non-modified 
	drift elements. If space charge was added first - that is a problem.
	The collimation should be added first.
	"""

	for node in lattice.getNodes()[node_start_ind:node_stop_ind + 1]:
		"""
		print "debug node = ", node.getName(), " type = ", node.getType(),\
		" L = ", node.getLength()
		"""
		if(not isinstance(node, DriftTEAPOT)):
			print "Non-drift node = ", node.getName(), " type = ",\
			node.getType(), " L = ", node.getLength()
			orbitFinalize("We have non-drift element at the place of \
				the error node! Stop!")
			"""
			if(node.getNumberOfChildren() != 4):
				print "Node = ", node.getName()," type = ", node.getType(),\
				" L = ", node.getLength(), " N child nodes = ",\
				node.getNumberOfChildren()
			orbitFinalize("Drift element was modified with additional \
			functionality (SC or something else)! Add collimation first! \
			Stop!")
			"""

	"""
	Make array of nodes with error node in the center and two possible drifts,
	if their length is more than length_tollerance [m]
	"""

	nodes_new_arr = [Error_Node,]
	drift_node_start = lattice.getNodes()[node_start_ind]
	drift_node_stop  = lattice.getNodes()[node_stop_ind]

	"""
	Now create two drift nodes: before and after the error node.
	If the length of either of these additional drifts is less than
	length_tollerance [m], skip this drift.
	"""

	if(position_start > lattice.getNodePositionsDict()[drift_node_start][0] +\
		length_tolerance):
		drift_node_start_new = DriftTEAPOT(drift_node_start.getName())
		drift_node_start_new.setLength(position_start -\
			lattice.getNodePositionsDict()[drift_node_start][0])
		nodes_new_arr.insert(0, drift_node_start_new)
	if(position_stop  < lattice.getNodePositionsDict()[drift_node_stop][1] -\
		length_tolerance):
		drift_node_stop_new = DriftTEAPOT(drift_node_stop.getName())
		drift_node_stop_new.setLength(lattice.getNodePositionsDict()[drift_node_stop][1] -\
			position_stop)
		nodes_new_arr.append(drift_node_stop_new)

	"""
	Now modify the lattice by replacing the old part with the new nodes
	"""

	lattice.getNodes()[node_start_ind:node_stop_ind + 1] = nodes_new_arr
	"""
	Initialize the lattice
	"""
	lattice.initialize()

def addErrorNodeAsChild(lattice, AccNode, Error_Node):
	AccNode.addChildNode(Error_Node, AccNode.BODY, 0, AccNode.BEFORE)
	lattice.initialize()

def addErrorNodeAsChild_I(lattice, AccNode, Error_Node):
	AccNode.addChildNode(Error_Node, AccNode.ENTRANCE)
	lattice.initialize()

def addErrorNodeAsChild_F(lattice, AccNode, Error_Node):
	AccNode.addChildNode(Error_Node, AccNode.EXIT)
	lattice.initialize()
