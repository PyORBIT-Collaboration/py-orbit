"""
Module. Includes a base function that modifies the accelerator lattice
by inserting SC accelerator nodes.
"""

# import auxiliary classes
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, \
	AccActionsContainer, AccNodeBunchTracker

def setSC_Centered_AccNodes(lattice, sc_path_length_min, \
	space_charge_calculator, SC_NodeConstructor):

	"""
	Inserts a set of space charge nodes into the lattice as child nodes
	of the first level accelerator nodes. The SC nodes will be inserted
	at chosen parts of first level AccNode elements such that space
	charge evaluation will be second order accurate.
	The distance between SC nodes must be more than sc_path_length_min.
	The function returns the array of SC nodes as a convenience for the
	user. This is a general function, and SC nodes will require further
	information regarding the nature of the SC nodes to be provided later.
	"""

	accNodes = lattice.getNodes()
	if(len(accNodes) == 0): return

	#-----------------------------------------------
	# nodes_arr[(accNode, part_index, position, path_length)]
	#-----------------------------------------------

	nodes_arr = [(accNodes[0], 0, 0.0, 0.0)]
	length_total = 0.0
	running_path = 0.0

	for accNode in accNodes:
		nParts = accNode.getnParts()
		for ip in range(nParts):
			part_length = accNode.getLength(ip)
			if(part_length > 1.0):
				print "Warning! Node ", accNode.getName(), " has length ", \
					part_length, "m which is greater than 1 m. Space charge \
					algorithm may be innacurate!"
			if(running_path > sc_path_length_min):
				nodes_arr.append((accNode, ip, length_total, running_path))
				running_path = 0.0
			running_path += part_length
			length_total += part_length

	if(len(nodes_arr) > 0):
		running_path = length_total - nodes_arr[len(nodes_arr) - 1][2]
	else:
		running_path = length_total

	#---------------------------------------------------
	# The final SC node in the lattice
	#---------------------------------------------------

	accNode = accNodes[len(accNodes) - 1]
	ip = accNode.getnParts() - 1
	nodes_arr.append((accNode, ip, length_total, running_path))

	#---------------------------------------------------
	# Now we put all SC nodes as childeren of accNodes
	#---------------------------------------------------

	scNodes_arr = []
	for inode in range(len(nodes_arr) - 1):
		(accNode, part_index, position, path_length) = nodes_arr[inode]
		(accNodeNext, part_indexNext, positionNext, path_lengthNext) = \
			nodes_arr[inode + 1]
		scNode = SC_NodeConstructor(space_charge_calculator, \
			accNode.getName() + ":" + str(part_index) + ":")
		NodeLength = (path_length + path_lengthNext) / 2.0
		scNode.setLengthOfSC(NodeLength)
		scNodes_arr.append(scNode)
		accNode.addChildNode(scNode, AccNode.BODY, part_index, AccNode.BEFORE)

	#---------------------------------------------------
	# Set the last SC node
	#---------------------------------------------------

	(accNode, part_index, position, path_length) = \
		nodes_arr[len(nodes_arr) - 1]
	scNode = SC_NodeConstructor(space_charge_calculator, \
		accNode.getName() + ":" + str(part_index) + ":")
	NodeLength = path_length / 2.0
	scNode.setLengthOfSC(NodeLength)
	scNodes_arr.append(scNode)
	accNode.addChildNode(scNode, AccNode.BODY, part_index, AccNode.AFTER)
	return scNodes_arr

