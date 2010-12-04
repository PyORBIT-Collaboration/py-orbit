"""
Module. Includes a base function that will modify the accelerator lattice by inserting the SC accelerator nodes.
"""

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

def setSC_General_AccNodes(lattice, sc_path_length_min, space_charge_calculator, SC_NodeConstructor):
	"""
	It will put a set of a space charge nodes into the lattice as child nodes of the first level accelerator nodes.
	The SC nodes will be inserted at the beginning of a particular part of the first level AccNode element.
	The distance between SC nodes should be more than sc_path_length_min. The function will return 
	the array of SC nodes as a convenience for the user. This is a general function, and SC nodes will need 
	specific information which will be provided for them later according to the specific nature of SC nodes.
	"""
	accNodes = lattice.getNodes()
	if(len(accNodes) == 0): return
	#-----------------------------------------------
	# nodes_arr[(accNode, part_index, position, path_length)] 
	#-----------------------------------------------
	nodes_arr = []
	length_total = 0.
	running_path = 0.
	for accNode in accNodes:
		nParts = accNode.getnParts()
		for ip in range(nParts):
			part_length = accNode.getLength(ip)
			if(running_path > sc_path_length_min):
				nodes_arr.append((accNode,ip,length_total,running_path))
				running_path = 0.
			running_path += part_length
			length_total += part_length
	rest_length = length_total - nodes_arr[len(nodes_arr) - 1][2]
	#the first SC node in the beginning of the lattice
	nodes_arr.insert(0,(accNodes[0],0,0.,rest_length))
	#---------------------------------------------------
	# Now we put all SC nodes as a childeren of accNodes
	#---------------------------------------------------
	scNodes_arr = []
	for inode in range(len(nodes_arr)-1):
		(accNode, part_index, position, path_length) = nodes_arr[inode]
		(accNodeNext, part_indexNext, positionNext, path_lengthNext) = nodes_arr[inode+1]
		scNode = SC_NodeConstructor(space_charge_calculator,accNode.getName()+":"+str(part_index)+":")
		scNode.setLengthOfSC(path_lengthNext)
		scNodes_arr.append(scNode)
		accNode.addChildNode(scNode,AccNode.BODY,part_index,AccNode.BEFORE)
	#set the last SC node
	(accNode, part_index, position, path_length) = nodes_arr[len(nodes_arr)-1]
	scNode = SC_NodeConstructor(space_charge_calculator,accNode.getName()+":"+str(part_index)+":")
	scNode.setLengthOfSC(rest_length)
	scNodes_arr.append(scNode)
	accNode.addChildNode(scNode,AccNode.BODY,part_index,AccNode.BEFORE)
	return scNodes_arr

