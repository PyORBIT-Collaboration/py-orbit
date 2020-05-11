#!/usr/bin/env python

#--------------------------------------------------------------
# The functions for the lattice modifications to replace simple 
# quads (class Quad) with objects that describe the fields of 
# several quads with overlapping fields (class OverlappingQuadsNode). 
# quadrupoles. 
# These combine fields can be chopped into parts to include the 
# elements with zero length. All others non-zero length elements
# can be only drifts.
#--------------------------------------------------------------

import math
import sys
import os
import time

# import from orbit Python utilities
from orbit.utils import orbitFinalize

from orbit.py_linac.lattice import Quad, Drift

from orbit.py_linac.lattice import OverlappingQuadsNode

from rf_quad_overlap_modifications_lib import GetEngeFunction

def Replace_Quads_to_OverlappingQuads_Nodes(\
	accLattice,\
	z_step,\
	accSeq_Names = [],\
	quad_Names = [], \
	EngeFunctionFactory = GetEngeFunction):
	"""
	Function will replace  Quad nodes by OverlappingQuadsNode nodes.
	The replacement will be performed only for specified sequences.
	If the quad names list is empty, all of them will be replaced!
	z_step defines the longitudinal step during the tracking through the 
	elements with overlapping fields.
	The magnetic field longitudinal dependency in quads will be described 
	by Enge Functions that will be produced by the GetEngeFunction function
	by default. The user can supply his/her own factory for these functions.
	"""
	#-----------------------------------------------------------------------------
	#---- if the new drift will be shorter than drift_length_tolerance we ingnore it
	drift_length_tolerance = 0.00000001
	node_pos_dict = accLattice.getNodePositionsDict()
	for accSeq_Name in accSeq_Names:
		accSeq = accLattice.getSequence(accSeq_Name)
		#print "debug ================== STAR seq=",accSeq.getName()
		if(accSeq == None):
			msg  = "The Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes Python function. "
			msg += os.linesep
			msg += "Cannot find the acc. sequence with this name in the lattice!"
			msg += os.linesep
			msg = msg + "accSeq name = " + 	accSeq_Name		
			msg = msg + os.linesep
			msg = msg + "lattice name = " + accLattice.getName()
			msg = msg + os.linesep
			orbitFinalize(msg)
		#--------------------------------------------------------------
		nodes = accSeq.getNodes()
		#---- the start and end positions of the accSeq in the lattice
		accSeq_z_start = node_pos_dict[nodes[0]][0]
		accSeq_z_end = node_pos_dict[nodes[len(nodes)-1]][1]
		#---- just for case: if the nodes are not in the right order
		nodes = sorted(nodes, key = lambda x: x.getPosition(), reverse = False)
		#---- the node to index in the AccSeq dictionary
		node_to_index_in_seq_dict = {}
		for ind in range(len(nodes)):
			node_to_index_in_seq_dict[nodes[ind]] = ind
		#---- quads list for replacement
		quads = accLattice.getQuads(accSeq)
		if(len(quad_Names) > 0):
			quads_tmp = []
			for quad in quads:
				if(quad.getName() in quad_Names):
					quads_tmp.append(quad)
			quads = quads_tmp
		#---- create Enge Functions' dictionary by using the usual quad nodes as keys
		enge_func_quad_dict = {}	
		for quad in quads:
			enge_func_quad_dict[quad] = EngeFunctionFactory(quad)
		#--------------------------------------------
		#---- let's create the array with groups of overlaping quads
		#---- quad_groups_and_ind_arr = [[quad0,quad1,..],pos_start,pos_end,ind_start,ind_end]
		quad_groups_and_ind_arr = Find_Groups_of_Quads(accLattice,accSeq,quads,enge_func_quad_dict,node_to_index_in_seq_dict)
		#--------------------------------------------	
		if(len(quad_groups_and_ind_arr) == 0): return
		#---- let's check that the ends of fields of groups of quads cover the drifts
		for group_ind in range(len(quad_groups_and_ind_arr)):
			[quads_arr,pos_start,pos_end,ind_start,ind_end] = quad_groups_and_ind_arr[group_ind]
			for node_ind in range(ind_start,ind_end+1):
				node = nodes[node_ind]
				if((group_ind == 0 or group_ind == (len(quad_groups_and_ind_arr)-1)) and node in quads_arr):
					(z0_min,z0_max) = enge_func_quad_dict[node].getLimitsZ()
					quad_pos = (node_pos_dict[node][0] + node_pos_dict[node][1])/2
					delta_pos_quad_start = quad_pos + z0_min - accSeq_z_start
					delta_pos_quad_end = quad_pos + z0_max - accSeq_z_end
					res = delta_pos_quad_start < (- drift_length_tolerance)
					res = res or (delta_pos_quad_end > drift_length_tolerance)
					if(res):
						msg  = "The  Replace_Quads_to_OverlappingQuads_Nodes function. STOP. "
						msg += os.linesep
						msg += "We have quad at the beginning or the end of the sequence! quad = "+node.getName()
						msg += os.linesep				
						msg += "quad field positions[m] [from,to] ="+str((delta_pos_quad_start,delta_pos_quad_end))
						msg += os.linesep						
						msg += "The quad's distributed field will be outside the sequence!"
						msg += os.linesep
						msg += "The code cannot handle this situation at this moment."
						msg += os.linesep						
						msg += "You have to change the EngeFunctionFactory to describe the right field function."
						msg += os.linesep	
						orbitFinalize(msg)
				if(node in quads_arr): continue
				if(node.getLength() != 0 and (not isinstance(node,Drift))):
					msg  = "The  Replace_Quads_to_OverlappingQuads_Nodes function. STOP. "
					msg += os.linesep
					msg += "The field of the group of quads covers the elements with L != 0 and different from Drift! "
					msg += os.linesep
					msg += "At this moment we do not know how to handle this situation."
					msg += os.linesep
					quad_names = ""
					for quad in quads_arr:
						quad_names += quad.getName()+" " 
					msg += "Quads = " + quad_names
					msg += os.linesep
					msg += "node=" + node.getName()
					msg += os.linesep				
					msg += "Type = "+ node.getType()
					orbitFinalize(msg)	
		#----------------------------------------------------------------------------
		#-----------------------------------------------------
		#---- Now we are going to build a new lattice with OverlappingQuadsNode
		#---- classes. The space covered by group of quads with overlapping fields
		#---- will be represented by one or several OverlappingQuadsNode instances.
		#---- We may need several OverlappingQuadsNode instances to include a zero 
		#---- length nodes. They will cut the covered space in parts. Each part will
		#---- be represented by one OverlappingQuadsNode instance with the same 
		#---- field source (all overlapping quads).
		#---- We will generate new nodes in an arbitrary order, but at the end 
		#---- we will sort them according to their position.
		#-----------------------------------------------------	
		new_nodes = []
		n_groups = len(quad_groups_and_ind_arr)
		#-------------------------------------------------------------------------
		#---- 1st STEP - cut the length of drift nodes at the beginning and the end
		#----            of quad groups
		for quad_group_ind in range(n_groups):	
			[quads_arr,pos_start,pos_end,ind_start,ind_end] = quad_groups_and_ind_arr[quad_group_ind]
			#---- if the node with index=ind_start or ind_end is a drift 
			#---- its length should be cut by the length of the overlapping 
			#---- region. These drifts are added to the new nodes
			node = nodes[ind_start]
			if(isinstance(node,Drift)):
				(node_pos_start,node_pos_end) = node_pos_dict[node]
				delta = node_pos_end - pos_start
				new_length = abs(node.getLength() - delta)
				#print "debug start group node = ",node.getName()," node.getLength() - delta =",node.getLength() - delta
				if(new_length > drift_length_tolerance):
					if(quad_group_ind > 0):
						#---- we have to check that this node is not the end node of the previous quad group
						#---- if it is true we have to skip it because we already accounted for this node
						#---- earlier when we considered node at the end of the previous group
						[quads_arr0,pos_start0,pos_end0,ind_start0,ind_end0] = quad_groups_and_ind_arr[quad_group_ind-1]
						if(ind_start != ind_end0):
							node_new = Drift(node.getName())
							node_new .setLength(new_length)
							node_new .setPosition(node.getPosition() - delta/2)
							new_nodes.append(node_new)
					else:
						node_new = Drift(node.getName())
						node_new .setLength(new_length)
						node_new .setPosition(node.getPosition() - delta/2)
						new_nodes.append(node_new)
			node = nodes[ind_end]
			if(isinstance(node,Drift)):
				(node_pos_start,node_pos_end) = node_pos_dict[node]
				delta = pos_end - node_pos_start 
				new_length = abs(node.getLength() - delta)
				#print "debug end   group node = ",node.getName()," node.getLength() - delta =",node.getLength() - delta
				if(new_length > drift_length_tolerance):
					if(quad_group_ind < n_groups - 1):
						#---- we have to check that this node is not the start node of the next quad group
						#---- if it is true we have to account for the cat in the length from this group also
						[quads_arr1,pos_start1,pos_end1,ind_start1,ind_end1] = quad_groups_and_ind_arr[quad_group_ind+1]
						if(ind_end == ind_start1):
							delta1 = node_pos_end - pos_start1
							new_length = abs(new_length - delta1)
							node_new = Drift(node.getName())
							node_new.setLength(new_length)
							node_new.setPosition(node.getPosition() + delta/2 - delta1/2 )	
							new_nodes.append(node_new)								
						else:
							node_new = Drift(node.getName())
							node_new.setLength(new_length)
							node_new.setPosition(node.getPosition() + delta/2)	
							new_nodes.append(node_new)							
					else:
						node_new = Drift(node.getName())
						node_new.setLength(new_length)
						node_new.setPosition(node.getPosition() + delta/2)	
						new_nodes.append(node_new)
		#---- 2st STEP - nodes from the beginning to the first group of quads
		[quads_arr,pos_start,pos_end,ind_start,ind_end] = quad_groups_and_ind_arr[0]
		for node_ind in range(0,ind_start):
			#--- we keep the old positions
			new_nodes.append(nodes[node_ind])
		#--------------------------------------------------------------------------
		#---- 3st STEP - nodes from the last quad group to the end of the Acc. Sequence
		[quads_arr,pos_start,pos_end,ind_start,ind_end] = quad_groups_and_ind_arr[n_groups-1]
		for node_ind in range(ind_end+1,len(nodes)):
			#--- we keep the old positions
			new_nodes.append(nodes[node_ind])
		#--------------------------------------------------------------------------
		#---- 4st STEP - nodes between the quad groups
		for quad_group_ind in range(n_groups-1):	
			[quads_arr0,pos_start0,pos_end0,ind_start0,ind_end0] = quad_groups_and_ind_arr[quad_group_ind]
			[quads_arr1,pos_start1,pos_end1,ind_start1,ind_end1] = quad_groups_and_ind_arr[quad_group_ind+1]
			for node_ind in range(ind_end0+1,ind_start1):
				#--- we keep the old positions
				new_nodes.append(nodes[node_ind])
		#--------------------------------------------------------------------------	
		#---- 5st STEP - create OverlappingQuadsNode nodes to cover all the quad groups
		for quad_group_ind in range(n_groups):	
			[quads_arr,group_pos_start,group_pos_end,ind_start,ind_end] = quad_groups_and_ind_arr[quad_group_ind]
			(quads_tmp,zero_length_nodes) = Get_quads_zeroLengthNodes_in_range(accSeq,ind_start,ind_end)
			node_pos_start = group_pos_start - accSeq_z_start
			node_pos_end = group_pos_end - accSeq_z_start			
			pos_start = node_pos_start
			pos_end = node_pos_end
			ovrlp_count = 0
			for ind in range(len(zero_length_nodes)+1):
				node = OverlappingQuadsNode()
				name = quads_arr[0].getName()+":group:"+str(ovrlp_count+1)+":"+node.getType()
				node.setName(name)
				if(ind == len(zero_length_nodes)):
					pos_end = node_pos_end
				else:
					pos_end = zero_length_nodes[ind].getPosition()
				length = pos_end - pos_start
				if(abs(length) < drift_length_tolerance):
					if(ind < len(zero_length_nodes)):
						new_nodes.append(zero_length_nodes[ind])
					pos_start = pos_end
					continue
				pos = (pos_end + pos_start)/2
				node.setLength(length)
				node.setPosition(pos)
				for quad in quads_arr:
					node.addQuad(quad,enge_func_quad_dict[quad],quad.getPosition()-pos_start)
				node.setZ_Step(z_step)
				nParts = int(length/z_step)+1
				node.setnParts(nParts)					
				new_nodes.append(node)
				ovrlp_count += 1
				if(ind < len(zero_length_nodes)):
					new_nodes.append(zero_length_nodes[ind])
				pos_start = pos_end
		#--------------------------------------------------------------------------
		new_nodes = sorted(new_nodes, key = lambda x: x.getPosition(), reverse = False)
		#--------------------------------------------------------------------------	
		#------------------------------------------
		#---- let's replace all nodes in the AccSeq by the new set
		accSeq.removeAllNodes()
		for node in new_nodes:
			accSeq.addNode(node)
	#---- new set of nodes for the lattice
	new_latt_nodes = []
	for accSeq in accLattice.getSequences():
		new_latt_nodes += accSeq.getNodes()
	accLattice.setNodes(new_latt_nodes)
	accLattice.initialize()
	#------- debug START printing of new nodes and their positions in the lattice
	"""
	node_pos_dict = accLattice.getNodePositionsDict()
	for accSeq_Name in accSeq_Names:
		accSeq = accLattice.getSequence(accSeq_Name)
		nodes = accSeq.getNodes()
		accSeq_z_start = node_pos_dict[nodes[0]][0]
		#--------------------------------------------------------------------------
		for node in nodes:
			pos = node.getPosition()
			(pos_start,pos_end) = node_pos_dict[node]
			delta = pos - ((pos_start+pos_end)/2 - accSeq_z_start)
			if(abs(delta) >  drift_length_tolerance):
				print "debug new node=",node.getName()," pos=",node.getPosition()," (pos_start,pos_end)=",node_pos_dict[node]," delta=",delta
	"""
	#------- debug STOP  printing of new nodes and their positions in the lattice
			
def Find_Groups_of_Quads(accLattice,accSeq,quads,enge_func_quad_dict,node_to_index_in_seq_dict):
	"""
	This function will find the group of quads with the fields overlapping each other.
	It returns the quads that should be taken into account.
	quad_groups_and_ind_arr = [[quad0,quad1,..],pos_start,pos_end,ind_start,ind_end]]
	"""
	nodes = accSeq.getNodes()
	node_pos_dict = accLattice.getNodePositionsDict()
	n_quads = len(quads)
	quad_groups_and_ind_arr = []
	quads_arr = []
	for quad_ind in range(n_quads):	
		quad0 = quads[quad_ind]
		(pos0_start,pos0_end) = node_pos_dict[quad0]
		(z0_min,z0_max) = enge_func_quad_dict[quad0].getLimitsZ()
		pos0_center = (pos0_start+pos0_end)/2 
		quads_arr.append(quad0)
		if(len(quads_arr) == 1):
			(node,index,posBefore,posAfter) = GetNodeInAccSeqForPosition(accLattice,nodes,pos0_center + z0_min)
			quad_groups_and_ind_arr.append([quads_arr,pos0_center+z0_min,0.,index,-1])	
		quad1 = quad0
		if(quad_ind != (n_quads-1)):
			quad1 = quads[quad_ind+1]
		(pos1_start,pos1_end) = node_pos_dict[quad1]
		(z1_min,z1_max) = enge_func_quad_dict[quad1].getLimitsZ()
		pos1_center = (pos1_start+pos1_end)/2 
		if(pos0_center+z0_max < pos1_center+z1_min or quad1 == quad0):
			(node,index,posBefore,posAfter) = GetNodeInAccSeqForPosition(accLattice,nodes,pos0_center+z0_max)
			[arr_tmp,pos_start,pos_end,ind_start,ind_end] = quad_groups_and_ind_arr[len(quad_groups_and_ind_arr)-1]
			pos_end = pos0_center + z0_max
			ind_end = index
			quad_groups_and_ind_arr[len(quad_groups_and_ind_arr)-1] = [arr_tmp,pos_start,pos_end,ind_start,ind_end] 
			quads_arr	= []			
	#----------------------------------- DEBUG PRINTING START
	"""
	for [quads_arr,pos_start,pos_end,ind_start,ind_end] in quad_groups_and_ind_arr:
		print "debug new group N quads=",len(quads_arr)," (pos_start,pos_end)=",(pos_start,pos_end)," (ind_start,ind_end)=",(ind_start,ind_end)
		for quad in quads_arr:
			(pos_start,pos_end) = node_pos_dict[quad]
			(z_min,z_max) = enge_func_quad_dict[quad].getLimitsZ()
			pos_center = (pos_start+pos_end)/2 
			print "debug          quad = ",quad.getName()," (pos_center+z_min,pos_center+z_max) =",(pos_center+z_min,pos_center+z_max)
	"""
	#----------------------------------- DEBUG PRINTING END		
	return quad_groups_and_ind_arr
				
				
def GetNodeInAccSeqForPosition(accLattice,nodes,z):
	"""
	It is a local convenience function. It returns the node in the AccSeq 
	which coordinates cover the z-position. The position z is a position
	in the lattice.
	This function is different from the AccLattice method getNodeForPosition(z),
	because it limits the nodes' array to speed up the finding process.
	"""
	node_pos_dict = accLattice.getNodePositionsDict()
	index0 = 0
	index1 = len(nodes) - 1
	(posBefore0, posAfter0) = node_pos_dict[nodes[index0]]
	(posBefore1, posAfter1) = node_pos_dict[nodes[index1]]
	index = 0
	if(z <= posBefore0): return (nodes[index0],index0,posBefore0,posAfter0)
	if(z >= posAfter1): return (nodes[index1],index1,posBefore1,posAfter1)
	while(index0 != index1):			
		index = (index0 + index1)/2
		#print "debug z=",z," index0=",index0," index1=",index1," index=",index," (posBefore0, posAfter0)=",(posBefore0, posAfter0)," (posBefore1, posAfter1)=",(posBefore1, posAfter1)
		(posBefore, posAfter) = node_pos_dict[nodes[index]]
		if(z < posBefore):
			index1 = index
			(posBefore1, posAfter1) = node_pos_dict[nodes[index1]]
		elif(z > posAfter):
			index0 = index
			(posBefore0, posAfter0) = node_pos_dict[nodes[index0]]
		elif(z >= posBefore and z <= posAfter):
			break
		if(z >= posBefore0 and z <= posAfter0):
			index =index0 
			break
		if(z >= posBefore1 and z <= posAfter1):
			index =index1
			break		
	#-------------------------------------------
	node = nodes[index]
	(posBefore, posAfter) = node_pos_dict[node]
	return (node,index,posBefore,posAfter)
		

def Get_quads_zeroLengthNodes_in_range(accSeq,node_ind_start,node_ind_end):
	"""
	Returns all quads and zero-length nodes in this index range.
	It also checks that all elements inside this range has zero length or they 
	are drifts of quads.
	"""
	nodes = accSeq.getNodes()
	zero_length_nodes = []
	child_nodes = []
	quads = []
	for node_ind in range(node_ind_start,node_ind_end+1):
		node = nodes[node_ind]
		children_arr = node.getBodyChildren()
		if(len(children_arr) > 0):
			#print "debug ========= parent node=",node.getName()," pos = ",node.getPosition()
			for child in children_arr:
				if(child.getLength() == 0.):
					child_nodes.append(child)
					#print "      debug child=",child.getName()," pos=",child.getPosition()
		length = node.getLength() 
		if(length == 0.):
			zero_length_nodes.append(node)
		else:
			if(isinstance(node,Quad)):
				quads.append(node)
			else:
				if(not isinstance(node,Drift)):
					msg  = "The Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes function. "
					msg += "This Acc. Sequence has an element that "
					msg += os.linesep
					msg += "1. has non-zero length"
					msg += os.linesep
					msg += "2. not a quad"
					msg += os.linesep
					msg += "3. not a drift"
					msg += os.linesep				
					msg += "This function does not know how to handle this elelement!"
					msg += os.linesep			
					msg = msg + "Acc Sequence =" + accSeq.getName()			
					msg = msg + os.linesep					
					msg = msg + "Acc element =" + node.getName()			
					msg = msg + os.linesep		
					msg = msg + "Acc element type =" + node.getType()			
					msg = msg + os.linesep
					orbitFinalize(msg)
	zero_length_nodes += child_nodes
	zero_length_nodes = sorted(zero_length_nodes, key = lambda x: x.getPosition(), reverse = False)
	return (quads,zero_length_nodes)
