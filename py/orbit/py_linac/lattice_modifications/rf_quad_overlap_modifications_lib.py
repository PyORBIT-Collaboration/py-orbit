#!/usr/bin/env python

#--------------------------------------------------------------
# The functions for the lattice modifications to replace the base 
# RF gap nodes with The Axis Fields nodes that can overlap with 
# quadrupoles. 
# The initial base RF gap nodes have zero length, and the length
# of new nodes is defined by the RF fields on the RF gap axis.
# We also have to adjust these lengths to avoid overlapping RF fields. 
# This correction will be on the level of microns, and it will 
# not change the strength of the fields.
# The RF fields can overlap the zero length elements like markers or 
# dipole correctors. The whole RF gap field will be chopped into pieces.
# The RF fields also can overlap the fields one or several quads, and 
# these quads' fields will be included into the RF gap node to use
# during the bunch tracking.
#--------------------------------------------------------------

import math
import sys
import os
import time

# import from orbit Python utilities
from orbit.utils import orbitFinalize

from orbit.py_linac.lattice import LinacApertureNode
from orbit.py_linac.lattice import Quad, Drift
from orbit.py_linac.lattice import BaseRF_Gap, AxisFieldRF_Gap
from orbit.py_linac.lattice import OverlappingQuadsNode
from orbit.py_linac.lattice import AxisField_and_Quad_RF_Gap

from orbit.py_linac.overlapping_fields import EngeFunction

from orbit_utils import Function

from rf_models_modifications_lib import Make_AxisFieldRF_Gaps_and_Find_Neihbor_Nodes

def GetEngeFunction(quad):
	"""
	This is an example of a EngeFunctionFactory function.
	It returns the EngeFunction function for the instance of Quad node.
	The quad should have the aperture parameter.
	It is used in the quad+rf lattice transformation by default.
	User can prepare his/her own EngeFunctionFactory function.
	"""
	length_param = quad.getLength()
	if(quad.hasParam("aperture")):
		acceptance_diameter_param = quad.getParam("aperture")
		cutoff_level = 0.001
		func = EngeFunction(length_param,acceptance_diameter_param,cutoff_level)
		return func
	else:
		msg  = "Inside the Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes Python function. "
		msg += os.linesep
		msg += "Cannot create the EngeFunction for the quad!"
		msg += os.linesep
		msg = msg + "quad name = " + 	quad.getName()	
		msg = msg + os.linesep		
		msg = msg + "It does not have the aperture parameter!"
		msg = msg + os.linesep		
		orbitFinalize(msg)
	return None

def Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(\
	accLattice,\
	z_step,\
	dir_location="",\
	accSeq_Names = [],\
	cavs_Names = [], \
	EngeFunctionFactory = GetEngeFunction):
	"""
	Function will replace  BaseRF_Gap nodes by AxisField_and_Quad_RF_Gap.
	It is assumed that AxisField_and_Quad_RF_Gap nodes do not overlap any
	others nodes (except Drifts). The location of the axis field 
	data files are specified by dir_location input variable.
	The replacement will be performed only for specified sequences.
	If the cavities list is empty, all of them will be replaced!
	If you want to replace the nodes in a particular cavity please specify it!
	The dir_location is the location of the directory with the axis field
	files.
	We assume that the RF gap field overlaps only with the fields of 
	two neighboring quads.
	z_step defines the longitudinal step during the tracking through the 
	elements with overlapping fields.
	The magnetic field longitudinal dependency in quads will be described 
	by Enge Functions that will be produced by the GetEngeFunction function
	by default. The user can supply his/her own factory for these functions.
	"""
	#-----------------------------------------------------------------------------
	"""
	node_pos_dict = accLattice.getNodePositionsDict()
	#for node_ind in range(len(accLattice.getNodes())):
	for node_ind in range(199,298):
		node = accLattice.getNodes()[node_ind]
		#if(not isinstance(node,Quad)): continue
		(pos_start,pos_end) = node_pos_dict[node]
		print "debug ind=",node_ind," node=",node.getName()," (pos_start,pos_end)=",(pos_start,pos_end)	
	"""
	#-----------------------------------------------------------------------------
	#---- rf_length_tolerance RF fields should not overlap more than this value
	rf_length_tolerance = 0.0001
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
		#--- just for case: if the nodes are not in the right order
		nodes = sorted(nodes, key = lambda x: x.getPosition(), reverse = False)		
		#---- create Enge Functions' dictionary by using the usual quad nodes as keys
		enge_func_quad_dict = {}		
		quads = accLattice.getQuads(accSeq)
		for quad in quads:
			enge_func_quad_dict[quad] = EngeFunctionFactory(quad)
		#--------------------------------------------
		cavs = accSeq.getRF_Cavities()
		if(len(cavs_Names) > 0):
			cavs_tmp = []
			for cav in cavs:
				if(cav.getName() in cavs_Names):
					cavs_tmp.append(cav)
			cavs = cavs_tmp
		#---- let's check that all rf gaps are BaseRF_Gap instances
		#---- and create the dictionaries to account for new rf gaps later
		#---- rf_gap_to_cavs_dict[rf_gap] = cav
		#---- new_rf_gaps_arr_in_cav_dict[cav] = [new_AxisField_and_Quad_RF_Gap_0,..]
		rf_gap_to_cavs_dict = {}
		new_rf_gaps_arr_in_cav_dict = {}
		for cav in cavs:
			rf_gaps = cav.getRF_GapNodes()
			new_rf_gaps_arr_in_cav_dict[cav] = []
			for rf_gap in rf_gaps:
				rf_gap_to_cavs_dict[rf_gap] = cav
				if(not isinstance(rf_gap,BaseRF_Gap)):
					msg  = "The Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes function. "
					msg += "You are trying to replace the RF gap which is not BaseRF_Gap instance!"
					msg += os.linesep
					msg = msg + "RF Gap =" + rf_gap.getName()				
					msg = msg + os.linesep
					msg = msg + "Type of gap node=" + rf_gap.getType()
					msg = msg + os.linesep
					orbitFinalize(msg)
		#---- af_rf_gap_dict[rf_gap] = AxisFieldRF_Gap(rf_gap)
		#---- rf_gap_ind_up_down_arr[[rf_gap,gap_ind,drift_down_ind,drift_up_ind],...]
		(af_rf_gap_dict,rf_gap_ind_up_down_arr)	 = Make_AxisFieldRF_Gaps_and_Find_Neihbor_Nodes(rf_length_tolerance,accLattice,accSeq,dir_location,cavs)
		if(len(rf_gap_ind_up_down_arr) == 0):
			msg  = "The Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes function. "
			msg += "This Acc. Sequence does not have BaseRF_Gaps!"
			msg += os.linesep
			msg += "It is better to use another method of the lattice modification!"
			msg += os.linesep			
			msg = msg + "Acc Sequence =" + accSeq.getName()			
			msg = msg + os.linesep
			orbitFinalize(msg)
		#-----------------------------------------------------
		#---- Now we are going to build a new lattice with AxisField_and_Quad_RF_Gap 
		#---- and OverlappingQuadsNode classes. Even drifts will be replaced with 
		#---- these nodes to simplify the lattice structure. The original lattice
		#---- should have only drifts, quads, and BaseRF_Gaps as elements with
		#---- non-zero length. If the user wants to include other sources of the EM
		#---- fields, he/she has to create or modify the classes responsible for
		#---- handling the overlapping EM fields.
		#---- We will generate new nodes in an arbitrary order, but at the end 
		#---- we will sort them according to their position.
		#-----------------------------------------------------	
		new_nodes = []
		#-------------------------------------------------------------------------
		#---- 1st STEP - nodes from the beginning to the first RF gap
		[rf_gap,gap_ind,drift_down_ind,drift_up_ind] = rf_gap_ind_up_down_arr[0]
		node_ind_start = 0
		node_ind_end = drift_down_ind
		(node_pos_start,pos_tmp) = node_pos_dict[nodes[node_ind_start]]
		accSeq_z_start = node_pos_start
		(rf_gap_pos_start, rf_gap_pos_end) = node_pos_dict[rf_gap]
		(z_gap_min,z_gap_max) = af_rf_gap_dict[rf_gap].getZ_Min_Max()
		(rf_gap_pos_start, rf_gap_pos_end) = (rf_gap_pos_start+z_gap_min, rf_gap_pos_end+z_gap_max)
		#print "debug rf gap=",rf_gap.getName()," (rf_gap_pos_start, rf_gap_pos_end)=",(rf_gap_pos_start, rf_gap_pos_end)
		#print "debug      (node_ind_start,node_ind_end) = ",(node_ind_start,node_ind_end)
		#for ind in range(node_ind_start,rf_gap_ind_up_down_arr[0][1]+1):
		#	print "debug    ind =",ind," node=",nodes[ind].getName()," (node_pos_start,node_pos_end)=",node_pos_dict[nodes[ind]]
		node_pos_start = node_pos_start - accSeq_z_start
		node_pos_end = rf_gap_pos_start - accSeq_z_start
		#---- if the length of the node is equal to 0. we will add only the zero length elements
		(local_quads,zero_length_nodes) = Get_quads_zeroLengthNodes_in_range(accSeq,node_ind_start,node_ind_end)
		if(abs(node_pos_start - node_pos_end) > drift_length_tolerance):
			#print "debug before n_quads=",len(local_quads)," n_markers=",len(zero_length_nodes)
			pos_start = node_pos_start
			pos_end = node_pos_end
			ovrlp_count = 0
			for ind in range(len(zero_length_nodes)+1):
				name = rf_gap.getName()+":Before:"+str(ovrlp_count+1)+":OVRLPQ"
				node = OverlappingQuadsNode(name)
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
				for quad in local_quads:
					node.addQuad(quad,enge_func_quad_dict[quad],quad.getPosition()-pos_start)
				node.setZ_Step(z_step)
				nParts = int(length/z_step)+1
				node.setnParts(nParts)				
				new_nodes.append(node)
				ovrlp_count += 1
				if(ind < len(zero_length_nodes)):
					new_nodes.append(zero_length_nodes[ind])
				pos_start = pos_end
		else:
			new_nodes += zero_length_nodes
		#--------------------------------------------------------------------------
		#---- 2st STEP - nodes from the last RF gap to the end of the Acc. Sequence
		[rf_gap,gap_ind,drift_down_ind,drift_up_ind] = rf_gap_ind_up_down_arr[len(rf_gap_ind_up_down_arr)-1]
		node_ind_start = drift_up_ind
		node_ind_end = len(nodes) - 1
		(rf_gap_pos_start, rf_gap_pos_end) = node_pos_dict[rf_gap]
		(z_gap_min,z_gap_max) = af_rf_gap_dict[rf_gap].getZ_Min_Max()
		(rf_gap_pos_start, rf_gap_pos_end) = (rf_gap_pos_start+z_gap_min, rf_gap_pos_end+z_gap_max)
		node_pos_start = rf_gap_pos_end
		(node_pos_tmp,node_pos_end) =  node_pos_dict[nodes[node_ind_end]]
		(local_quads,zero_length_nodes) = Get_quads_zeroLengthNodes_in_range(accSeq,node_ind_start,node_ind_end)
		if(abs(node_pos_start - node_pos_end) > drift_length_tolerance):
			#print "debug after n_quads=",len(local_quads)," n_markers=",len(zero_length_nodes)
			pos_start = node_pos_start - accSeq_z_start
			pos_end = node_pos_end - accSeq_z_start
			ovrlp_count = 0
			for ind in range(len(zero_length_nodes)+1):
				name = rf_gap.getName()+":After:"+str(ovrlp_count+1)+":OVRLPQ"
				node = OverlappingQuadsNode(name)
				if(ind == len(zero_length_nodes)):
					pos_end = node_pos_end - accSeq_z_start
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
				for quad in local_quads:
					node.addQuad(quad,enge_func_quad_dict[quad],quad.getPosition()-pos_start)
				node.setZ_Step(z_step)
				nParts = int(length/z_step)+1
				node.setnParts(nParts)					
				new_nodes.append(node)
				ovrlp_count += 1
				if(ind < len(zero_length_nodes)):
					new_nodes.append(zero_length_nodes[ind])
				pos_start = pos_end
		else:
			new_nodes += zero_length_nodes
		#--------------------------------------------------------------------------
		#---- 3st STEP - nodes of AxisField_and_Quad_RF_Gap type - new RF gaps
		#---- rf_gap_to_cavs_dict[rf_gap] = cav
		#---- new_rf_gaps_arr_in_cav_dict[cav] = [new_AxisField_and_Quad_RF_Gap_0,..]		
		n_rf_gaps = len(rf_gap_ind_up_down_arr)
		for local_gap_ind in range(n_rf_gaps):
			[rf_gap,gap_ind,drift_down_ind,drift_up_ind] = rf_gap_ind_up_down_arr[local_gap_ind]
			cav = rf_gap_to_cavs_dict[rf_gap]
			new_rf_gaps_arr = new_rf_gaps_arr_in_cav_dict[cav]
			axisFieldRF_Gap = af_rf_gap_dict[rf_gap]
			(local_quads,zero_length_nodes) = Get_quads_zeroLengthNodes_in_range(accSeq,drift_down_ind,drift_up_ind)
			local_quads = Find_Neihbor_Quads(accLattice,accSeq,rf_gap,gap_ind,af_rf_gap_dict,enge_func_quad_dict)
			(rf_gap_pos_start, rf_gap_pos_end) = node_pos_dict[rf_gap]
			rf_gap_pos_zero = rf_gap_pos_start - accSeq_z_start
			(z_gap_min,z_gap_max) = axisFieldRF_Gap.getZ_Min_Max()
			(rf_gap_pos_start, rf_gap_pos_end) = (rf_gap_pos_start+z_gap_min, rf_gap_pos_end+z_gap_max)
			pos_start = rf_gap_pos_start - accSeq_z_start
			pos_end = rf_gap_pos_end - accSeq_z_start
			#--------------------------------------------------------
			ovrlp_count = 0
			for ind in range(len(zero_length_nodes)+1):
				name = rf_gap.getName()+":GAP&Q:"+str(ovrlp_count+1)
				axisField_and_Quad_RF_Gap = AxisField_and_Quad_RF_Gap(axisFieldRF_Gap)
				new_rf_gaps_arr.append(axisField_and_Quad_RF_Gap)
				axisField_and_Quad_RF_Gap.setName(name)
				if(ind == len(zero_length_nodes)):
					pos_end = rf_gap_pos_end - accSeq_z_start
				else:
					pos_end = zero_length_nodes[ind].getPosition()
				length = pos_end - pos_start
				if(abs(length) < drift_length_tolerance):
					if(ind < len(zero_length_nodes)):
						new_nodes.append(zero_length_nodes[ind])
					pos_start = pos_end
					continue	
				pos = (pos_end + pos_start)/2
				z_min = pos_start - rf_gap_pos_zero
				z_max = pos_end - rf_gap_pos_zero
				axisField_and_Quad_RF_Gap.setZ_Min_Max(z_min,z_max)
				axisField_and_Quad_RF_Gap.setPosition(pos)
				for quad in local_quads:
					axisField_and_Quad_RF_Gap.addQuad(quad,enge_func_quad_dict[quad],quad.getPosition() - rf_gap_pos_zero)
				axisField_and_Quad_RF_Gap.setZ_Step(z_step)
				new_nodes.append(axisField_and_Quad_RF_Gap)
				ovrlp_count += 1
				if(ind < len(zero_length_nodes)):
					new_nodes.append(zero_length_nodes[ind])
				pos_start = pos_end
		#------------------------------------------------------------------------------------------------
		#---- sometimes there will be duplicate nodes: if RF filed covers the beginning or end of accSeq
		#---- we will remove the duplicate nodes
		new_nodes_tmp = []
		new_nodes_dict_tmp = {}
		for node in new_nodes:
			if(not new_nodes_dict_tmp.has_key(node)):
				new_nodes_tmp.append(node)
				new_nodes_dict_tmp[node] = None
		new_nodes = new_nodes_tmp
		#--------------------------------------------------------------------------
		#---- 4st STEP - nodes between two RF gaps
		for local_gap_ind in range(n_rf_gaps-1):
			[rf_gap0,gap0_ind,drift_down0_ind,drift_up0_ind] = rf_gap_ind_up_down_arr[local_gap_ind]
			[rf_gap1,gap1_ind,drift_down1_ind,drift_up1_ind] = rf_gap_ind_up_down_arr[local_gap_ind+1]
			(rf_gap0_pos_start, rf_gap0_pos_end) = node_pos_dict[rf_gap0]
			(z_gap0_min,z_gap0_max) = af_rf_gap_dict[rf_gap0].getZ_Min_Max()
			(rf_gap0_pos_start, rf_gap0_pos_end) = (rf_gap0_pos_start+z_gap0_min, rf_gap0_pos_end+z_gap0_max)
			(rf_gap1_pos_start, rf_gap1_pos_end) = node_pos_dict[rf_gap1]
			(z_gap1_min,z_gap1_max) = af_rf_gap_dict[rf_gap1].getZ_Min_Max()
			(rf_gap1_pos_start, rf_gap1_pos_end) = (rf_gap1_pos_start+z_gap1_min, rf_gap1_pos_end+z_gap1_max)
			(local_quads,zero_length_nodes) = Get_quads_zeroLengthNodes_in_range(accSeq,drift_up0_ind,drift_down1_ind)
			node_pos_start = rf_gap0_pos_end - accSeq_z_start
			node_pos_end = rf_gap1_pos_start - accSeq_z_start			
			if(abs(node_pos_start - node_pos_end) > drift_length_tolerance):
				pos_start = node_pos_start
				pos_end = node_pos_end
				ovrlp_count = 0
				for ind in range(len(zero_length_nodes)+1):
					name = rf_gap0.getName()+":After:"+str(ovrlp_count+1)+":OVRLPQ"
					node = OverlappingQuadsNode(name)
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
					for quad in local_quads:
						node.addQuad(quad,enge_func_quad_dict[quad],quad.getPosition()-pos_start)
					node.setZ_Step(z_step)
					nParts = int(length/z_step)+1
					node.setnParts(nParts)					
					new_nodes.append(node)
					ovrlp_count += 1
					if(ind < len(zero_length_nodes)):
						new_nodes.append(zero_length_nodes[ind])
					pos_start = pos_end
			else:
				new_nodes += zero_length_nodes
		#--------------------------------------------------------------------------
		new_nodes = sorted(new_nodes, key = lambda x: x.getPosition(), reverse = False)
		#--------------------------------------------------------------------------
		#---- rf_gap_to_cavs_dict[rf_gap] = cav
		#---- new_rf_gaps_arr_in_cav_dict[cav] = [new_AxisField_and_Quad_RF_Gap_0,..]
		#---- the new istances of the AxisField_and_Quad_RF_Gap has setAsFirstRFGap(False) by default 
		for cav in cavs:
			new_rf_gaps_arr = new_rf_gaps_arr_in_cav_dict[cav]
			cav.removeAllGapNodes()
			for axisField_and_Quad_RF_Gap in new_rf_gaps_arr:
				cav.addRF_GapNode(axisField_and_Quad_RF_Gap)
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
		if(not isinstance(node,BaseRF_Gap)):
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
	
def Find_Neihbor_Quads(accLattice,accSeq,rf_gap,gap_ind,af_rf_gap_dict,enge_func_quad_dict):
	"""
	This function will find the quads with the fields overlapping RF gap fields.
	It returns the quads that should be taken into account.
	gap_ind is an index of rf_gap in the original accLattice.getNodes().
	"""
	nodes = accSeq.getNodes()
	node_pos_dict = accLattice.getNodePositionsDict()
	(z_min,z_max) = af_rf_gap_dict[rf_gap].getZ_Min_Max()
	(gap_pos_start,gap_pos_end) = node_pos_dict[rf_gap]
	(gap_pos_start,gap_pos_end) = (gap_pos_start+z_min,gap_pos_end+z_max)	
	quads = []
	#---- let's go down
	node_ind = gap_ind
	while(node_ind >= 0):
		node = nodes[node_ind]
		if(isinstance(node,Quad)):
			(node_pos_start,node_pos_end) = node_pos_dict[node]
			(z_min_node,z_max_node) = enge_func_quad_dict[node].getLimitsZ()
			(node_pos_start,node_pos_end) = ((node_pos_start+node_pos_end)/2 + z_min_node,(node_pos_start+node_pos_end)/2 + z_max_node)
			delta = node_pos_end - gap_pos_start
			if(delta > 0.):
				quads.append(node)
			else:
				break
		node_ind -= 1
	#---- let's go up
	node_ind = gap_ind
	while(node_ind <= len(nodes) -1):
		node = nodes[node_ind]
		if(isinstance(node,Quad)):
			(node_pos_start,node_pos_end) = node_pos_dict[node]
			(z_min_node,z_max_node) = enge_func_quad_dict[node].getLimitsZ()
			(node_pos_start,node_pos_end) = ((node_pos_start+node_pos_end)/2 + z_min_node,(node_pos_start+node_pos_end)/2 + z_max_node)
			delta = node_pos_start - gap_pos_end
			if(delta < 0.):
				quads.append(node)
			else:
				break
		node_ind += 1
	quads = sorted(quads, key = lambda x: x.getPosition(), reverse = False)
	return quads