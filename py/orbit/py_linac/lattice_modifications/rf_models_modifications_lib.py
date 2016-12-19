#!/usr/bin/env python

#--------------------------------------------------------------
# The functions for the lattice modifications to replace the base 
# RF gap nodes with 
# 1. The Axis Fields nodes. 
#    The initial base RF gap nodes have zero length, and the length
#    of new nodes is defined by the RF fields on the RF gap axis.
#    We also have to adjust these lengths to avoid overlapping. 
#    This correction will be on the level of microns, and it will 
#    not change the strength of the fields.
#    We assume that the non-zero field does not overlap any other 
#    elements of the lattice. So we have to modify the drifts only.
#--------------------------------------------------------------

import math
import sys
import os
import time

#---- MPI environment
import orbit_mpi
from orbit_mpi import mpi_comm

# import from orbit Python utilities
from orbit.utils import orbitFinalize

from orbit.py_linac.lattice import LinacApertureNode
from orbit.py_linac.lattice import Quad, Drift
from orbit.py_linac.lattice import BaseRF_Gap, AxisFieldRF_Gap

from orbit_utils import Function
from orbit_utils import SplineCH
from orbit_utils import GaussLegendreIntegrator


def Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location="",accSeq_Names = [],cavs_Names = []):
	"""
	Function will replace  BaseRF_Gap nodes by AxisFieldRF_Gap.
	It is assumed that AxisFieldRF_Gap nodes do not overlap any
	others nodes (except Drifts).
	The replacement will be performed only for specified sequences.
	If the cavities list is empty, all of them will be replaced!
	If you want to replace the nodes in a particular cavity please specify it!
	The dir_location is the location of the directory with the axis field
	files.
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
	drift_length_tolerance = 0.000000001
	for accSeq_Name in accSeq_Names:
		accSeq = accLattice.getSequence(accSeq_Name)
		if(accSeq == None):
			msg  = "The Replace_BaseRF_Gap_to_AxisField_Nodes Python function. "
			msg += "Cannot find the acc. sequence with this name in the lattice!"
			msg += os.linesep
			msg = msg + "accSeq name = " + 	accSeq_Name		
			msg = msg + os.linesep
			msg = msg + "lattice name = " + accLattice.getName()
			msg = msg + os.linesep
			orbitFinalize(msg)			
		#print "debug ======================================== start seq=",accSeq_Name," L=",accSeq.getLength()," n nodes=",len(accSeq.getNodes())		
		nodes = accSeq.getNodes()
		node_pos_dict = accLattice.getNodePositionsDict()
		#--- just for case: if the nodes are not in the right order
		nodes = sorted(nodes, key = lambda x: x.getPosition(), reverse = False)		
		cavs = accSeq.getRF_Cavities()
		if(len(cavs_Names) > 0):
			cavs_tmp = []
			for cav in cavs:
				if(cav.getName() in cavs_Names):
					cavs_tmp.append(cav)
				cavs = cavs_tmp
		#---- let's check that all rf gaps are BaseRF_Gap instances
		for cav in cavs:
			rf_gaps = cav.getRF_GapNodes()
			for rf_gap in rf_gaps:
				if(not isinstance(rf_gap,BaseRF_Gap)):
					msg  = "The Replace_BaseRF_Gap_to_AxisField_Nodes Python function. "
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
		for rf_gap in af_rf_gap_dict.keys():
			af_rf_gap_dict[rf_gap].setZ_Step(z_step)
		#---- check that all elements covered by the axis rf fields are drifts
		for [rf_gap,gap_ind,drift_down_ind,drift_up_ind] in rf_gap_ind_up_down_arr:
			#print "debug rf_gap=",rf_gap.getName()," gap_ind=",gap_ind," drift_down_ind=",drift_down_ind," drift_up_ind=",drift_up_ind
			ind_arr = []
			if(drift_down_ind < gap_ind):
				ind_arr += range(drift_down_ind,gap_ind)
			if(drift_up_ind > gap_ind):
				ind_arr += range(gap_ind+1,drift_up_ind+1)
			for node_ind in ind_arr:
				node = nodes[node_ind]
				if(not isinstance(node,Drift)):
					(gap_pos_start,gap_pos_end) = node_pos_dict[rf_gap]
					(z_min,z_max) = af_rf_gap_dict[rf_gap].getZ_Min_Max()
					(node_pos_start,node_pos_end) = node_pos_dict[node]
					msg  = "The Replace_BaseRF_Gap_to_AxisField_Nodes Python function. "
					msg += os.linesep
					msg += "The RF gap field covers the non-drift node! Stop!"
					msg += os.linesep
					msg = msg + "RF gap =" + rf_gap.getName()				
					msg = msg + os.linesep
					msg = msg + "Type of gap node= " + rf_gap.getType()
					msg = msg + os.linesep
					msg = msg + "Pos of the gap (s_start,s_stop)= " + str((gap_pos_start+z_min,gap_pos_start+z_max))
					msg = msg + os.linesep				
					msg = msg + "Near non drift element= " + node.getName()
					msg = msg + os.linesep	
					msg = msg + "Type of non drift element= " + node.getType()
					msg = msg + os.linesep	
					msg = msg + "Pos of non drift element (s_start,s_stop)= " + str((node_pos_start,node_pos_end))
					msg = msg + os.linesep						
					orbitFinalize(msg)
		#-----------------------------------------------------------------
		#---- let's change the length of the drifts that overlaps the ends 
		#---- of the axis rf fields. This will not change the positions of the 
		#---- drifts and rf_gaps in the node_pos_dict because the lattice is
		#---- not initialized yet.
		#dbg_chng_length_drift_arr = []
		for [rf_gap,gap_ind,drift_down_ind,drift_up_ind] in rf_gap_ind_up_down_arr:
			#print "debug modify rf gap=",rf_gap.getName()
			(gap_pos_start,gap_pos_end) = node_pos_dict[rf_gap]
			(z_min,z_max) = af_rf_gap_dict[rf_gap].getZ_Min_Max()
			#-----------------------------
			drift = nodes[drift_down_ind]
			"""
			if(len(dbg_chng_length_drift_arr) == 0):
				dbg_chng_length_drift_arr.append(drift)
			if(drift != dbg_chng_length_drift_arr[len(dbg_chng_length_drift_arr)-1]):
				dbg_chng_length_drift_arr.append(drift)
			"""
			(drift_pos_start,drift_pos_end) = node_pos_dict[drift]
			if(gap_pos_start + z_min >= drift_pos_start and gap_pos_start + z_min <= drift_pos_end):
				delta = drift_pos_end - (gap_pos_start + z_min)
				length = drift.getLength()
				drift.setLength(length - delta)
				drift.setPosition(drift.getPosition() - delta/2)
			#------------------------------
			drift = nodes[drift_up_ind]
			"""
			if(len(dbg_chng_length_drift_arr) > 0 and drift != dbg_chng_length_drift_arr[len(dbg_chng_length_drift_arr)-1]):
				dbg_chng_length_drift_arr.append(drift)
			"""
			(drift_pos_start,drift_pos_end) = node_pos_dict[drift]
			if(gap_pos_end + z_max >= drift_pos_start and gap_pos_end + z_max <= drift_pos_end):
				delta = (gap_pos_end + z_max) - drift_pos_start
				length = drift.getLength()
				drift.setLength(length - delta)	
				drift.setPosition(drift.getPosition() + delta/2)
		#------------------------------------------
		"""
		print "debug n_drifts =",len(dbg_chng_length_drift_arr)
		for drift in dbg_chng_length_drift_arr:
			length = drift.getLength()
			if(math.fabs(length) < rf_length_tolerance):
				print "debug drift = ",drift.getName()," L=",length
		"""
		#------------------------------------------
		#---- let's create the new set of nodes with the new type of rf gaps
		run_index = 0
		nodes_new = []
		for [rf_gap,gap_ind,drift_down_ind,drift_up_ind] in rf_gap_ind_up_down_arr:
			for node_ind in range(run_index,drift_down_ind):
				nodes_new.append(nodes[node_ind])
			if(nodes[drift_down_ind].getLength() > drift_length_tolerance):
				nodes_new.append(nodes[drift_down_ind])
			nodes_new.append(af_rf_gap_dict[rf_gap])
			if(nodes[drift_up_ind].getLength() > drift_length_tolerance):
				nodes_new.append(nodes[drift_up_ind])			
			run_index = drift_up_ind + 1
		for node_ind in range(run_index,len(nodes)):
			nodes_new.append(nodes[node_ind])
		#---- lets replace old fashion gaps in the cavity with the new ones
		for cav in cavs:
			rf_gaps = cav.getRF_GapNodes()
			cav.removeAllGapNodes()
			new_rf_gaps = []
			for rf_gap in rf_gaps:
				new_rf_gap = af_rf_gap_dict[rf_gap]
				cav.addRF_GapNode(new_rf_gap)
		#------------------------------------------
		"""
		for cav in cavs:
			rf_gaps = cav.getRF_GapNodes()
			print "debug cav=",cav.getName()," n_gaps=",len(rf_gaps )
		total_length = 0.
		for node_ind in range(len(nodes_new)):
			total_length +=  nodes_new[node_ind].getLength()
		print "total new L=",total_length, " n nodes=",len(nodes_new)
		"""
		#------------------------------------------
		#---- let's replace all nodes in the AccSeq by the new set
		accSeq.removeAllNodes()
		for node in nodes_new:
			accSeq.addNode(node)
		#---- new accSeq is ready. Now let's set up the lattice with the new nodes
	#---- new set of nodes for the lattice
	new_latt_nodes = []
	for accSeq in accLattice.getSequences():
		new_latt_nodes += accSeq.getNodes()
	accLattice.setNodes(new_latt_nodes)
	accLattice.initialize()
	#----------------------------------------------------------
	"""
	node_pos_dict = accLattice.getNodePositionsDict()
	#for node_ind in range(len(accLattice.getNodes())):
	for node_ind in range(73,109):
		node = accLattice.getNodes()[node_ind]
		#if(not isinstance(node,Quad)): continue
		(pos_start,pos_end) = node_pos_dict[node]
		print "debug ind=",node_ind," node=",node.getName()," (pos_start,pos_end)=",(pos_start,pos_end)
	"""
	#-----------------------------------------------------------
		
def Make_AxisFieldRF_Gaps_and_Find_Neihbor_Nodes(rf_length_tolerance,accLattice,accSeq,dir_location,cavs):
	"""
	It returns (af_rf_gap_dict,rf_gap_ind_up_down_arr).
	This function analyzes the nodes in the accSeq and creates a dictionary and 
	an array:
	af_rf_gap_dict[rf_gap] = AxisFieldRF_Gap(rf_gap)
	and
	rf_gap_ind_up_down_arr[[rf_gap,gap_ind,drift_down_ind,drift_up_ind],...]
	where rf_gap is a BaseRF_Gap instance, and indexes drift_down_ind and 
	drift_up_ind are the indexes covering the edges of the axis filed of the
	particular AxisFieldRF_Gap.
	"""
	rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
	nodes = accSeq.getNodes()
	node_pos_dict = accLattice.getNodePositionsDict()
	#--------------------------------------------------
	#---- let's create the new AxisFieldRF_Gap instances
	af_rf_gap_dict = {}
	for cav in cavs:
		rf_gaps = cav.getRF_GapNodes()
		for rf_gap in rf_gaps:
			af_rf_gap = AxisFieldRF_Gap(rf_gap)
			af_rf_gap.readAxisFieldFile(dir_location)
			af_rf_gap_dict[rf_gap] = af_rf_gap
	#--------------------------------------------------
	#---- Let's fix the length of the axis fields to avoid the fields overlaping
	for cav in cavs:
		rf_gaps = cav.getRF_GapNodes()
		for gap_ind in range(len(rf_gaps) - 1):
			rf_gap0 = rf_gaps[gap_ind]
			rf_gap1 = rf_gaps[gap_ind+1]
			(gap0_pos_start,gap0_pos_end) = node_pos_dict[rf_gap0]
			(gap1_pos_start,gap1_pos_end) = node_pos_dict[rf_gap1]
			gap0_pos_end += af_rf_gap_dict[rf_gap0].getZ_Min_Max()[1]
			gap1_pos_start += af_rf_gap_dict[rf_gap1].getZ_Min_Max()[0]
			delta_z = gap0_pos_end - gap1_pos_start
			if(math.fabs(delta_z) < rf_length_tolerance):
				(z_min,z_max) = af_rf_gap_dict[rf_gap0].getZ_Min_Max()
				z_max -= delta_z
				af_rf_gap_dict[rf_gap0].setZ_Min_Max(z_min,z_max)
				#print "debug gap0=",rf_gap0.getName()," gap1=",rf_gap1.getName()," delta=",delta_z
			else:
				if(delta_z > 0.):
					msg  = "The Replace_BaseRF_Gap_to_AxisField_Nodes Python function. "
					msg += os.linesep
					msg += "The RF gap field overlaps more than rf_length_tolerance[mm]= "+str(1000.*rf_length_tolerance)
					msg += os.linesep
					msg = msg + "RF gap 0 = " + rf_gap0.getName()				
					msg = msg + os.linesep
					msg = msg + "RF gap 1 = " + rf_gap1.getName()				
					msg = msg + os.linesep
					(z_min,z_max) = af_rf_gap_dict[rf_gap0].getZ_Min_Max()
					(pos_start,pos_stop) = (gap0_pos_start+z_min,gap0_pos_start+z_max) 
					msg = msg + "Gap 0 (pos_start,pos_stop)= " + str((pos_start,pos_stop))
					msg = msg + os.linesep
					(z_min,z_max) = af_rf_gap_dict[rf_gap1].getZ_Min_Max()
					(pos_start,pos_stop) = (gap1_pos_end+z_min,gap1_pos_end+z_max) 
					msg = msg + "Gap 1 (pos_start,pos_stop)= " + str((pos_start,pos_stop))
					msg = msg + os.linesep	
					msg = msg + "Gap 0 limits (z_min,z_max)= " + str(af_rf_gap_dict[rf_gap0].getZ_Min_Max())
					msg = msg + os.linesep	
					msg = msg + "Gap 1 limits (z_min,z_max)= " + str(af_rf_gap_dict[rf_gap1].getZ_Min_Max())
					msg = msg + os.linesep
					msg = msg + "Overlapping delta= " + str(delta_z)
					msg = msg + os.linesep
					orbitFinalize(msg)
	#--------------------------------------------------------------------------------------
	#---- array with [rf_gap, drift indexes for the gap and drifts before and after the gap]
	#---- Here we will go through all rf gaps and find indexes of the drifts before (down) 
	#---- and after (up). These drifts should be replaced with the shorter drifts.
	#---- The drifts that covered by the RF gap field completely should be removed.
	#---------------------------------------------------------------------------------------
	#---- to speed up indexing let's build rf gaps vs. index dictionary
	rf_gap_ind_dict = {}
	for node_ind in range(len(nodes)):
		node = nodes[node_ind]
		if(isinstance(node,BaseRF_Gap)):
			rf_gap_ind_dict[node] = node_ind
	#-------------------------------------
	rf_gap_ind_up_down_arr = []
	for cav in cavs:
		rf_gaps = cav.getRF_GapNodes()
		for rf_gap in rf_gaps:
			gap_ind = rf_gap_ind_dict[rf_gap]
			(gap_pos_start,gap_pos_end) = node_pos_dict[rf_gap]
			drift_down_ind = gap_ind
			drift_up_ind = gap_ind
			(z_min,z_max) = af_rf_gap_dict[rf_gap].getZ_Min_Max()
			#---- let's find the next upstream node covering the edge of the field
			drift_down_ind = gap_ind - 1
			node =  nodes[drift_down_ind]
			(node_pos_start,node_pos_end) = node_pos_dict[node]
			while(node_pos_end > gap_pos_start + z_min):
				drift_down_ind = drift_down_ind - 1
				#--- if it is the beginning of sequence - we are done!
				if(drift_down_ind < 0):
					if(gap_pos_start + z_min < node_pos_start):
						node = nodes[drift_down_ind+1]
						#---- by default gap_pos_start=gap_pos_end for rf gap with length=0
						(gap_pos_start,gap_pos_end) = node_pos_dict[rf_gap]
						(z_min,z_max) = af_rf_gap_dict[rf_gap].getZ_Min_Max()
						(gap_pos_start,gap_pos_end) = (gap_pos_start+z_min,gap_pos_end+z_max) 
						(pos_start,pos_end) = node_pos_dict[node]
						func = af_rf_gap_dict[rf_gap].getAxisFieldFunction()
						delta_cut = pos_start - gap_pos_start
						func_new = RenormalizeFunction(func,z_min+delta_cut,z_max)
						af_rf_gap_dict[rf_gap].setAxisFieldFunction(func_new)
						(z_min_new,z_max_new) = (func_new.getMinX(),func_new.getMaxX())
						af_rf_gap_dict[rf_gap].setZ_Min_Max(z_min_new,z_max_new)
						msg  = "debug =============== WARNING  START ================ RF Gap="+rf_gap.getName()	
						msg += os.linesep
						msg += "Inside the Replace_BaseRF_Gap_to_AxisField_Nodes Python function. "
						msg += os.linesep
						msg += "The RF gap field overlaps the first element in AccSequence."
						msg += os.linesep
						msg += "It means that the field goes outside the AccSequence."
						msg += os.linesep
						msg += "That is wrong! The field will be cut shorter and re-normalized!"
						msg += os.linesep
						msg += "node    = " + node.getName()			
						msg += os.linesep
						msg += "node (pos_start,pos_end)   = " + 	str((pos_start,pos_end))
						msg += os.linesep
						msg += "rf_gap (pos_start,pos_end)   = " + 	str((gap_pos_start,gap_pos_end))
						msg += os.linesep
						msg += "old rf gap (z_min,z_max) = " + str((z_min,z_max))
						msg += os.linesep
						msg += "new rf gap (z_min,z_max) = " + str((z_min_new,z_max_new))
						msg += os.linesep
						msg += "debug =============== WARNING  END ================"
						if(rank == 0): print msg
						break
					else:
						break
				#---------------------------------------------------------------------
				node =  nodes[drift_down_ind]
				(node_pos_start,node_pos_end) = node_pos_dict[node]
			drift_down_ind = drift_down_ind + 1
			#---------------------------------
			drift_up_ind = gap_ind + 1
			node =  nodes[drift_up_ind]
			(node_pos_start,node_pos_end) = node_pos_dict[node]
			while(node_pos_start < gap_pos_start + z_max):	
				drift_up_ind = drift_up_ind + 1
				#--- if it is the end of sequence - we are done!
				if(drift_up_ind > len(nodes) -1):
					if(gap_pos_start + z_max > node_pos_end):
						node = nodes[drift_up_ind-1]
						msg  = "The Replace_BaseRF_Gap_to_AxisField_Nodes Python function. "
						msg += os.linesep
						msg += "The RF gap field overlaps the last element in AccSequence."
						msg += os.linesep
						msg += "It means that the field goes outside the AccSequence."
						msg += os.linesep
						msg += "That is wrong! Stop! Please check the lattice!"
						msg += os.linesep
						msg += "RF gap  = " + rf_gap.getName()				
						msg += os.linesep
						msg += "node    = " + node.getName()			
						msg += os.linesep
						(gap_pos_start,gap_pos_end) = node_pos_dict[rf_gap]
						(z_min,z_max) = af_rf_gap_dict[rf_gap].getZ_Min_Max()
						(gap_pos_start,gap_pos_end) = (gap_pos_end+z_min,gap_pos_end+z_max) 
						(pos_start,pos_end) = node_pos_dict[node]	
						msg += "node (pos_start,pos_end)   = " + 	str((pos_start,pos_end))
						msg += os.linesep
						msg += "rf_gap (pos_start,pos_end)   = " + 	str((gap_pos_start,gap_pos_end))
						msg += os.linesep
						orbitFinalize(msg)
					else:
						break
				node =  nodes[drift_up_ind]
				(node_pos_start,node_pos_end) = node_pos_dict[node]
			drift_up_ind = drift_up_ind - 1
			rf_gap_ind_up_down_arr.append([rf_gap,gap_ind,drift_down_ind,drift_up_ind])
	#----------------------------------------------------------------------
	"""
	#---- Debug printing part of the code ---------------START-------------
	for node_ind in range(len(nodes)):
		node = nodes[node_ind]
		if(af_rf_gap_dict.has_key(node)):
			(pos_start,pos_stop) = node_pos_dict[node]
			pos_start += af_rf_gap_dict[node].getZ_Min_Max()[0]
			pos_stop += af_rf_gap_dict[node].getZ_Min_Max()[1]
			print "debug node_ind=",node_ind," node=",node.getName()," (pos_start,pos_end)=",(pos_start,pos_stop)
		else:
			print "debug node_ind=",node_ind," node=",node.getName()," (pos_start,pos_end)=",node_pos_dict[node]
	for [rf_gap,gap_ind,drift_down_ind,drift_up_ind] in rf_gap_ind_up_down_arr:
		print "debug gap=",rf_gap.getName()," gap_ind=",gap_ind," drift_down_ind=",drift_down_ind," drift_up_ind=",drift_up_ind
	#---- Debug printing part of the code ---------------STOP--------------
	"""
	#----------------------------------------------------------------------
	return (af_rf_gap_dict,rf_gap_ind_up_down_arr)
	
def RenormalizeFunction(func,z_min,z_max):
	"""
	It re-normalizes the Function in the new limits (z_min,z_max).
	We assume that region of the function definition will be cut not extended.
	"""
	spline = SplineCH()					
	spline.compile(func)	
	integrator = GaussLegendreIntegrator(500)
	integrator.setLimits(z_min,z_max)
	integral = integrator.integral(spline)
	n_points = func.getSize()
	step = (z_max - z_min)/(n_points-1)
	new_func = Function()
	for i in range(n_points):
		x = z_min + step*i
		y = spline.getY(x)/integral
		new_func.add(x,y)
	new_func.setConstStep(1)
	return new_func
	