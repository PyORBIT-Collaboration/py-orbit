#!/usr/bin/env python

#--------------------------------------------------------------
# Function to replace overlapping quads with the special nodes
# This function will work for forward and backward lattices
# This is an example for SNS MEBT and can be used as a template 
# for other accelerator or sequence.
#--------------------------------------------------------------

import math
import sys
import os


from overlapping_quad_fields_lib import EngeFunction
from overlapping_quad_fields_lib import OverlappingQuadsNode
from overlapping_quad_fields_lib import OverlappingQuadsController

def SNS_MEBT_OverlappingQuadsSubst(accLattice):
	quads = accLattice.getQuads()
	quad_ovrlp_names = ["MEBT_Mag:QH05","MEBT_Mag:QV06","MEBT_Mag:QH07"]
	quads_enge_fields_arr = [] 
	for quad in quads:
		if(not quad.getName() in quad_ovrlp_names): continue
		length_param = 0.066
		acceptance_diameter_param = 0.0363
		cutoff_level = 0.01	
		func = EngeFunction(length_param,acceptance_diameter_param,cutoff_level)	
		quads_enge_fields_arr.append([quad,func])
	
	ovrlp_quads1 = OverlappingQuadsController("MEBT:QUADS05t07")
	ovrlp_quads1.addQuadsAndFields(accLattice,quads_enge_fields_arr)
	ovrlp_quads1.setMaxPartsDistance(0.01)

	quad_ovrlp_names = ["MEBT_Mag:QH08","MEBT_Mag:QV09","MEBT_Mag:QH10"]
	quads_enge_fields_arr = []
	for quad in quads:
		if(not quad.getName() in quad_ovrlp_names): continue
		length_param = 0.066
		acceptance_diameter_param = 0.0363
		cutoff_level = 0.01	
		func = EngeFunction(length_param,acceptance_diameter_param,cutoff_level)	
		quads_enge_fields_arr.append([quad,func])
	
	ovrlp_quads2 = OverlappingQuadsController("MEBT:QUADS09t10")
	ovrlp_quads2.addQuadsAndFields(accLattice,quads_enge_fields_arr)
	ovrlp_quads2.setMaxPartsDistance(0.01)
	accLattice.initialize()
	return (ovrlp_quads1,ovrlp_quads2)