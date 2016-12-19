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
from overlapping_quad_fields_lib import SimpleQuadFieldFunc

def SNS_MEBT_OverlappingQuadsSubst(accLattice):
	"""
	This function transform the MEBT central 6 quads to the 
	overlapping quads. It is SNS specific.
	"""
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
	
def SNS_EngeFunctionFactory(quad):
	"""
	It generates the Enge's Function for SNS quads. For some of the quads in the SNS
	lattice we know the Enge's parameters. So, it is a SNS specific function.
	"""
	name = quad.getName()
	if(name.find("MEBT") >= 0):
		names_042mm = ["MEBT_Mag:QH05","MEBT_Mag:QV06","MEBT_Mag:QH07"]
		names_042mm += ["MEBT_Mag:QH08","MEBT_Mag:QV09","MEBT_Mag:QH10"]
		if(name in names_042mm):
			length_param = 0.066
			acceptance_diameter_param = 0.0363
			cutoff_level = 0.001	
			func = EngeFunction(length_param,acceptance_diameter_param,cutoff_level)
			return func
		names_032mm = ["MEBT_Mag:QH01","MEBT_Mag:QV02","MEBT_Mag:QH03","MEBT_Mag:QV04"]
		names_032mm += ["MEBT_Mag:QV11","MEBT_Mag:QH12","MEBT_Mag:QV13","MEBT_Mag:QH14"]
		if(name in names_032mm):
			length_param = 0.061
			acceptance_diameter_param = 0.029
			cutoff_level = 0.001	
			func = EngeFunction(length_param,acceptance_diameter_param,cutoff_level)
			return func
	if(name.find("DTL") >= 0):
		if(name == "DTL_Mag:PMQH100"):
			func = SimpleQuadFieldFunc(quad)
			return func
		#----- the SNS DTL quads are all the same L= 35 mm D = 25 mm
		length_param = 0.0408
		acceptance_diameter_param = 0.0192
		cutoff_level = 0.001	
		func = EngeFunction(length_param,acceptance_diameter_param,cutoff_level)
		return func
	#----- general Enge's Function
	length_param = quad.getLength()
	if(quad.hasParam("aperture")):
		acceptance_diameter_param = quad.getParam("aperture")
		cutoff_level = 0.001	
		func = EngeFunction(length_param,acceptance_diameter_param,cutoff_level)
		return func
	else:
		msg  = "SNS_EngeFunctionFactory Python function. "
		msg += os.linesep
		msg += "Cannot create the EngeFunction for the quad!"
		msg += os.linesep
		msg = msg + "quad name = " + 	quad.getName()	
		msg = msg + os.linesep		
		msg = msg + "It does not have the aperture parameter!"
		msg = msg + os.linesep		
		orbitFinalize(msg)
	return None	
