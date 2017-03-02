#!/usr/bin/env python

#--------------------------------------------------------------
# This is a Enge Function Factory specific for the SNS. Some 
# Enge's function parameters are found by fitting the measured or
# calculated field distributions. Others are generated with quad's
# length and beam pipe diameter. 
# These parameters for SNS with the specific quads' names.
# For your accelerator you should create your own factory.
#--------------------------------------------------------------

import math
import sys
import os

from overlapping_quad_fields_lib import EngeFunction
from overlapping_quad_fields_lib import SimpleQuadFieldFunc
	
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
		#---- for the DTL_Mag:PMQH100 we have to cut the field because it leaks to the MEBT section
		if(name == "DTL_Mag:PMQH100"):
			func = SimpleQuadFieldFunc(quad)
			return func
		#---- for the DTL_Mag:PMQV300 we have to cut the field because it leaks to the DTL2 section
		if(name == "DTL_Mag:PMQV300"):
			func = SimpleQuadFieldFunc(quad)
			return func			
		#---- for the DTL_Mag:PMQH500 we have to cut the field because it leaks to the DTL4 section
		if(name == "DTL_Mag:PMQH500"):
			func = SimpleQuadFieldFunc(quad)
			return func
		#----- the SNS DTL quads are all the same L= 35 mm D = 25 mm
		length_param = 0.0408
		acceptance_diameter_param = 0.0192
		cutoff_level = 0.001	
		func = EngeFunction(length_param,acceptance_diameter_param,cutoff_level)
		return func
	if(name.find("CCL") >= 0):
		#---- for the CCL_Mag:QH00 we have to cut the field because it leaks to the DTL6 section
		if(name == "CCL_Mag:QH00"):
			func = SimpleQuadFieldFunc(quad)
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
