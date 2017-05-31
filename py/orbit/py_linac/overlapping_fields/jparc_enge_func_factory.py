#!/usr/bin/env python

#--------------------------------------------------------------
# This is a Enge Function Factory specific for the J-PARC. Some 
# Enge's function parameters are defined by the aperture and length,
# and others are defined by the field distribution formula from Trace3D
# documentation.
#--------------------------------------------------------------

import math
import sys
import os

from overlapping_quad_fields_lib import PMQ_Trace3D_Function
from overlapping_quad_fields_lib import EngeFunction
from overlapping_quad_fields_lib import SimpleQuadFieldFunc
	
def JPARC_EngeFunctionFactory(quad):
	"""
	It generates the Enge's Function for J-PARC quads.
	"""
	name = quad.getName()
	length_param = quad.getLength()		
	#---- general PMQ function described in Trace3D documentation
	if(quad.hasParam("radIn") and quad.hasParam("radOut")):
		radIn = quad.getParam("radIn")
		radOut = quad.getParam("radOut")
		cutoff_level = 0.01
		if(name == "LI_DTL1:DTQ01"): cutoff_level = 0.02
		func = PMQ_Trace3D_Function(length_param,radIn,radOut,cutoff_level)
		return func
	#----- general Enge's Function
	if(quad.hasParam("aperture")):
		acceptance_diameter_param = quad.getParam("aperture")
		cutoff_level = 0.01	
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
