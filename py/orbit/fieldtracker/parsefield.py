#!/usr/bin/env python

"""
This is not a parallel version! 
"""

import math
import random
import sys
from spacecharge import Grid3D

class FieldParser3D:
	""" 
	This class parses a file with 3D field and returns it as a 3D grid object
	"""
	def __init__(self):
			
	def parseFile(self, filename, xsize, ysize, zsize):
		
		infile = open(filename,"r")
		
		fieldgrid3D = Grid3D(sizeX, sizeY, sizeZ)
		
		for line in infile.readlines():
			print, 'parsing ', line
	
		return fieldgrid3D


		   
