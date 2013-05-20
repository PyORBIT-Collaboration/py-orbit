import os
import sys
import re
import math
from spacecharge import Grid3D

class Field_Parser3D:
	""" 3D field parser """

	def __init__(self):
		""" Create instance of the Field_Parser3D class """
		self.__lines =  []
	
	def __del__(self):
		del self.__lines

	def parse(self, filename, xsize, ysize, zsize):
		
		infile = open(filename,"r")
		
		fieldgrid3D = Grid3D(xsize, ysize, zsize)
		
		for line in infile.readlines():
			print 'parsing ', line
		
		return fieldgrid3D
		