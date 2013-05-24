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
		
		
	def dumb(self, i):
		return 3 + i
		
		
	def getLimits(self,filename):
		
		infile = open(filename,"r")
		# 		Resize the grid so more files can be processed
		numLines = 0
		xmin,xmax,ymin,ymax,zmin,zmax = 0,0,0,0,0,0
		for line in infile.readlines():
			splitline = line.split()
			value = map(float, splitline)
			
			# finding minimums and maximums of data to compute range
			xmax = max(xmax,value[0])
			ymax = max(ymax,value[1])
			zmax = max(zmax,value[2])
			xmin = min(xmin,value[0])
			ymin = min(ymin,value[1])
			zmin = min(zmin,value[2])
			
			numLines += 1
		
		xmax = max(xmax,0)
		ymax = max(ymax,0)
		zmax = max(zmax,0)
		xmin = min(xmin,0)
		ymin = min(ymin,0)
		zmin = min(zmin,0)
		
		print "Min and Max values: " , xmin, xmax, ymin, ymax, zmin, zmax, numLines, "\n" 

		limits = [xmin, xmax, ymin, ymax, zmin, zmax, numLines]
			
		return limits	
	#limits is a list of format [xmin,xmax,ymin,ymax,zmin,zmax]
	def getRange(self,limits):
		xmin = limits[0]
 		xmax = limits[1]
 		ymin = limits[2]
 		ymax = limits[3]
 		zmin = limits[4]
 		zmax = limits[5]
 		
 	
 		Xrange = xmax-xmin
		Yrange = ymax-ymin
		Zrange = zmax-zmin
		
		range = [Xrange,Yrange,Zrange]
		print "Range of X,Y,Z values in data: ", Xrange, " ", Yrange, " ", Zrange
		
		return range
	
	def getGridSize(self,range, step):
		
		for i in xrange(3):
			range[i] = int(range[i]*1/step)
		gridSize = [range[0],range[1], range[2]]
		print "Grid Size [x,y,z]: " , gridSize
		
		return gridSize
	
	def getCoordinates(self, gridSize, step,rawNumbers, limits):
		
		coordinates = [rawNumbers[0] ,rawNumbers[1],rawNumbers[2]]
		for i in xrange(len(coordinates)):
			coordinates[i] = coordinates[i]*(1/step)
			coordinates[i] = coordinates[i]+2*-limits[2*i]
		coordinates = map(int, coordinates)
		
		return coordinates
	
	
###############################################################################
# Parameters 
# filename: name of the text file to be processed
# xsize: size of the grid in the x diminsion
# ysize: size of the grid in the y diminsion
# zsize: size of the grid in the z diminsion
# All Grid sizes are user defined.
###############################################################################
	def parse(self, filename, xsize, ysize, zsize):
				
		limits = self.getLimits(filename)
		
		range = self.getRange(limits)
		step = 0.5
		
		gridSize = self.getGridSize(range, step)

 		numLines = limits[6]
 		
 		
 		
		print "Number of lines in the file: ",numLines , "\n"
		
		#for now we will say that the size of the grid encompasses all datapoints
		print "GridSize " , gridSize[0],gridSize[1],gridSize[2]
  		fieldgrid3DBx = Grid3D(gridSize[0]+1,gridSize[1]+1,gridSize[2]+1)
  		fieldgrid3DBy = Grid3D(gridSize[0]+1,gridSize[1]+1,gridSize[2]+1)
  		fieldgrid3DBz = Grid3D(gridSize[0]+1,gridSize[1]+1,gridSize[2]+1)
  		fieldgrid3DMag = Grid3D(gridSize[0]+1,gridSize[1]+1,gridSize[2]+1)
  		
  		setGridX = (limits[0],limits[1])
  		setGridY = (limits[2],limits[3])
  		setGridZ = (limits[4],limits[5])
  		
  		
	
		print
	
		##
		# Maps values from file to grid.
		##
		infile1 = open(filename,"r")
		i = 0;
		for line in infile1.readlines():
			splitLine = line.split()
			rawNumbers = map(float, splitLine)
			coordinates = self.getCoordinates(gridSize,step,rawNumbers, limits)
# 			print coordinates
 			fieldgrid3DBx.setValue(rawNumbers[3]/10000.0, coordinates[0], coordinates[1], coordinates[2])
 			fieldgrid3DBy.setValue(rawNumbers[4]/10000.0, coordinates[0], coordinates[1], coordinates[2])
 			fieldgrid3DBz.setValue(rawNumbers[5]/10000.0, coordinates[0], coordinates[1], coordinates[2])
 			
 			getMag = ((rawNumbers[3]**2.0+rawNumbers[4]**2.0+rawNumbers[5]**2.0)**0.5)/10000.0
 		 	
 		 	fieldgrid3DMag.setValue(getMag, coordinates[0], coordinates[1], coordinates[2])
#  		 	print fieldgrid3Dx.getValueOnGrid(coordinates[0], coordinates[1], coordinates[2])
 		
 		MagList = [fieldgrid3DBx,fieldgrid3DBy,fieldgrid3DBz,fieldgrid3DMag]
 		
		return MagList
	