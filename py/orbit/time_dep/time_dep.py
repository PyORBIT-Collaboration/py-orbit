"""
Module. Includes classes for all time dependent lattice.
"""
import sys
import os
import math

from orbit.teapot import TEAPOT_Lattice
from orbit.parsers.mad_parser import MAD_Parser, MAD_LattLine
from orbit.lattice import AccNode
from orbit.time_dep import waveform

class TIME_DEP_Lattice(TEAPOT_Lattice):
	"""
	"""
	def __init__(self, name = "no name"):
		TEAPOT_Lattice.__init__(self,name)		
		self.__mad_file = None
		self.__mad_line = None
		self.__latticeDict = {}
		
	def readMAD(self,mad_file_name,linename):
		TEAPOT_Lattice.readMAD(self,mad_file_name,linename)
		self.__mad_file = mad_file_name
		self.__mad_line = linename
			
	def setLatticeOrder(self):		
		parser = MAD_Parser()
		parser.parse(self.__mad_file)
		accLines = parser.getMAD_LinesDict()	
		accMAD_Line = accLines[self.__mad_line]
		accMADElements = accMAD_Line.getElements()
		elemInLine = {}
		for i in range(len(accMADElements)):
			elem = accMADElements[i]			
			elemname = elem.getName()
			if(elemInLine.has_key(elemname)):
				elemInLine[elemname] += 1
			else:	elemInLine[elemname] = 1
			node = self.getNodes()[i]
			node.setParam("order",elemInLine[elemname])
			
	def getTimeDepNode(self, elem, order):
		flag = 0
		for node in self.getNodes():
			if (elem == node.getName()) and (order == node.getParam("order")):
				flag = 1
				return node
		if not flag:
			print "The",order,"th element",elem,"is not found."
			sys.exit(1)
		
	def setTimeDepNode(self, elem, order, waveform):
		node = self.getTimeDepNode(elem, order)	
		waveformType = waveform.getType()		
		if waveformType == "kicker waveform":
			if node.getType() == "kick teapot":
				self.setParam(node,"kx",waveform.getKx())
				self.setParam(node,"ky",waveform.getKy())
			else: print "No kicker waveform added. Please check node type."
		elif waveformType == "magnet waveform":
			strength = waveform.getStrength()
			if node.getType() == "multipole teapot":		
				self.setParam(node,"kls",strength)
			elif node.getType() == "quad teapot":
				self.setParam(node,"kls",strength)
				self.setParam(node,"kq",strength)
			elif node.getType() == "solenoid teapot":
				self.setParam(node,"B",strength)
			else: print "No magnet waveform added. Please check node type."  

	def setParam(self, node, Kparam, strength):
		if node.hasParam(Kparam):	
			paramval = node.getParam(Kparam)
			if Kparam == "kls":
				newparamval = []
				for i in range(len(paramval)):					
					newparamval.append(paramval[i]*strength)
				paramval = newparamval
			else:paramval = paramval*strength
			node.setParam(Kparam,paramval)
