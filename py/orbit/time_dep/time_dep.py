"""
Module. Includes classes for all time dependent lattice.
"""
import sys
import os
import math

from orbit.teapot import TEAPOT_Lattice
from orbit.parsers.mad_parser import MAD_Parser, MAD_LattLine
from orbit.lattice import AccNode, AccActionsContainer
from orbit.time_dep import waveform

class TIME_DEP_Lattice(TEAPOT_Lattice):
	"""
	The subclass of the TEAPOT_Lattice.
	TIME_DEP_Lattice has the ability to set time dependent parameters to the Lattice.
	Multi-turn track also available. 
	"""
	def __init__(self, name = "no name"):
		TEAPOT_Lattice.__init__(self,name)		
		self.__latticeDict = {}
		self.__TDNodeDict = {}
		self.__turns = 1
		
	def setLatticeOrder(self):
		"""
		Sets the time dependent lattice names to the lattice. 
		"""
		accNodes = self.getNodes()
		elemInLine = {}
		for i in range(len(accNodes)):
			elem = accNodes[i]			
			elemname = elem.getName()
			if(elemInLine.has_key(elemname)):
				elemInLine[elemname] += 1
			else:	elemInLine[elemname] = 1
			node = self.getNodes()[i]
			node.setParam("TPName",node.getName()+"_"+str(elemInLine[elemname]))
			#node.setParam("sequence",i+1)
			#print "debug node",node.getName(),node.getParamsDict()
     
	def setTimeDepNode(self, TPName, waveform):
		"""
		Sets the waveform function to the TP node before track.
		"""
        	flag = 0
                for node in self.getNodes():
                        if (TPName == node.getParam("TPName")):
                                flag = 1
                                node.setParam("waveform",waveform)
                                self.__TDNodeDict[TPName] = node
                if not flag:
                        print "The",TPName,"is not found."
                        sys.exit(1)
               
	def setTimeDepStrength(self, time):
		"""
		Set strength to the TP node while running. 
		"""
		NodeDict = self.__TDNodeDict		
		for i in NodeDict.keys():
			node = NodeDict[i]
			waveform = node.getParam("waveform")			
			waveform.calc(time)		
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

	def trackBunchTurns(self, bunch):
		"""
		It tracks the bunch through the lattice with multi-turn.
		"""
		turns = self.__turns
		#start
		for i in range(turns-1):			
			self.trackBunch(bunch)	
			syncPart = bunch.getSyncParticle()
			time = syncPart.time()
			self.setTimeDepStrength(time)
			print "debug trackBunchTurns time",time,"in",i,"turn"
		#getsublattice
		#sublattice.trackBunch(bunch)
		
	def setTurns(self, turns, startPosition = 0, endPosition = -1):
		"""
		Sets the turns and start end position before track.
		"""
		startNode = StartNode("start node")
		endNode = EndNode("end node")	
		self.addNode(startNode, startPosition)
		self.addNode(endNode, endPosition)
		self.__turns = turns	
		#print self.getNodes()

class StartNode(AccNode):
	def __init__(self, name = "no name"):
		AccNode.__init__(self,name)
		self.setType("start node")
		
	def track(self, paramsDict):
		bunch = paramsDict["bunch"]
		#bunch.getSyncParticle().time(0.)		
		
class EndNode(AccNode):
	def __init__(self, name = "no name"):
		AccNode.__init__(self,name)
		self.setType("end node")
		
	def track(self, paramsDict):
		pass
