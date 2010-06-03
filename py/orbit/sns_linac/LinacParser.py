"""
The SimplifiedLinacParser reads the input xml data file and creates LinacStructureTree
that consists of sequences (class LinacStructureSeq) which in turn consist of nodes
(class LinacStuctureNode). These objects have limited functionality. They are just 
containers of dictionaries. Later the LinacStructureTree will be used to build the 
Accelerator Lattice instance.
"""

import sys,os

#import python XML DOM parser
import xml.dom.minidom

# import pyORBIT Python utilities classes for objects with names, types, and dictionary parameters
from orbit.utils import orbitFinalize
from orbit.utils import NamedObject, TypedObject, ParamsDictObject

class LinacStuctureNode(NamedObject, TypedObject, ParamsDictObject):
	"""
	The node that keeps the information from xml file.
	Parameters includes position length and others type dependent.
	"""
	def __init__(self,name = "node"):
		NamedObject.__init__(self,name)
		TypedObject.__init__(self,type_in = "none")
		ParamsDictObject.__init__(self)	
		self.setParam("length",0.)
		
	def getLength(self):
		"""
		Returns the total length of the node[m].
		"""		
		return self.getParam("length")
		
	def setLength(self, length):
		"""
		Sets the total length of the node [m].
		"""		
		return self.setParam("length",length)			

		
class LinacStructureSeq(LinacStuctureNode):
	"""	
	The linac sequence. It includes the LinacStuctureNodes. It has the length parameter
	"""
	def __init__(self, name = "None"):
		LinacStuctureNode.__init__(self,name)
		self.setType("sequence")
		self.nodes = []
		
	def addNode(self,node):
		if(not isinstance(node,LinacStuctureNode)):
			msg = "LinacStructureSeq: cannot add node. It is not LinacStuctureNode instance!"
			msg = msg + os.linesep
			msg = msg + "========================================="
			msg = msg + os.linesep
			orbitFinalize(msg)	
		self.nodes.append(node)

	def getNodes(self):
		"""
		Returns the array with linac nodes. It is a reference to the inner array, so
		user can modify it on his/her own risk.
		"""		
		return self.nodes

class LinacStructureTree(NamedObject, TypedObject, ParamsDictObject):
	"""	
	The linac lattice. It includes set of LinacStructureSeq.
	"""
	def __init__(self, name = "None"):
		NamedObject.__init__(self,name)
		TypedObject.__init__(self,type_in = "linac")
		ParamsDictObject.__init__(self)		
		self.seqs = []
		self.length = 0.
		
	def addSeq(self,node):
		if(not isinstance(node,LinacStructureSeq)):
			msg = "LinacStructureTree: cannot add node. It is not LinacStructureSeq instance!"
			msg = msg + os.linesep
			msg = msg + "========================================="
			msg = msg + os.linesep
			orbitFinalize(msg)	
		self.seqs.append(node)
		self.length = self.length + node.getParam("length")

	def getSeqs(self):
		"""
		Returns the array with sequences. It is a reference to the inner array, so
		user can modify it on his/her own risk.
		"""
		return self.seqs

	def getLength(self):
		"""
		Returns the total length of the linac [m].
		"""
		return self.length

class SimplifiedLinacParser:
	"""
	This is a parser for simplified XML file with a linac structure.
	The linac structure has sequences with elements. RF gaps has the name of a 
	RF cavity as an envelop structure.
	"""
	def __init__(self,xml_file_name):
		self.dom_doc = xml.dom.minidom.parse("sns_linac.xml")
		if(len(self.dom_doc.childNodes) != 1):
			msg = "SimplifiedLinacParser: input xml file has a wrong structure!"
			msg = msg + os.linesep
			msg = msg + "File: " + xml_file_name
			msg = msg + os.linesep
			msg = msg + "========================================="
			msg = msg + os.linesep
			orbitFinalize(msg)		
		self.domLinac = self.dom_doc.childNodes[0]
		self.linacTree = LinacStructureTree(name = self.domLinac.localName)		
		domSequences = self._stripDOMtoElements(self.domLinac)
		for domSeq in domSequences:
			linacSeq = LinacStructureSeq(name = domSeq.localName)
			print "debug name=",domSeq.localName
			seqParamDict = {}
			for i in range(domSeq.attributes.length):
				seqParamDict[domSeq.attributes.item(i).name] = domSeq.attributes.item(i).value
			self._transformDict(seqParamDict)
			linacSeq.setLength(seqParamDict["length"])
			linacSeq.setParam("rfFrequency",seqParamDict["rfFrequency"])
			linacSeq.setParam("bpmFrequency",seqParamDict["bpmFrequency"])
			domNodes = self._stripDOMtoElements(domSeq)
			for domNode in domNodes:
				nNodeParam = domNode.attributes.length
				paramDict = {}
				for i in range(nNodeParam):
					paramDict[domNode.attributes.item(i).name] = domNode.attributes.item(i).value
					#print "i=",i," name=",domNode.attributes.item(i).name," val=",domNode.attributes.item(i).value
				#---- add parameters from <parameters> child domNode
				domParameters = self._stripDOMtoElements(domNode)
				#print "domNode =",domNode.localName ," domParameters=",domParameters
				if(len(domParameters) > 1):
					msg = "SimplifiedLinacParser:more than 1 child in accNode!"
					msg = msg + os.linesep	
					orbitFinalize(msg)
				if(len(domParameters) == 1):
					domParameter = domParameters[0]
					for i in range(domParameter.attributes.length):
						paramDict[domParameter.attributes.item(i).name] = domParameter.attributes.item(i).value
				nodeName = paramDict["name"]
				#the name is not in the dictionary of parameters
				del paramDict["name"]			
				linacNode = LinacStuctureNode(name = nodeName)
				linacNode.setType(type_in = paramDict["type"])
				#the type is not in the dictionary of parameters
				del paramDict["type"]
				self._transformDict(paramDict)
				linacNode.setParamsDict(paramDict)
				linacSeq.addNode(linacNode)
			self.linacTree.addSeq(linacSeq)
			#print "name=",linacSeq.getName()," type=",linacSeq.getType()," length=",linacSeq.getLength()

	def _stripDOMtoElements(self,domNode):
		"""
		Removes all DOM Document componenets that are not DOM Elements and
		returns the array with dom elements.
		"""
		domChildren = []
		for child in domNode.childNodes:
			if(child.nodeType == child.ELEMENT_NODE):
				domChildren.append(child)
		return domChildren	

	def _transformDict(self, paramDict):
		"""
		It will replace {key:string} the float or integer number instead of string if it is possible or
		it will keep the orginal string.
		"""
		for key in paramDict.keys():
			s_in = paramDict[key]
			if(isinstance(s_in, basestring)):
				if(s_in.isdigit()):
					paramDict[key] = int(s_in)
				else:
					try:
						paramDict[key] = float(s_in)
					except ValueError:
						continue 
				


	def getLinacStructureTree(self):
		"""
		Returns the linac structure tree. It will be used to build the linac accelerator lattice.
		"""
		return self.linacTree


