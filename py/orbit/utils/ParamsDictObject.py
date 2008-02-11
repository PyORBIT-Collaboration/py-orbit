import os

from orbit.utils.orbitFinalize import orbitFinalize
from orbit.utils.NamedObject import NamedObject
from orbit.utils.TypedObject import TypedObject

class ParamsDictObject:
	"""
	An object that has a parameters dictionary.
	"""
	
	def __init__(self):
		"""
		Contsructor for an object that has a parameters dictionary.
		"""
		self.__paramsDict = {}

	def addParam(self, key, value):
		"""
		It sets the parameter of the node
		"""
		self.__paramsDict[key] = value

	def setParamsDict(self, params):
		"""
		It sets an external dictionary as a parameter
		dictionary for this element.
		"""
		self.__paramsDict = params

	def updateParamsDict(self, params):
		"""
		It updates the dictionary with an external dictionary data.
		"""
		self.__paramsDict.update(params)

	def getParam(self, key):
		"""
		It returns the parameters of the node
		"""
		if(not self.hasParam(key)):
			msg = "The object does not have a parameter for the key you requested!"
			msg = msg + os.linesep
			msg = msg + "method getParam(self, key)"
			msg = msg + os.linesep
			if(isinstance(self,NamedObject) == True):
				msg = msg + "Name of element = " + self.getName()
				msg = msg + os.linesep
			if(isinstance(self,TypedObject) == True):
				msg = msg + "Type of element = " + self.getType()
				msg = msg + os.linesep
			msg = msg + "key = " + str(key)
			orbitFinalize(msg)
		return self.__paramsDict[key]

	def getParamsDict(self):
		"""
		It returns the whole parameters dictionary.
		"""
		return self.__paramsDict

	def hasParam(self, key):
		"""
		It returns True if the node has a parameter for this key, or 
		returns False otherwise.
		"""
		return self.__paramsDict.has_key(key)


