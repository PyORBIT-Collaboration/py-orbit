class NamedObject:
	"""
	Class. An object that has a name.
	"""

	def __init__(self, name = "no name" ):
		"""
		Contsructor. Object that has a name.
		"""
		self.__name = name

	def setName(self, name = "no name"):
		"""
		Method. Sets the name.
		"""
		self.__name = name

	def getName(self):
		"""
		Method. Returns the name.
		"""
		return self.__name
