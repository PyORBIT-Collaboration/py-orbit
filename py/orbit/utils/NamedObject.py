class NamedObject:
	"""
	An object that has a name.
	"""
	
	def __init__(self, name = "no name" ):
		"""
		Contsructor for an object that has a name.
		"""
		self.__name = name

	def setName(self, name = "no name"):
		"""
		It sets the name.
		"""
		self.__name = name

	def getName(self):
		"""
		It returns the name.
		"""
		return self.__name
