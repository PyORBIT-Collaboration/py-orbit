class TypedObject:
	"""
	An object that has a type.
	"""
	
	def __init__(self, type_in = "no type" ):
		"""
		Contsructor for an object that has a type.
		"""
		self.__type = type_in

	def setType(self, type_in = "no type"):
		"""
		It sets the name.
		"""
		self.__type = type_in

	def getType(self):
		"""
		It returns the type.
		"""
		return self.__type
