class TypedObject:
	"""
	Class. Object that has a type.
	"""

	def __init__(self, type_in = "no type" ):
		"""
		Constructor. Object that has a type.
		"""
		self.__type = type_in

	def setType(self, type_in = "no type"):
		"""
		Method. Sets the type.
		"""
		self.__type = type_in

	def getType(self):
		"""
		Method. Returns the type.
		"""
		return self.__type
