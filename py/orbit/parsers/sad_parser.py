import os
import sys
import re
import math

class _possibleElementType:
	""" This class keeps all possible element's types """

	def __init__(self):
		self.__names_type = []
		self.__names_type.append("MARK")
		self.__names_type.append("MONI")
		self.__names_type.append("DRIFT")
		self.__names_type.append("APERT")
		self.__names_type.append("BEND")
		self.__names_type.append("SEXT")
		self.__names_type.append("QUAD")
		self.__names_type.append("MULT")
		self.__names_type.append("CAVI")

	def __del__(self):
		del self.__names_type

	def checkType(self,type_in):
		type = type_in.upper()
		if self.__names_type.count(type) == 0:
			print "Error of creating lattice element."
			print "There can not be an element with type:",type
			print "Stop."
			sys.exit (0)
		return type

#===============================================================

class SAD_LattElement:
	""" An Arbitrary Element in the Lattice """
	_typeChecker = _possibleElementType()

	def __init__(self,name,Tname):
		""" Create instance with name, typeName and type """
		self.__name = name
		self.__type = self._typeChecker.checkType(Tname)
		self.__par = {}

	def __del__(self):
		del self.__par

	def getName(self):
		""" Returns name of the element """
		return self.__name

	def getType(self):
		""" Returns type of the element """
		return self.__type

	def setType(self,tp):
		""" Sets the type of the element without checking """
		self.__type = tp

	def addParameter(self,nameOfPar,parVal):
		self.__par[nameOfPar] = parVal

	def getParameter(self,nameOfPar):
		if self.__par.has_key(nameOfPar) == 0:
			print "class SAD_LattElement, method getParameter"
			print "The name of Element =",self.__name
			print "The type of Element =",self.__type
			print "The Element's key-val =",self.__par
			print "This Element does not have Parameter=", nameOfPar
			print "Stop."
			sys.exit (0)
		return self.__par[nameOfPar]

	def hasParameter(self,nameOfPar):
		return self.__par.has_key(nameOfPar)

	def getParameters(self):
		""" Returns the dictionary with (key=name, val) pairs """
		return self.__par

	def getElements(self):
		""" Returns list of elements (only one here) """
		elements = []
		elements.append(self)
		return elements

#====================================================================

class SAD_LattLine:
	"""
	An Arbitrary Line in the Lattice.
	It includes SAD Lattice Lines and Elements
	"""

	def __init__(self,name):
		""" Create instance with list of lines or elements """
		self.__name = name
		self.__items = []
		self.__signDict = {}

	def __del__(self):
		del self.__name
		del self.__items

	def getType(self):
		return "LINE"

	def getName(self):
		""" Returns name of the line """
		return self.__name

	def getDirection(self,item):
		""" Returns direction for particular Item. """
		return self.__signDict[item.getName()]

	def addItem(self,item,sign = +1):
		""" Adds a line or element to this line with cetain direction."""
		self.__items.append(item)
		self.__signDict[item.getName()] = sign

	def getItems(self):
		""" Returns list with elements and lines """
		return self.__items

	def getLinesDict(self):
		""" Returns the dictionary with all lattice lines inside, recursive. """
		dict = {}
		for item in self.__items:
			if(item.getType() == self.getType()):
				dict[item.getName()] = item
		return dict

	def getElements(self):
		""" Returns list of elements """
		elements = []
		for item in self.__items:
			sign = self.__signDict[item.getName()]
			elems = item.getElements()
			if(sign == -1): elems.reverse()
			for el in elems:
				elements.append(el)
		return elements

class _SAD_String:
	"""
	The line of a SAD file. It could be one of three types:
	variable, accelerator line, or accelerator element.
	"""

	def __init__(self):
		"""
		Constructor. Creates the SAD file line class instance.
		"""
		self.__line = ""
		self.__type = None

	def setLine(self,line):
		self.__line = line

	def getLine(self):
		return self.__line

	def setType(self,type):
		self.__type = type

	def getType(self):
		return self.__type

class _variable:
	"""
	The SAD variable class. It keeps initial string (line) from SAD file.
	"""

	def __init__(self):
		self._name = None
		self._expression = ""
		self._value = None

	def getType():
		"""
		Method. It is static method of this class.
		It returns the name of the type.
		"""
		return "variable"

	getType = staticmethod(getType)


	def getValue(self):
		"""
		Method. It returns the numerical value of this variable.
		"""
		return self._value

	def setValue(self, val):
		"""
		Method. It sets the numerical value of this variable.
		"""
		self._value = val

	def getName(self):
		return self._name

	def getExpression(self):
		return self._expression

	def parseLine(self,line_init):
		"""
		Method. It does the first parsing of the initial string.
		"""
		#divide string onto two parts: name of value and value
		names_vals = re.compile(r'\A\s*([^=\s]+)\s*={1}(.*)').findall(line_init)
		if(len(names_vals) != 1 and len(names_vals[0]) != 2 ):
			print "The SAD file line=",line_init
			print "It cannot be parsed. It is not variable declaration"
			print "Stop."
			sys.exit (0)
		self._name = names_vals[0][0]
		#replace all math cos,sin, etc by math.cos, math.sin, GEV, DEG etc
		self._expression = StringFunctions.replaceMath(names_vals[0][1])

class _element:
	"""
	The SAD element class. It also keeps initial string (line) from SAD file.
	"""

	def __init__(self):
		self._name = None
		self._elementType = None
		self._expression = ""
		self._parameters = {}
		self._numParameters = {}

	def getType():
		"""
		Method. It is static method of this class.
		It returns the name of the type.
		"""
		return "element"

	getType = staticmethod(getType)

	def getName(self):
		return self._name

	def getExpression(self):
		return self._expression

	def getElementType(self):
		return self._elementType

	def getParameters(self):
		"""
		Method. It returns the dictionary with key:string pairs
		"""
		return self._parameters

	def getNumParameters(self):
		"""
		Method. It returns the dictionary with {key:(numeric value)} pairs
		"""
		return self._numParameters

	def parseLine(self,line_init):
		"""
		Method. It does the first parsing of the initial string.
		"""
		type_name = re.compile(r'\A\s*(\S+)\s+([^=\s]+)\s*={1,1}').findall(line_init)[0]
		type = type_name[0]
		name = type_name[1]
		self._name = name
		self._elementType = type
		self._expression = line_init
		#==================================================
		#let's get key-val pairs
		#==================================================
		indStart = line_init.find(name)
		line_init = line_init[indStart:]
		(ind0,ind1) = StringFunctions.findBracketsPos(line_init)
		elemLine = line_init[ind0+1:ind1]
		keys = re.compile(r'\s*([^=\s]+)\s*={1,1}').findall(elemLine)
		for i in xrange(len(keys)):
			i_start = elemLine.find(keys[i])
			i_stop = len(elemLine)
			if(i != (len(keys) -1)):
				i_stop = elemLine.find(keys[i+1])
			vals = re.compile(r'\s*[^=\s]+\s*={1,1}(.*)').findall(elemLine[i_start:i_stop])
			if(len(vals) != 1):
				print "The SAD file line:",line_init
				print "The sub-string for parsing:",elemLine[i_start:i_stop]
				print "It cannot be parsed. Cannot parse the key ",keys[i]," value"
				print "Stop."
				sys.exit (0)
			self._parameters[keys[i].upper()] = vals[0]
			self._numParameters[keys[i].upper()] = vals[0]

class _accLine:
	"""
	The SAD accelerator class. It also keeps initial string (line) from SAD file.
	The component lines could have reverse order of element. We keep sign in
	compSign array.
	"""

	def __init__(self):
		self._name = None
		self._expression = ""
		self._components = []

	def getType():
		"""
		Method. It is static method of this class.
		It returns the name of the type.
		"""
		return "accLine"

	getType = staticmethod(getType)

	def getName(self):
		return self._name

	def getExpression(self):
		return self._expression

	def getComponents(self):
		"""
		Method. It returns the set with (components' names, sign)
		The sign define the order of elements from sub-line.
		The sign +- 1. The -1 means the revers order.
		"""
		return self._components

	def parseLine(self,line_init):
		"""
		Method. It does the first parsing of the initial string.
		Parses the SAD file line with lattice line definition.
		"""
		#the name was found before
		name = re.compile(r'\A\s*([^=\s]+)\s*={1,1}').findall(line_init)[0]
		self._name = name
		(ind0,ind1) = StringFunctions.findBracketsPos(line_init)
		lineLine = line_init[ind0+1:ind1]
		self._expression = lineLine
		#define names of components (lines or elements)
		line_names = re.compile(r'\s*(\S+)\s*').findall(lineLine)
		#=========================================
		#deal with the N*name expressions
		#=========================================
		line_names_new = []
		for it_in in line_names:
			sign = +1
			it = it_in
			if(it.find("-") == 0):
				sign = -1
				it = it_in[1:]
			n_name = re.compile(r'([\d]+?)\*{1}(.*)').findall(it)
			if len(n_name) > 0:
				n_rep = int(n_name[0][0])
				it_new = n_name[0][1]
				#print "debug n=",n_rep," name=",name," rep_name=",it_new
				for i in range(1,n_rep+1):
					line_names_new.append((it_new,sign))
			else:
				line_names_new.append((it,sign))
		self._components = line_names_new

class StringFunctions:
	"""
	This class defines the set of static string functions.
	"""

	def calculateString(self,str_in, localDict):
		"""
		Method. It returns a tuple (True,value) if
		the expression can be evaluated and (False,None) otherwise.
		"""
		try:
			val = eval(self.replaceMath(str_in),globals(),localDict)
			return (True, val)
		except:
			return (False, None)

	calculateString = classmethod(calculateString)

	def replaceMath(self,str_in):
		"""
		Method. It replaces math symbols
		to make them readable for python eval().
		"""
		#replace .e by .0e
		str_out = re.sub("\.e",".0e",str_in)
		str_out = str_out.strip()
		#check the math operatons
		str_out = re.sub("sin\(","math.sin(",str_out)
		str_out = re.sub("SIN\(","math.sin(",str_out)
		str_out = re.sub("cos\(","math.cos(",str_out)
		str_out = re.sub("COS\(","math.cos(",str_out)
		str_out = re.sub("tan\(","math.tan(",str_out)
		str_out = re.sub("TAN\(","math.tan(",str_out)
		str_out = re.sub("exp\(","math.exp(",str_out)
		str_out = re.sub("EXP\(","math.exp(",str_out)
		str_out = re.sub("log\(","math.log(",str_out)
		str_out = re.sub("LOG\(","math.log(",str_out)
		str_out = re.sub("acos\(","math.acos(",str_out)
		str_out = re.sub("ACOS\(","math.acos(",str_out)
		str_out = re.sub("asin\(","math.asin(",str_out)
		str_out = re.sub("ASIN\(","math.asin(",str_out)
		str_out = re.sub("atan\(","math.atan(",str_out)
		str_out = re.sub("ATAN\(","math.atan(",str_out)
		str_out = re.sub("sqrt\(","math.sqrt(",str_out)
		str_out = re.sub("SQRT\(","math.sqrt(",str_out)
		#in SAD file we can have GEV and DEG
		res = re.compile(r'(.*)\s+GEV\s*').findall(str_out)
		if(len(res) == 1):
			str_out = res[0]
		res = re.compile(r'(.*)\s+gev\s*').findall(str_out)
		if(len(res) == 1):
			str_out = res[0]
		res = re.compile(r'(.*)\s+DEG\s*').findall(str_out)
		if(len(res) == 1):
			str_out = "("+res[0]+")*math.pi/180."
		res = re.compile(r'(.*)\s+deg\s*').findall(str_out)
		if(len(res) == 1):
			str_out = "("+res[0]+")*math.pi/180."
		return str_out

	replaceMath = classmethod(replaceMath)

	def findLastBracketPos(self,line,nPos,nOpen):
		""" The method will find a last closing bracket. """
		for ip in xrange(nPos,len(line)):
			if(line[ip] == "("):
				nOpen = nOpen + 1
			if(line[ip] == ")"):
				nOpen = nOpen - 1
			if(nOpen == 0):
				return ip
		return -1

	findLastBracketPos = classmethod(findLastBracketPos)

	def findFirstBracketPos(self,line,nPos):
		""" The method will find a first openning bracket. """
		for ip in xrange(nPos,len(line)):
			if(line[ip] == "("):
				return ip
		return -1

	findFirstBracketPos = classmethod(findFirstBracketPos)

	def findBracketsPos(self,line):
		""" The method returns the tuple with ( and ) positions."""
		i0 = self.findFirstBracketPos(line,0)
		if(i0 == -1):
				return (-1,-1)
		i1 = self.findLastBracketPos(line,i0,0)
		if(i1 == 1):
				return (-1,-1)
		return (i0,i1)

	findBracketsPos = classmethod(findBracketsPos)

#====================================================================

class SAD_Parser:
	""" SAD parser """
	_typeChecker = _possibleElementType()

	def __init__(self):
		""" Create instance of the SAD_Parser class """
		self.__SAD_Strings =  []
		self.__accValues = []
		self.__accElements = []
		self.__accLines = []
		self.__unknownLines = []
		self.__SADFilePath = ""
		#old style lattice elements and lines
		self.__lattElems = []
		self.__lattLines = []

	def __del__(self):
		del self.__SAD_Strings

	def parse(self,SADfileName):
		#1-st stage read SAD file into the lines array
		self.__SADFilePath = os.path.dirname(SADfileName)
		fileName = os.path.basename(SADfileName)
		#the initilize can be recursive if there are nested SAD files
		self.__init__()
		self.initialize(fileName)

		#-----------------------------------------------------------
		#The SADLines has all lines
		#Let's create values, elements, and accelerator lines arrays
		#-----------------------------------------------------------
		for line in self.__SAD_Strings:
			if(line.getType() == _accLine.getType()):
				accLine = _accLine()
				accLine.parseLine(line.getLine())
				self.__accLines.append(accLine)
				continue
			if(line.getType() == _variable.getType()):
				var = _variable()
				var.parseLine(line.getLine())
				self.__accValues.append(var)
				continue
			if(line.getType() == _element.getType()):
				elem = _element()
				elem.parseLine(line.getLine())
				self.__accElements.append(elem)
				continue
			print "Error in parsing the SAD file."
			print "Cannot parse the line:",line.getLine()
			print "Stop."
			sys.exit (0)
		#print "debug size values=",len(self.__accValues)
		#print "debug size elements=",len(self.__accElements)
		#print "debug size accLines=",len(self.__accLines)
		#----------------------------------------------------------
		# Check if there is a redefinition of elements or variables
		#----------------------------------------------------------
		dict = {}
		for var in self.__accValues:
			if(dict.has_key(var.getName())):
				print "Warning the variable:",var.getName(), " was redefined. Are you are sure?"
			dict[var.getName()] = var
		dict = {}
		for elem in self.__accElements:
			if(dict.has_key(elem.getName())):
				print "Warning the element:",elem.getName(), " was redefined. Are you are sure?"
			dict[elem.getName()] = elem
		dict = {}
		for line in self.__accLines:
			if(dict.has_key(line.getName())):
				print "Warning the line:",line.getName(), " was redefined. Are you are sure?"
			dict[line.getName()] = line
		#-------------------------------------------------------
		#Then let's calculate numerical values.
		#They can be defined recursivelly, so we need iterations
		#-------------------------------------------------------
		accVarDict = {}
		for var in self.__accValues:
			accVarDict[var.getName()] = var
		localValDict = {}
		doNotStop = True
		while(doNotStop):
			doNotStop = False
			accVarDictInner = accVarDict.copy()
			#print "debug variables dict. size=",len(accVarDictInner)
			for name,var in accVarDictInner.iteritems():
				str_in = var.getExpression()
				res,val = StringFunctions.calculateString(str_in.lower(),localValDict)
				if(res):
					localValDict[name.lower()] = val
					var.setValue(val)
					del accVarDict[name]
				else:
					doNotStop = True
			if(len(accVarDictInner) == len(accVarDict)):
				print "=========== Unresolved Variables============"
				for name,var in accVarDictInner.iteritems():
					print "name=",name,"  str=",var.getExpression()
				print "=========== SAD File Problem ==============="
				print "=================STOP======================="
				sys.exit(1)
		"""
		#debug printing
		for var in self.__accValues:
			print "debug val_name=",var.getName()," val=",var.getValue()
		"""
		#-------------------------------------------
		# Now calculate all parameters in key,string_value
		# for accelerator elements
		#--------------------------------------------
		for accElem in self.__accElements:
			kvs = accElem.getParameters()
			kvNums = accElem.getNumParameters()
			for key,val in kvs.iteritems():
				val_out = None
				if val != None:
					res,val_out = StringFunctions.calculateString(val.lower(),localValDict)
					if(not res):
						print "=============SAD File problem ==============",
						print "Problem with acc. element:",accElem.getName()
						print "SAD file line:",accElem.getExpression()
						print "Parameter name:",key
						print "Can not calculate string:",val
						print "============ STOP =========================="
						sys.exit(1)
				kvNums[key] = val_out
		"""
		#debug printing
		i = 0
		for accElem in self.__accElements:
			i = i + 1
			print "i=",i,"type=",accElem.getElementType()," name=",accElem.getName(),
			for key,val in accElem.getNumParameters().iteritems():
				print " ",key,":",val," ",
			print " "
		sys.exit(1)
		"""
		#---------------------------------------------
		#Let's create all lattice elements (old style)
		#---------------------------------------------
		for accElem in self.__accElements:
			lattElem = SAD_LattElement(accElem.getName(),accElem.getElementType())
			self.__lattElems.append(lattElem)
			kvs = accElem.getNumParameters()
			for key,val in kvs.iteritems():
				lattElem.addParameter(key,val)
		#----------------------------------------------
		#Let's create lattice lines (old style)
		#----------------------------------------------
		lattElemDict = {}
		for elem in self.__lattElems:
			lattElemDict[elem.getName()] = elem
		accLineDict = {}
		for accLine in self.__accLines:
			accLineDict[accLine.getName()] = accLine
		lattLineDict = {}
		for accLine in self.__accLines:
			lattLine = SAD_LattLine(accLine.getName())
			self.__lattLines.append(lattLine)
			lattLineDict[accLine.getName()] = lattLine
		for lattLine in self.__lattLines:
			name = lattLine.getName()
			accLine = accLineDict[name]
			childs = accLine.getComponents()
			for (child,sign) in childs:
				if(lattElemDict.has_key(child) and accLineDict.has_key(child)):
					print "=========== SAD File Problem ==============="
					print "Accelerator line and element have the same name:",child
					print "=================STOP======================="
				if((not lattElemDict.has_key(child)) and (not accLineDict.has_key(child))):
					print "=========== SAD File Problem ==============="
					print "Can not find Accelerator line and element with name:",child
					print "=================STOP======================="
				if(lattElemDict.has_key(child)):
					lattLine.addItem(lattElemDict[child])
				else:
					lattLine.addItem(lattLineDict[child], sign)

	def initialize(self,SAD_file_name):
		""" This method will do the preliminary editing the input file."""
		fl = open(os.path.join(self.__SADFilePath, SAD_file_name))
		sad_strings = []
		str_local = ""
		for str_in in fl.readlines():
			#check if the line is a comment
			str0 = ""
			if(str_in.find("!") == 0):
				continue
			if(str_in.find("!") > 0):
				str0 = str_in[0:str_in.find("!")].strip()
			else:
				str0 = str_in.strip()
			if(str0.rfind(";") >= 0):
				str_local = str_local + " " + str0
				strs = re.compile(r'([^;]*);{1,1}').findall(str_local)
				for str1 in strs:
					str1 = str1.strip()
					if(len(str1) > 0):
						sad_strings.append(str1.expandtabs(1))
				str_local = ""
			else:
				str_local = str_local + " " + str0
		#--------------------------------------------
		"""
		#debug printing the input file after initial modifications
		i = 0
		for str in sad_strings:
			i = i + 1
			print "i=",i," str=",str
		"""
		fl.close()
		# the sad_strings list should be transformed to
		# the list of _SAD_String with types
		for str0 in sad_strings:
			#check if it is a line
			type_name = re.compile(r'\A\s*(\S*)\s*').findall(str0)
			if(len(type_name) == 1 and type_name[0].upper() == "LINE"):
				#find the first line name
				name = re.compile(r'\A\s*\S+\s+(\S+)\s*={1,1}').findall(str0)[0]
				#cut the LINE from the beggining of the string
				str_rest = str0[str0.find(name):]
				while(str_rest != None and len(str_rest) > 0):
					#loop through all lines NAME = (...)
					(ind0,ind1) = StringFunctions.findBracketsPos(str_rest)
					if(ind0 < 0 or ind1 < 0): break
					sadString = _SAD_String()
					sadString.setLine(str_rest[0:ind1+1])
					sadString.setType(_accLine.getType())
					self.__SAD_Strings.append(sadString)
					str_rest = str_rest[ind1+1:].strip()
				continue
			#check if it is element - str = TYPE NAME =
			type_name = re.compile(r'\A\s*(\S+)\s+\S+\s*={1,1}').findall(str0)
			if(len(type_name) == 1 and self._typeChecker.checkType(type_name[0]) != None):
				type0 = self._typeChecker.checkType(type_name[0])
				#find the first line name
				name = re.compile(r'\A\s*\S+\s+(\S+)\s*={1,1}').findall(str0)[0]
				#cut the LINE from the beggining of the string
				str_rest = str0[str0.find(name):]
				while(str_rest != None and len(str_rest) > 0):
					#loop through all lines NAME = (...)
					(ind0,ind1) = StringFunctions.findBracketsPos(str_rest)
					if(ind0 < 0 or ind1 < 0): break
					sadString = _SAD_String()
					sadString.setLine(type0+" "+str_rest[0:ind1+1])
					sadString.setType(_element.getType())
					self.__SAD_Strings.append(sadString)
					str_rest = str_rest[ind1+1:].strip()
				continue
			#check the valriables
			name = re.compile(r'\A\s*(\S+)\s*={1,1}.*').findall(str0)
			if(len(name) == 1):
				sadString = _SAD_String()
				sadString.setLine(str0)
				sadString.setType(_variable.getType())
				self.__SAD_Strings.append(sadString)
				continue
			self.__unknownLines.append(str0)
			continue
		#--------------------------------------------
		#SADLine list is ready for sorting according
		# the types
		#--------------------------------------------

	def getUnknownLines(self):
		"""
		Returns the list of unreccongized string lines in the file.
		"""
		return self.__unknownLines

	def getSAD_Lines(self):
		"""
		Method. It returns the list of the lattice lines
		that are defined in the SAD file.
		"""
		return self.__lattLines

	def getSAD_Elements(self):
		"""
		Method. It returns the list of the lattice elements
		that are defined in the SAD file.
		"""
		return self.__lattElems

	def getSAD_Variables(self):
		"""
		Method. It returns the list of the variables
		that are defined in the SAD file.
		"""
		return self.__accValues

	def getSAD_LinesDict(self):
		"""
		Method. It returns the dictionary of the lattice lines
		that are defined in the SAD file.
		"""
		dict = {}
		for lattLine in self.__lattLines:
			dict[lattLine.getName()] = lattLine
		return dict

	def getSAD_ElementsDict(self):
		"""
		Method. It returns the dictionary of the lattice elements
		that are defined in the SAD file.
		"""
		dict = {}
		for lattElem in self.__lattElems:
			dict[lattElem.getName()] = lattElem
		return dict

	def getSAD_VariablesDict(self):
		"""
		Method. It returns the dictionary of the variables
		that are defined in the SAD file.
		"""
		dict = {}
		for var in self.__accValues:
			dict[var.getName()] = var.getValue()
		return dict
