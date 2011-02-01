import os
import sys
import re
import math

#===============================================================

class _possibleElementType:
	"""
	Class. Specifies all possible element types
	"""
	def __init__(self):
		"""
		Constructor. Creates list of element types.
		"""
		self.__names_type = []
		self.__names_type.append("drift")
		self.__names_type.append("sbend")
		self.__names_type.append("rbend")
		self.__names_type.append("quad")
		self.__names_type.append("quadrupole")
		self.__names_type.append("sextupole")
		self.__names_type.append("octupole")
		self.__names_type.append("multipole")
		self.__names_type.append("solenoid")
		self.__names_type.append("kicker")
		self.__names_type.append("hkicker")
		self.__names_type.append("vkicker")
		self.__names_type.append("hkick")
		self.__names_type.append("vkick")
		self.__names_type.append("rfcavity")
		self.__names_type.append("rcollimator")
		self.__names_type.append("marker")
		self.__names_type.append("monitor")

	def __del__(self):
		"""
		Method. Deletes element.
		"""
		del self.__names_type

	def checkType(self, name_in):
		"""
		Method. Confirms validity of element type.
		"""
		name = name_in.lower()
		if self.__names_type.count(name) == 0:
			print "Error creating lattice element:"
			print "There is no element with type: ", name
			print "Stop."
			sys.exit (0)
		return name_in

#===============================================================

class MAD_LattElement:
	"""
	Class. Represents an arbitrary element in the lattice
	"""
	_typeChecker = _possibleElementType()

	def __init__(self, name, Typename):
		"""
		Constructor. Creates element with name, type,
		and parameter dictionary.
		"""
		self.__name = name
		self.__type = self._typeChecker.checkType(Typename)
		self.__par = {}

	def __del__(self):
		"""
		Method. Deletes parameters.
		"""
		del self.__par

	def getName(self):
		"""
		Method. Returns name of element
		"""
		return self.__name

	def setType(self, tp):
		"""
		Method. Sets the type of element without checking.
		"""
		self.__type = tp

	def getType(self):
		"""
		Method. Returns type of element
		"""
		return self.__type

	def addParameter(self, nameOfPar, parVal):
		"""
		Method. Adds parameter and value to element.
		"""
		self.__par[nameOfPar] = parVal

	def getParameter(self, nameOfPar):
		"""
		Method. Returns name of parameter.
		"""
		if self.__par.has_key(nameOfPar) == 0:
			print "class MAD_LattElement, method getParameter"
			print "The name of Element = ", self.__name
			print "The type of Element = ", self.__type
			print "The Element's key-val = ", self.__par
			print "This Element does not have Parameter = ", nameOfPar
			print "Stop."
			sys.exit (0)
		return self.__par[nameOfPar]

	def getParameters(self):
		"""
		Method. Returns parameter dictionary.
		"""
		return self.__par

	def getElements(self):
		"""
		Method. Returns list of elements (only one here)
		"""
		elements = []
		elements.append(self)
		return elements

#====================================================================

class MAD_LattLine:
	"""
	Class. Represents an arbitrary line in the lattice.
	"""
	def __init__(self, name):
		"""
		Constructor. Creates line with list of lines and/or elements.
		"""
		self.__name = name
		self.__items = []
		self.__sign_arr = []
		
	def __del__(self):
		"""
		Method. Deletes instance with list of lines and/or elements.
		"""
		del self.__name
		del self.__items

	def getName(self):
		"""
		Method. Returns name of the line.
		"""
		return self.__name

	def getType(self):
		"""
		Method. Returns type: "LINE".
		"""
		return "LINE"

	def addItem(self, item, sign = +1):
		"""
		Method. Adds a line or element to this line.
		"""
		self.__items.append(item)
		self.__sign_arr.append(sign)

	def getItems(self):
		"""
		Method. Returns list with elements and lines.
		"""
		return self.__items

	def getLinesDict(self):
		"""
		Method. Returns a dictionary
		containing all lattice lines, recursive.
		"""
		dict = {}
		for item in self.__items:
			if(item.getType() == self.getType() or item.getType() == self.getType().lower()):
				dict[item.getName()] = item
		return dict

	def getElements(self):
		""" Returns list of elements """
		elements = []
		#print "debug ================= name line = ",self.getName()
		for i in range(len(self.__items)):
			item = self.__items[i]
			sign = self.__sign_arr[i]
			#print "           debug item=",item.getName()," sign=",sign
			elems = item.getElements()
			if(sign == -1): elems.reverse()
			for el in elems:
				elements.append(el)
		return elements

#====================================================================

class _madLine:
	"""
	Class: A line of a MAD file. It can be one of three types:
	variable, accelerator line, or accelerator element.
	"""
	def __init__(self):
		"""
		Constructor. Creates a blank MAD file line class instance.
		"""
		self.__line = ""
		self.__type = None

	def setLine(self, line):
		"""
		Method. Sets a MAD file line class instance.
		"""
		self.__line = line

	def getLine(self):
		"""
		Method. Returns a MAD file line class instance.
		"""
		return self.__line

	def setType(self, ltype):
		"""
		Method. Sets a MAD file line type.
		"""
		self.__type = ltype

	def getType(self):
		"""
		Method. Returns a MAD file line type.
		"""
		return self.__type

#====================================================================

class _variable:
	"""
	Class. Holds MAD variables.
	"""
	def __init__(self):
		"""
		Constructor. Creates empty MAD variable.
		"""
		self._name = None
		self._expression = ""
		self._value = None

	def getType():
		"""
		Method. Static method of this class.
		Returns the name of the type.
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

	def setName(self, name):
		self._name = name

	def getName(self):
		return self._name

	def setExpression(self,line):
		self._expression = line

	def getExpression(self):
		return self._expression

	def parseLine(self,line_init):
		"""
		Method. It does the first parsing of the initial string.
		"""
		#divide string onto two parts: name of value and value
		patt = re.compile(r'(?<=:=).*')
		s_val = re.findall(patt,line_init)
		if(len(s_val) > 0):
			patt = re.compile(r'.*(?=:=)')
			s_name  = re.findall(patt,line_init)
			s_val  = s_val[0]
			s_name = s_name[0]
		else:
			#deal with const defenition like AA = 1.2
			patt = re.compile(r'(?<==).*')
			s_val = re.findall(patt,line_init)
			patt = re.compile(r'.*(?==)')
			s_name  = re.findall(patt,line_init)
			s_val  = s_val[0]
			s_name = s_name[0]
		self.setName(s_name)
		self.setExpression(s_val)

#====================================================================

class _element:
	"""
	The MAD element class. It also keeps initial string (line) from MAD file.
	"""

	def __init__(self):
		self._name = None
		self._elementType = None
		self._parameters = {}
		self._numParameters = {}

	def getType():
		"""
		Method. It is static method of this class.
		It returns the name of the type.
		"""
		return "element"

	getType = staticmethod(getType)

	def setName(self, name):
		self._name = name

	def getName(self):
		return self._name

	def setElementType(self, etype):
		self._elementType = etype

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
		#divide string onto three parts: name,type and key-values of element
		patt = re.compile(r'[\w]+(?=:)')
		name = re.findall(patt,line_init)[0]
		patt = re.compile(r'(?<=:)[\w]+?(?=\s|,)')
		type = re.findall(patt,line_init+",")[0]
		self.setName(name)
		self.setElementType(type)
		patt = re.compile(r'(?<=\w),.*')
		s_key_vals = re.findall(patt,line_init)
		#==================================================
		#now we have name and type, let's get key-val pairs
		#==================================================
		if len(s_key_vals) < 1:
			return
		s_key_vals = s_key_vals[0] + ","
		patt = re.compile(r'(?<=,).+?(?=,)')
		s_key_vals = re.findall(patt,s_key_vals)
		for s_key_val in s_key_vals:
			patt = re.compile(r'[\w]+?(?==)')
			key =  re.findall(patt,s_key_val)
			if(len(key) == 0):
				#we have stand alone key, e.g. T1
				key = s_key_val
				val = None
			else:
				key = key[0]
				patt = re.compile(r'(?<==).*')
				val = re.findall(patt,s_key_val)[0]
			self._parameters[key] = val
			self._numParameters[key] = val


class _accLine:
	"""
	The MAD accelerator class. It also keeps initial string (line) from MAD file.
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

	def setName(self, name):
		self._name = name

	def getName(self):
		return self._name

	def setExpression(self,line):
		self._expression = line

	def getExpression(self):
		return self._expression

	def getComponents(self):
		"""
		Method. It returns the set with tuples (components name, sign)
		The sign define the order of elements from sub-line.
		The sign +- 1. The -1 means the revers order.
		"""
		return self._components

	def parseLine(self,line_init):
		"""
		Method. It does the first parsing of the initial string.
		Parses the MAD file line with lattice line definition.
		The sign is +- 1. The -1 means the revers order for the elements
		inside the component.
		"""
		#define name of new line
		patt = re.compile(r'[\w]+(?=:)')
		name = re.findall(patt,line_init)[0]
		self.setName(name)
		patt = re.compile(r'(?<==).*')
		s_def = re.findall(patt,line_init)[0]
		self.setExpression(s_def)
		patt = re.compile(r'(?<=\(|,).+?(?=,|\))')
		item_names = re.findall(patt,s_def)
		#=========================================
		#deal with the N*name expressions
		#=========================================
		item_names_new = []
		for it_in in item_names:
			sign = +1
			it = it_in
			if(it.find("-") == 0):
				sign = -1
				it = it_in[1:]			
			patt = re.compile(r'[\d]+?(?=\*)')
			n_rep = re.findall(patt,it)
			if len(n_rep) > 0:
				n_rep = int(n_rep[0])
				patt = re.compile(r'(?<=\*)[\w]+')
				s=re.findall(patt,it)
				s=s[0]
				for i in range(1,n_rep+1):
					item_names_new.append((s,sign))
			else:
				item_names_new.append((it,sign))
		self._components = item_names_new

class StringFunctions:
	"""
	This class defines the set of static string functions.
	"""

	def replaceElementKeys(self, str_in, elem, key, value):
		"""
		Method. It will replace elem[key] in the string expression of this variable.
		"""
		new_val = r'(' + str(value) + r')'
		s = elem+r'\['+key+r'\]'
		patt = re.compile(s)
		str_out = re.sub(patt,new_val,str_in)
		return str_out

	replaceElementKeys = classmethod(replaceElementKeys)

	def getElementKeys(self, str_in):
		"""
		Method. It returns the set of [element,key] pairs for input string.
		"""
		res = []
		patt=re.compile(r'[\w]*\[[\w]*\]')
		s_name_key = re.findall(patt,str_in)
		if len(s_name_key) == 0:
			return res
		patt_elem = re.compile(r'[\w]*(?=\[)')
		patt_key = re.compile(r'(?<=\[)[\w]*(?=\])')
		for s in s_name_key:
			elem = re.findall(patt_elem,s)[0]
			key = re.findall(patt_key,s)[0]
			res.append([elem,key])
		return res

	getElementKeys = classmethod(getElementKeys)

	def calculateString(self, str_in, localDict):
		"""
		Method. It returns a tuple (True,value) if
		the expression can be evaluated and (False,None) otherwise.
		"""
		try:
			val = eval(str_in,globals(),localDict)
			return (True, val)
		except:
			return (False, None)

	calculateString = classmethod(calculateString)

	def replaceMath(self, str_in):
		"""
		Method. It replaces math symbols
		to make them readable for python eval().
		"""
		#replace .e by .0e
		str_out = re.sub("\.e",".0e",str_in)
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
		return str_out

	replaceMath = classmethod(replaceMath)

#====================================================================

class MAD_Parser:
	""" MAD parser """

	def __init__(self):
		""" Create instance of the MAD_Parser class """
		self.__madLines =  []
		self.__accValues = []
		self.__accElements = []
		self.__accLines = []
		self.__madFilePath = ""
		#old style lattice elements and lines
		self.__lattElems = []
		self.__lattLines = []
		#the lines to ignore will start with these words
		self.__ingnoreWords = ["title","beam"]

	def __del__(self):
		del self.__madLines

	def parse(self,MADfileName):
		self.__init__()
		#1-st stage read MAD file into the lines array
		self.__madFilePath = os.path.dirname(MADfileName)
		fileName = os.path.basename(MADfileName)
		#the initialize can be recursive if there are nested MAD files
		self.initialize(fileName)
		#-----------------------------------------------------------
		#The madLines has all lines
		#Let's create values, elements, and accelerator lines arrays
		#-----------------------------------------------------------
		for line in self.__madLines:
			if(line.getType() == _accLine.getType()):
				accLine = _accLine()
				accLine.parseLine(line.getLine())
				self.__accLines.append(accLine)
			if(line.getType() == _variable.getType()):
				var = _variable()
				var.parseLine(line.getLine())
				self.__accValues.append(var)
			if(line.getType() == _element.getType()):
				elem = _element()
				elem.parseLine(line.getLine())
				self.__accElements.append(elem)
		#print "debug size values=",len(self.__accValues)
		#print "debug size elements=",len(self.__accElements)
		#print "debug size accLines=",len(self.__accLines)
		#-----------------------------------------------
		#replace all elem[key] substrings in elements by
		#variables
		#-----------------------------------------------
		accElemDict = {}
		for accElem in self.__accElements:
			accElemDict[accElem.getName()] = accElem
		accElemDictInit = accElemDict.copy()
		doNotStop = True
		while(doNotStop):
			doNotStop = False
			accElemDictCp = accElemDict.copy()
			#print "debug dict size=",len(accElemDictCp)
			for name,accElem in accElemDictCp.iteritems():
				kvs = accElem.getParameters()
				for key,val in kvs.iteritems():
					if val != None:
						resArr = StringFunctions.getElementKeys(val)
						if(len(resArr) == 0 and accElemDict.has_key(name)):
							del accElemDict[name]
						for [el,k] in resArr:
							doNotStop = True
							accElemInside = accElemDictInit[el]
							replVal = accElemInside.getParameters()[k]
							val = StringFunctions.replaceElementKeys(val,el,k,replVal)
					kvs[key] = val
			if(len(accElemDictCp) == len(accElemDict)):
				print "=========== Unresolved AccElements============"
				for name,accElem in accElemDictCp.iteritems():
					print "name=",name,"  params=",accElem.getParameters()
				print "=========== MAD File Problem ==============="
				print "=================STOP======================="
				sys.exit(1)
		#---------------------------------------------------------
		#Elements are ready!
		#Now let's substitute elements parameters in variables' expression.
		#---------------------------------------------------------
		for var in self.__accValues:
			val = var.getExpression()
			resArr = StringFunctions.getElementKeys(val)
			for [el,k] in resArr:
				accElemInside = accElemDictInit[el]
				replVal = accElemInside.getParameters()[k]
				val = StringFunctions.replaceElementKeys(val,el,k,replVal)
			var.setExpression(val)
		#-----------------------------------------------
		#now let's calculate all variables
		#-----------------------------------------------
		#replace all math cos,sin, etc by math.cos, math.sin, etc
		for var in self.__accValues:
			val = var.getExpression()
			val = StringFunctions.replaceMath(val)
			var.setExpression(val)
		#Then let's calculate numerical values.
		#They can be defined recursivelly, so we need iterations
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
			if(len(accVarDictInner) == len(accVarDict) and len(accVarDict) > 0):
				print "=========== Unresolved Variables============"
				for name,var in accVarDictInner.iteritems():
					print "name=",name,"  str=",var.getExpression()
				print "=========== MAD File Problem ==============="
				print "=================STOP======================="
				sys.exit(1)
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
						print "=============MAD File problem ==============",
						print "Problem with acc. element:",accElem.getName()
						print "Parameter name:",key
						print "Can not calculate string:",val
						print "============ STOP =========================="
						sys.exit(1)
				kvNums[key] = val_out
		#---------------------------------------------
		#Let's create all lattice elements (old style)
		#---------------------------------------------
		for accElem in self.__accElements:
			lattElem = MAD_LattElement(accElem.getName(),accElem.getElementType())
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
			lattLine = MAD_LattLine(accLine.getName())
			self.__lattLines.append(lattLine)
			lattLineDict[accLine.getName()] = lattLine
		for lattLine in self.__lattLines:
			name = lattLine.getName()
			accLine = accLineDict[name]
			childs = accLine.getComponents()
			for (child,sign) in childs:
				if(lattElemDict.has_key(child) and accLineDict.has_key(child)):
					print "=========== MAD File Problem ==============="
					print "Accelerator line and element have the same name:",child
					print "=================STOP======================="
					print "Lattice Line name=",name
					sys.exit(1)
				if((not lattElemDict.has_key(child)) and (not accLineDict.has_key(child))):
					print "=========== MAD File Problem ==============="
					print "Can not find Accelerator line and element with name:",child
					print "Lattice Line name=",name
					print "=================STOP======================="
					sys.exit(1)
				if(lattElemDict.has_key(child)):
					lattLine.addItem(lattElemDict[child],sign)
				else:
					lattLine.addItem(lattLineDict[child],sign)

	def initialize(self,mad_file_name):
		fl = open(os.path.join(self.__madFilePath, mad_file_name))
		str_local = ""
		for str in fl.readlines():
			#check if the line is a comment
			if str.find("!") == 0:
				str_local = ""
				continue
			#STOP parsing a MAD file if there is a start of mad commands
			if(str.strip().lower().find("use") == 0):
				break
			str_local = "".join([str_local,str.strip()])
			#this part deal with a comment "!" at the end of line
			str0 = str_local
			if str0.rfind("!") > 0:
				str_local = ""
				for i in xrange(str0.rfind("!")):
					str_local = "".join([str_local,str0[i]])
			str_local.strip()
			#this part deal with a comment ";" at the end of line
			str0 = str_local
			if str0.rfind(";") > 0:
				str_local = ""
				for i in xrange(str0.rfind(";")):
					str_local = "".join([str_local,str0[i]])
			str_local.strip()
			#this part deal with a continue sign at the end of line
			if str_local.endswith("&"):
				str0 = str_local
				if str0.rfind("&") > 0:
					str_local = ""
					for i in xrange(str0.rfind("&")):
						str_local = "".join([str_local,str0[i]])
				str_local.strip()
				continue
			else:
				#check the empty line
				if str_local == "":
					str_local = ""
					continue
				#now we have the line to parse
				self._init_string(str_local)
				str_local = ""
		fl.close()

	def _init_string(self,str0):
		""" The method initializes the one string """
		#Now here 5 types of string
		# 0 - unknown type of the line
		# 1 - variables calculations
		# 2 - element definition
		# 3 - MAD line definition
		# 4 - call another nested MAD file
		#Delete spaces
		str0=re.sub(r'[ ]',"",str0)
		tp = self._findLineType(str0)
		if tp == 0:
			#print "StrType =0 :",str0
			return
		if tp == 1:
			#print "StrType =1 :",str0
			madLine = _madLine()
			madLine.setLine(str0)
			madLine.setType("variable")
			self.__madLines.append(madLine)
		if tp == 2:
			#print "StrType =2 :",str0
			madLine = _madLine()
			madLine.setLine(str0)
			madLine.setType("element")
			self.__madLines.append(madLine)
		if tp == 3:
			#print "StrType =3 :",str0
			madLine = _madLine()
			madLine.setLine(str0)
			madLine.setType("accLine")
			self.__madLines.append(madLine)
		if tp == 4:
			#print "StrType =4 :",str0
			fl = self._parseNestedFileLine(str0)
			self.initialize(fl)

	def _findLineType(self,line):
		""" Return type of the string """
		#Now here 5 types of string
		# 0 - unknown type of the line
		# 1 - variables calculations
		# 2 - element definition
		# 3 - MAD line definition
		# 4 - call another nested MAD file
		stype = 0
		for word in self.__ingnoreWords:
			if(line.lower().find(word) == 0):
				return stype
		t_match = re.search(r'.*:=.*',line)
		if t_match:
			stype = 1
		t_match = re.search(r'.*=.*',line)
		if stype != 1:
			t_match = re.search(r'[\w]* *:.*',line)
			if t_match:
				stype = 2
		t_match = re.search(r'.*:.*line.*=',line.lower())
		if t_match:
			stype = 3
		t_match = re.search(r' *call *file *=',line.lower())
		if t_match:
			stype = 4
		#deal with constant defenition like AAA = 1.0
		if(stype == 0):
			t_match = re.search(r'.*=.*',line)
			if t_match:
				stype = 1
		return stype

	def _parseNestedFileLine(self, line):
		""" Returns the name of the nested file"""
		#Define delimiter
		dl="'"
		if line.find(dl) < 0:
			dl = "\""
		ind0  = line.find(dl)
		ind0 += 1
		ind1  = line.rfind(dl)
		str_res=""
		if ind0 >= ind1 or ind0 < 0 or ind1 < 0:
			print "Wrong line in the MAD file"
			print "line Call file= defines wrong name of file format"
			print "Should be : Call file = 'name of file'"
			print "Line:",line
			print "Stop."
			sys.exit (0)
		for i in range(ind0,ind1):
			str_res = "".join([str_res,line[i]])
		return str_res

	def getMAD_Lines(self):
		"""
		Method. It returns the list of the lattice lines
		that are defined in the MAD file.
		"""
		return self.__lattLines

	def getMAD_Elements(self):
		"""
		Method. It returns the list of the lattice elements
		that are defined in the MAD file.
		"""
		return self.__lattElems

	def getMAD_Variables(self):
		"""
		Method. It returns the list of the variables
		that are defined in the MAD file.
		"""
		return self.__accValues

	def getMAD_LinesDict(self):
		"""
		Method. It returns the dictionary of the lattice lines
		that are defined in the MAD file.
		"""
		dict = {}
		for lattLine in self.__lattLines:
			dict[lattLine.getName()] = lattLine
		return dict

	def getMAD_ElementsDict(self):
		"""
		Method. It returns the dictionary of the lattice elements
		that are defined in the MAD file.
		"""
		dict = {}
		for lattElem in self.__lattElems:
			dict[lattElem.getName()] = lattElem
		return dict

	def getMAD_VariablesDict(self):
		"""
		Method. It returns the dictionary of the variables
		that are defined in the MAD file.
		"""
		dict = {}
		for var in self.__accValues:
			dict[var.getName()] = var.getValue()
		return dict
