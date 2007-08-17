import os
import sys
import re
import math

class _possibleElementType:
	""" This class keeps all possible element's types """

	def __init__(self):
		self.__names_type = []
		self.__names_type.append("drift")
		self.__names_type.append("quad")
		self.__names_type.append("sbend")
		self.__names_type.append("hkicker")
		self.__names_type.append("vkicker")
		self.__names_type.append("hkick")
		self.__names_type.append("vkick")
		self.__names_type.append("marker")
		self.__names_type.append("sextupole")
		self.__names_type.append("monitor")
		self.__names_type.append("octupole")
		self.__names_type.append("quadrupole")
		self.__names_type.append("rfcavity")
		self.__names_type.append("solenoid")
		self.__names_type.append("multipole")
		self.__names_type.append("rbend")
		self.__names_type.append("rcollimator")
		self.__names_type.append("kicker")

	def __del__(self):
		del self.__names_type

	def checkType(self,name_in):
		name = name_in.lower()
		if self.__names_type.count(name) == 0:
			print "Error of creating lattice element."
			print "There can not be an element with type:",name
			print "Stop."
			sys.exit (0)
		return name_in

#===============================================================

class MAD_LattElement:
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
		""" Sets the type of the element without checking"""
		self.__type = tp

	def addParameter(self,nameOfPar,parVal):
		self.__par[nameOfPar] = parVal

	def getParameter(self,nameOfPar):
		if self.__par.has_key(nameOfPar) == 0:
			print "class MAD_LattElement, method getParameter"
			print "The name of Element =",self.__name
			print "The type of Element =",self.__type
			print "The Element's key-val =",self.__par
			print "This Element does not have Parameter=", nameOfPar
			print "Stop."
			sys.exit (0)
		return self.__par[nameOfPar]

	def getParameters(self):
		return self.__par

	def getElements(self):
		""" Returns list of elements (only one here) """
		elements = []
		elements.append(self)
		return elements

#====================================================================

class MAD_LattLine:
	""" An Arbitrary Line in the Lattice """

	def __init__(self,name):
		""" Create instance with list of lines or elements """
		self.__name = name
		self.__items = []

	def __del__(self):
		del self.__name
		del self.__items

	def getName(self):
		""" Returns name of the line """
		return self.__name

	def getType(self):
		return "LINE"

	def addItem(self,item):
		""" Adds a line or element to this line"""
		self.__items.append(item)

	def getItems(self):
		""" Returns list with elements and lines """
		return self.__items

	def getLinesDic(self):
		""" Returns the dictionary with all lattice lines inside, recursive. """
		dic = {}
		for item in self.__items:
			if(item.getType() == self.getType() or item.getType() == self.getType().lower()):
				dic[item.getName()] = item
		return dic

	def getElements(self):
		""" Returns list of elements """
		elements = []
		for item in self.__items:
			elems = item.getElements()
			for el in elems:
				elements.append(el)
		return elements

class _madLine:
	"""
	The line of a MAD file. It could be one of three types:
	variable, accelerator line, or accelerator element.
	"""

	def __init__(self):
		"""
		Constructor of the MAD file line class instance.
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
	The MAD variable class. It keeps initial string (line) from MAD file.
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

	def setElementType(self, type):
		self._elementType = type

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
		Method. It returns the set with components' names
		"""
		return self._components

	def parseLine(self,line_init):
		"""
		Method. It does the first parsing of the initial string.
		Parses the MAD file line with lattice line definition.
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
		for it in item_names:
			patt = re.compile(r'[\d]+?(?=\*)')
			n_rep = re.findall(patt,it)
			if len(n_rep) > 0:
				n_rep = int(n_rep[0])
				patt = re.compile(r'(?<=\*)[\w]+')
				s=re.findall(patt,it)
				s=s[0]
				for i in range(1,n_rep+1):
					item_names_new.append(s)
			else:
				item_names_new.append(it)
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

	def __del__(self):
		del self.__madLines

	def parse(self,MADfileName):
		self.__init__()
		#1-st stage read MAD file into the lines array
		self.__madFilePath = os.path.dirname(MADfileName)
		fileName = os.path.basename(MADfileName)
		#the initilize can be recursive if there are nested MAD files
		self.initilize(fileName)
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
		accElmDic = {}
		for accElm in self.__accElements:
			accElmDic[accElm.getName()] = accElm
		accElmDicInit = accElmDic.copy()
		doNotStop = True
		while(doNotStop):
			doNotStop = False
			accElmDicCp = accElmDic.copy()
			#print "debug dic size=",len(accElmDicCp)
			for name,accElm in accElmDicCp.iteritems():
				kvs = accElm.getParameters()
				for key,val in kvs.iteritems():
					if val != None:
						resArr = StringFunctions.getElementKeys(val)
						if(len(resArr) == 0 and accElmDic.has_key(name)):
							del accElmDic[name]
						for [el,k] in resArr:
							doNotStop = True
							accElmInside = accElmDicInit[el]
							replVal = accElmInside.getParameters()[k]
							val = StringFunctions.replaceElementKeys(val,el,k,replVal)
					kvs[key] = val
			if(len(accElmDicCp) == len(accElmDic)):
				print "=========== Unresolved AccElements============"
				for name,accElm in accElmDicCp.iteritems():
					print "name=",name,"  params=",accElm.getParameters()
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
				accElmInside = accElmDicInit[el]
				replVal = accElmInside.getParameters()[k]
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
		for accElm in self.__accElements:
			kvs = accElm.getParameters()
			kvNums = accElm.getNumParameters()
			for key,val in kvs.iteritems():
				val_out = None
				if val != None:
					res,val_out = StringFunctions.calculateString(val.lower(),localValDict)
					if(not res):
						print "=============MAD File problem ==============",
						print "Problem with acc. element:",accElm.getName()
						print "Parameter name:",key
						print "Can not calculate string:",val
						print "============ STOP =========================="
						sys.exit(1)
				kvNums[key] = val_out
		#---------------------------------------------
		#Let's create all lattice elements (old style)
		#---------------------------------------------
		for accElm in self.__accElements:
			lattElm = MAD_LattElement(accElm.getName(),accElm.getElementType())
			self.__lattElems.append(lattElm)
			kvs = accElm.getNumParameters()
			for key,val in kvs.iteritems():
				lattElm.addParameter(key,val)
		#----------------------------------------------
		#Let's create lattice lines (old style)
		#----------------------------------------------
		lattElemDict = {}
		for elm in self.__lattElems:
			lattElemDict[elm.getName()] = elm
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
			for child in childs:
				if(lattElemDict.has_key(child) and accLineDict.has_key(child)):
					print "=========== MAD File Problem ==============="
					print "Accelerator line and element have the same name:",child
					print "=================STOP======================="
				if((not lattElemDict.has_key(child)) and (not accLineDict.has_key(child))):
					print "=========== MAD File Problem ==============="
					print "Can not find Accelerator line and element with name:",child
					print "=================STOP======================="
				if(lattElemDict.has_key(child)):
					lattLine.addItem(lattElemDict[child])
				else:
					lattLine.addItem(lattLineDict[child])

	def initilize(self,mad_file_name):
		fl = open(os.path.join(self.__madFilePath, mad_file_name))
		str_local = ""
		for str in fl.readlines():
			#check if the line is a comment
			if str.find("!") == 0:
				str_local = ""
				continue
			#STOP parsing a MAD file if there is a start of mad commands
			if(re.findall("(Use)",str)):
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
		""" The method initilizes the one string """
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
			self.initilize(fl)

	def _findLineType(self,line):
		""" Return type of the string """
		#Now here 5 types of string
		# 0 - unknown type of the line
		# 1 - variables calculations
		# 2 - element definition
		# 3 - MAD line definition
		# 4 - call another nested MAD file
		type = 0
		t_match = re.search(r'.*:=.*',line)
		if t_match:
			type = 1
		t_match = re.search(r'.*=.*',line)
		if type != 1:
			t_match = re.search(r'[\w]* *:.*',line)
			if t_match:
				type = 2
		t_match = re.search(r'.*:.*line.*=',line.lower())
		if t_match:
			type = 3
		t_match = re.search(r' *call *file *=',line.lower())
		if t_match:
			type = 4
		#deal with constant defenition like AAA = 1.0
		if(type == 0):
			t_match = re.search(r'.*=.*',line)
			if t_match:
				type = 1
		return type

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

	def getMAD_LinesDic(self):
		"""
		Method. It returns the dictionary of the lattice lines
		that are defined in the MAD file.
		"""
		dic = {}
		for lattLine in self.__lattLines:
			dic[lattLine.getName()] = lattLine
		return dic

	def getMAD_ElementsDic(self):
		"""
		Method. It returns the dictionary of the lattice elements
		that are defined in the MAD file.
		"""
		dic = {}
		for lattElem in self.__lattElems:
			dic[lattElem.getName()] = lattElem
		return dic

	def getMAD_VariablesDic(self):
		"""
		Method. It returns the dictionary of the variables
		that are defined in the MAD file.
		"""
		dic = {}
		for var in self.__accValues:
			dic[var.getName()] = var.getValue()
		return dic
