import os
import sys
import re
import math

from orbit.utils   import orbitFinalize

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
		self.__names_type.append("aperture")
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
		self.__names_type.append("tkicker") 
		self.__names_type.append("hkick")
		self.__names_type.append("vkick")
		self.__names_type.append("rfcavity")
		self.__names_type.append("rcollimator")
		self.__names_type.append("collimator")
		self.__names_type.append("marker")
		self.__names_type.append("monitor")
		self.__names_type.append("hmonitor")
		self.__names_type.append("vmonitor")
		self.__names_type.append("dipedge")
		self.__names_type.append("elseparator")
	
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
		if(self.__names_type.count(name) == 0):
			msg = "Error creating lattice element:"
			msg = msg + os.linesep
			msg = msg + "There is no element with type:" + name
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)
		return name_in
#===============================================================

class MADX_LattElement:
	"""Class. Represents an arbitrary element in the lattice"""
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
	
	def hasParameter(self, nameOfPar):
	
		if self.__par.has_key(nameOfPar) == 0:
			return 0

		else:
			return 1
	
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
class _variable:
	"Class. Holds MADX variables."
	def __init__(self):
		"Constructor. Creates empty MAD variable."
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
		if "." in s_name:
			s_name = "".join(s_name.split("."))

		self.setName(s_name)
		self.setExpression(s_val)

#====================================================================
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
		str_out = re.sub("cos\(","math.cos(",str_out)
		str_out = re.sub("tan\(","math.tan(",str_out)
		str_out = re.sub("exp\(","math.exp(",str_out)
		str_out = re.sub("log\(","math.log(",str_out)
		str_out = re.sub("acos\(","math.acos(",str_out)
		str_out = re.sub("asin\(","math.asin(",str_out)
		str_out = re.sub("atan\(","math.atan(",str_out)
		str_out = re.sub("sqrt\(","math.sqrt(",str_out)
		str_out = re.sub("pi","math.pi",str_out)
		return str_out

	replaceMath = classmethod(replaceMath)
#====================================================================

class MADX_Parser:
	""" MADX parser """
	
	def __init__(self):
		""" Create instance of the MAD_Parser class """
		self._madxLines = [] # lines of madx file-(s)
		self._accElemDict = {}
		self._varDict = {} # dict of variables
		self._sequencename = ""
		self._sequencelength = ""
		self._sequencelist = []
		#the lines to ignore will start with these words
		self.__ingnoreWords = ["title","beam", "none", "initial"]

	def collect_madx_lines(self,fileName,madFilePath):

		#1-st stage read MAD file into the lines array
		#the initialize can be recursive if there are nested MADX files
		fl = open(os.path.join(madFilePath, fileName))

		for str in fl:
			#check if the line is a comment
			if str.find("!") == 0:
				str_local = ""
				continue
			#check if the line is a comment
			if str.find("//") == 0:
				str_local = ""
				continue
			#check the empty line
			if not str:
				str_local = ""
				continue
			# remove spaces, capital letters, words,...
			tmp=str[:].lower().split()
			str0="".join([x for x in tmp if x not in ["real","const"]]) # removes "const" and "real" from var definitions
			str0 = "".join(str0.split()) # removes all spaces

			# check if the line with ignore word
			ignore_word = [word for word in self.__ingnoreWords if word in str0 and str0.lstrip(word)!=str0]
			if ignore_word:
				print("The line starting with {} word found. line [{}] is ignored".format(ignore_word,str0))
				str_local = ""
				continue
			#take off the ";" at end of each line
			if str0.rfind(";") > 0:
				str_local = ""
				for i in range(str0.rfind(";")):
					str_local = "".join([str_local,str0[i]])
				str_local.strip()
			# deal with multi-line definitions in madx file
			# is it better to preliminary collect all the lines, join and split by ";"?
			else:
				pass
			#check if there is a nested file
			if "call" in str.lower() and "file" in str.lower():
				new_file = self._parseNestedFileLine(str)
				self.collect_madx_lines(new_file,self.madxFilePath)
			else:
				self._madxLines.append(str_local)
		

	def parse(self,MADXfileName):
		
		self.__init__()
		
		str_local = ""
		aper_warning = 0

		self.madxFilePath = os.path.dirname(MADXfileName)
		fileName = os.path.basename(MADXfileName)

		self.collect_madx_lines(fileName,self.madxFilePath)


		for str_local in self._madxLines:
			loc_flag_elem = True
		
			#Parse element/variable":="/sequence definition
			if "=" in str_local and ":" not in str_local.split("=")[0]: # here we avoid the variable parsing twice 
				loc_flag_elem = False # in order to avoid the parsing as an element
				if "at=" not in str_local:
					# we have a MADX variable here!
					var = _variable()
					var.parseLine(str_local)
					self._varDict[var.getName()] = var
				else:
					# we have the defined elem with additional params at the positioning !! not a variable, but can't be parsed as an elem !!
					tokens = str_local.split(",")
					name = tokens[0]
					elem_tmp = self._accElemDict[name]
					elem = self.parseParameters(str_local,elem_tmp)		

			if re.search(r'[\w]* *:.*',str_local):
				if(str_local.rfind("sequence") >=0):

					tokens = str_local.split(":") # is it possible to have more than two tokens here?
					self._sequencename = tokens[0]
					tmp_str = ":".join(tokens[1:])
					aux = [x.split("l=")[-1] for x in tmp_str.split(",") if "l=" in x]

					var = _variable()
					var.parseLine("SEQUENCE_LEN={}".format(aux[0]))
					self._varDict[var.getName()] = var

					localValDict =self.calculateVariables() # all variables are ready, now we can recalculate them  
					self._sequencelength = self._varDict["SEQUENCE_LEN"].getValue()

					aux = [x.split("refer=")[-1] for x in tmp_str.split(",") if "refer=" in x]
					if not aux:
						elementRefer = "centre"
					else:
						elementRefer = aux[0]

				else:
					if ":=" in str_local and "at=" not in str_local:
						if ":" not in str_local.split(":=")[0]:
							# we have a MADX variable here!
							var = _variable()
							var.parseLine(str_local)
							self._varDict[var.getName()] = var
							loc_flag_elem = False # in order to avoid the parsing as an element
					# element parsing
					if loc_flag_elem:
						if "at=" in str_local: 
							tokens = str_local.split(",") 
							tmp = ",".join([x for x in tokens if "at=" not in x]) # remove location from the string (is it necessary?)
						else:
							tmp = str_local
						elem = self.parseElem(tmp) # elem can be defined directly at the location!
						self._accElemDict[elem.getName()] = elem

			#Add elements to the sequence list (and secondary elem params to an elem).
			# if the positioning started, we can't have any variable definition (can we?)
			# we can have the definition of elements at the locations 
			if(str_local.rfind("at=") >= 0):
				if(self._sequencename==""):
					print "Warning, adding elements to sequence with no name."
				if "," in str_local.split(":")[0]:
					tokens = str_local.split(",") 
					elem_name = tokens[0]	# elem name 100% is the first position, otherwise the user did a mistake
					tmp = tokens[1:]
					aux = [x.split("at=")[-1] for x in tmp if "at=" in x]				
					position = eval(aux[0]) 
				else:	
					tokens = str_local.split(":")
					elem_name = tokens[0]	
					tmp_str = "".join(tokens[1:])
					tmp = tmp_str.split(",") 
					aux = [x.split("at=")[-1] for x in tmp if "at=" in x]
					position = eval(aux[0]) 
	
				latt_elem = self._accElemDict[elem_name]
				# he have the element, let's replace variables in parameters by numerical values here
				latt_elem = self.recalculateParameters(latt_elem,localValDict)

				# all monitors in PyORBIT have zero len by definition				
				if "monitor" in latt_elem.getType() and latt_elem.getParameter("l"):
					latt_elem.addParameter("l", 0.0)

				length = latt_elem.getParameter("l")

				if "from" in latt_elem.getParameters().keys():
					refer_elem_name = latt_elem.getParameter("from")
					refer_elem = self._accElemDict[refer_elem_name]
					position += refer_elem.getParameter("position")

				latt_elem.addParameter("position", position)
				if latt_elem.hasParameter("apertype"):
					latt_aper_entry = self.makeAperture(latt_elem)
					latt_aper_entry.addParameter("position", position-length/2.0)
					latt_aper_exit = self.makeAperture(latt_elem)
					latt_aper_exit.addParameter("position", position+length/2.0)
					latt_drift = self.makeDrift(latt_elem,elementRefer)
					self._sequencelist.append(latt_drift)
					self._sequencelist.append(latt_aper_entry)
					self._sequencelist.append(latt_elem)
					self._sequencelist.append(latt_aper_exit)
					aper_warning = aper_warning + 2
				else:
					latt_drift = self.makeDrift(latt_elem,elementRefer)
					self._sequencelist.append(latt_drift)
					self._sequencelist.append(latt_elem)

			if(str_local.rfind("endsequence") >= 0):
				#If the last element is not at the end of the lattice, make a drift		
				if not len(self._sequencelist):
					print "Warning: Creating empty lattice."
					sys.exit(1)
				else:
					#If the last element is not at the end of the lattice, make a drift
					lattEnd = MADX_LattElement("lattEnd", "marker")
					endPos = float(self._sequencelength)
					lattEnd.addParameter("position", endPos)
					lattEnd.addParameter("l", 0)
					latt_drift = self.makeDrift(lattEnd,elementRefer)
					self._sequencelist.append(latt_drift)
					self._sequencelist.append(lattEnd)
					
				endlength = float(self._sequencelength) - (float(self._sequencelist[-1].getParameter("position")) + 0.5*float(self._sequencelist[-1].getParameter("l")))
				if(endlength < -1e-10):
					print "Warning: Lattice parsing resulted in a lattice with length longer than specified by sequence command."

				
		if aper_warning >= 1:
			print "Warning, adding", aper_warning ,"aperture nodes to the teapot lattice. That will slow down the simluation."
			print "If the lost of particles on the aperture is not necessary, please use a madx file without the aperture labels."
		


	def calculateVariables(self):
		#---------------------------------------------------------
		#Now let's substitute elements parameters in variables' expression.
		#---------------------------------------------------------
		for name,var in self._varDict.iteritems():
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
		for name,var in self._varDict.iteritems():
			val = var.getExpression()
			val = StringFunctions.replaceMath(val)
			var.setExpression(val)
		#Then let's calculate numerical values.
		#They can be defined recursivelly, so we need iterations
		accVarDict = {}
		for name,var in self._varDict.iteritems():
			accVarDict[var.getName()] = var
		localValDict = {}
		doNotStop = True
		while(doNotStop):
			doNotStop = False
			accVarDictInner = accVarDict.copy()
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
				print "=========== MADX File Problem ==============="
				print "=================STOP======================="
				sys.exit(1)

		return localValDict

	def recalculateParameters(self,accElem,localValDict):
		#-----------------------------------------------
		#replace all elem[key] substrings in element by variables
		#-----------------------------------------------
		name = accElem.getName()
		accElemDictInit = self._accElemDict.copy()
		accElemDictCp = self._accElemDict.copy()
		doNotStop = True
		while(doNotStop):
			doNotStop = False
			kvs = accElem.getParameters()
			for key,val in kvs.iteritems():
				if val != None:
					tmp = "{}".format(val)
					if not tmp.split(","):
						resArr = StringFunctions.getElementKeys(tmp)
					else:
						resArr = []
						for x in tmp.split(","):
							resArr += StringFunctions.getElementKeys(x)
					if(len(resArr) == 0 and accElemDictInit.has_key(name)):
						del accElemDictInit[name]
					for [el,k] in resArr:
						doNotStop = True
						accElemInside = accElemDictCp[el]
						replVal = accElemInside.getParameters()[k]
						val = StringFunctions.replaceElementKeys(val,el,k,replVal)
				kvs[key] = val
		if(len(accElemDictCp) == len(accElemDictInit)):
			print "=========== Unresolved AccElements============"
			for name,accElem in accElemDictCp.iteritems():
				print "name=",name,"  params=",accElem.getParameters()
			print "=========== MADX File Problem ==============="
			print "=================STOP======================="
			sys.exit(1)

		#-------------------------------------------
		# Now calculate all parameters in key,string_value
		# for accelerator elements
		#--------------------------------------------

#		for name,accElem in self._accElemDict.iteritems():
		kvs = accElem.getParameters()
		for key,val in kvs.iteritems():
			val_out = None
			if val != None and key!="apertype" and key!="from":
				tmp = ("{}".format(val)).split(",")
				out = []
				for aux in tmp:
					auxRepl = StringFunctions.replaceMath(aux)
					res,val_out = StringFunctions.calculateString(auxRepl,localValDict)
					if not res:
						val_out = 0.0
						print "============= MADX File problem ==============",
						print "Problem with acc. element:",accElem.getName()
						print "Parameter name:",key,out
						print "Can not calculate string:",val
						print "Set variable", aux, "to zero"
						print "================ CONTINUE ===================\n"
					if len(tmp) == 1:
						accElem.getParameters()[key] = val_out # directly
					else:
						out+=[val_out]
				if out:
					accElem.getParameters()[key] = out

		return accElem

		
						
		
	def makeDrift(self, downstreamelem,elementRefer):
	
		# Now we have to create a drift between elements 
		if elementRefer == "entry":
			refer = [1.0,0.0]
		if elementRefer in ["centre", "center"]:
			refer = [0.5,0.5]
		if elementRefer == "exit":
			refer = [0.0,1.0]

		if elementRefer not in ("entry","centre","center","exit"):
			print "=============MADX File problem ==============",
			print "Problem with sequence, refer is:",elementRefer
			print "Refer has to be one of :", ("entry","centre","center","exit")
			print "============ STOP =========================="
			sys.exit(1)		
	
		lenUp = 0.0 # upstream
		lenDown = 0.0 # downstream
		if not self._sequencelist:
			posUp = 0.0
			lenDown = downstreamelem.getParameter("l")
		else:
			upstreamelem = self._sequencelist[-1]
			lenUp = upstreamelem.getParameter("l")
			lenDown = downstreamelem.getParameter("l")
			posUp = upstreamelem.getParameter("position")
		posDown = downstreamelem.getParameter("position")
		driftlength = abs(posDown - posUp) - refer[0]*lenUp -refer[1]*lenDown

		name = "Drift{}".format(len(self._sequencelist)//2)
		type_local = "drift"
		
		if driftlength < 0:
			print "Warning: Drift between {} and {} has negative length, value = {}".format(upstreamelem.getName(), downstreamelem.getName(),driftlength)
			print "Setting length to zero."
			lattElem = MADX_LattElement(name, type_local)
		else:
			lattElem = MADX_LattElement(name, type_local)
			lattElem.addParameter("l", driftlength)

		return lattElem
#-------------------------------------------------------------------
	def parseParameters(self,line_init, lattElem):
		"""
			Method. It finds all parameters and applies to the existing element.
		"""
		length = 0.0
		strength = 0.0 
		# search for knl, ksl, APERTURE, ...  !! :={} or = {}
		line_tmp = line_init[:]
		if "{" in line_init and "}" in line_init:

			line_init =line_init.replace("}", "!{")
			line_init="".join([x if "!" not in x else "!" for x in line_init.split("{")])
			line_init= ",".join([x for x in line_init.split(",") if "!" not in x])
			
			if "aperture" in line_tmp.lower():
				tmp = [s for s  in line_tmp.split(",") if "apertype" in s.lower()]
				if tmp:
					apertype = tmp[0]
				else:
					apertype = 'ellipse' # default					
				tmp = [i for i,s  in enumerate(line_tmp.split(",")) if "aperture" in s.lower()]
				aper = line_tmp.split(",")[tmp[0]:tmp[0]+2] 
				dim="{},{}".format(aper[0].split("{")[-1], aper[1].split("}")[0]) # super straightforward, can cause some troubles, if, for example, aperture = {a,...,b}
			

		tokens = line_init.split(",")
		for token in tokens[1:]: # name+type are already parsed, aperture and knl/ksl are excluded
			subtokens = token.split("=")
			keyParam,valParam = subtokens[0].rstrip(":"),subtokens[1]

			if "." in valParam:
				# delete dots from names of variables in references to variables
				valTmp = valParam[:]
				for x in re.split(r'[)(*/+-]', valTmp):
					if "." in x:
						tmp = "".join(x.split("."))
						if sum([s.isdigit() for s in tmp]) != len(tmp):
							aux = valParam.split(x)
							x =  "".join(x.split("."))
							valParam =x.join(aux)
				

			if "." in keyParam:
				# delete dots from names of params
				keyParam =  "".join(keyParam.split("."))

			lattElem.addParameter(keyParam,valParam)
	
		if "{" in line_tmp and "}" in line_tmp:
			if "knl" in line_tmp.lower() or "ksl" in line_tmp.lower():
				kls, poles, skews = self.makeMultiPols(line_tmp)
				lattElem.addParameter("poles", poles)
				lattElem.addParameter("skews", skews)
				lattElem.addParameter("kls", kls)
			if "aper" in line_tmp.lower():
				apertype = apertype.split("=")
				lattElem.addParameter("aperture", dim)
				lattElem.addParameter("apertype", apertype[1])

		#Make sure there is always a length parameter.
		if not lattElem.hasParameter("l"):
			lattElem.addParameter("l","0.0")
				
		return lattElem			
#-------------------------------------------------------------------			
	def parseElem(self,line_init):
		"""
			Method. It does the first parsing of the initial string.
		"""
		name = ""
		type = ""
		length = 0.0
		strength = 0.0

		tokens = line_init.split(",")
		subtokens = tokens[0].split(":")
		name = subtokens[0]
		type_local = subtokens[1]

		# elem can be defined as a child node
		paramsDict = {}
		if type_local not in self._accElemDict.keys():			
			type_upd = type_local
		else:	
			type_upd = self._accElemDict[type_local].getType()
			paramsDict = self._accElemDict[type_local].getParameters()

		if "marker" in type_upd.lower():
			type_upd = "marker"
		if "collimator" in type_upd.lower():
			type_upd = "drift"
		if "monitor" in type_upd.lower():
			type_upd = "monitor"
		if "kicker" in type_upd.lower():
			type_upd = "kicker"
		if "elseparator" in type_upd.lower():
			type_upd = "drift"
		if "instrument" in type_upd.lower():
			type_upd = "marker"
		if "rfcavity" in type_upd.lower(): # matrix doesnt exist
			type_upd = "drift"

		lattElem = MADX_LattElement(name, type_upd)
		for key,val in paramsDict.iteritems():
			lattElem.addParameter(key,val)
				
		if not name:
			lattElem = MADX_LattElement(name, type_local)
			print "Warning: Empty lattice element type."

		lattElem_upd = self.parseParameters(line_init, lattElem)
				
		return lattElem_upd

		
	def makeMultiPols(self,line_init):
		kls = []
		poles = []
		skews = []
		line_init =line_init.replace("}", "!{")
		line_init="".join([x if "!" not in x else "!".join(x.split(",")) for x in line_init.split("{")])
		tokens = [",".join(x.rstrip("!").split("!")) for x in line_init.split(",") if "!" in x and "aper" not in x]

		for token in tokens:
			subtokens = token.split("=")
			name = subtokens[0].rstrip(":")
			kls = subtokens[1].split(',')
			if len(kls)==1:
				kls+=['0.0']

			if name == "knl":
			    skews = ["0" for x in kls]
			if name == "ksl":
			    skews = ["1" for x in kls]
			poles = ["{}".format(i) for i,x in enumerate(kls)]

		#return kls, poles, skews
		return ",".join(kls), ",".join(poles), ",".join(skews)


	def makeAperture(self, downstreamelem):
	
		# Now we have to create an aperture before and after the element with the MADX label aperture
		type_local = "apertype"
		name = "aperture"
		lattElem = MADX_LattElement("Aperture", name)
		lattElem.addParameter("l", 0.0)
		dim = downstreamelem.getParameter(name)
		shape_type = downstreamelem.getParameter(type_local)
		if shape_type == "circle":
			shape = 1
		elif shape_type == "ellipse":
			shape = 2
		elif shape_type == "rectangle":
			shape = 3
		else:
			print "======== Can not create elementwith type:",shape_type
			print "You have to fix the _teapotFactory, aperture class and madx_parser."
			print "Stop."
			sys.exit(1)
		lattElem.addParameter(name, dim)
		lattElem.addParameter(type_local, shape)
		return lattElem
			
	def getSequenceName(self):
		"""
			Method. Returns name of the sequence
			"""
		return self._sequencename

	def getSequenceList(self):
		"""
			Method. Returns list of elements in the sequence 
			"""
		return self._sequencelist

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
			print "Wrong line in the MADX file"
			print "line Call file= defines wrong name of file format"
			print "Should be : Call file = 'name of file'"
			print "Line:",line
			print "Stop."
			sys.exit (0)
		for i in range(ind0,ind1):
			str_res = "".join([str_res,line[i]])
		return str_res
#STOP parsing a MADX file if there is a start of madx commands
