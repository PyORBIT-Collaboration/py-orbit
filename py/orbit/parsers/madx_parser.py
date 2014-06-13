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

class MADX_LattElement:
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

class MADX_Parser:
	""" MAD parser """
	
	def __init__(self):
		""" Create instance of the MAD_Parser class """
		self._accElemDict = {}
		self._sequencename = ""
		self._sequencelength = ""
		self._sequencelist = []
		#the lines to ignore will start with these words
		self.__ingnoreWords = ["title","beam", "none"]

	def parse(self,MADXfileName):
		
		self.__init__()
		#1-st stage read MAD file into the lines array
		self.__madFilePath = os.path.dirname(MADXfileName)
		fileName = os.path.basename(MADXfileName)
		#the initialize can be recursive if there are nested MAD files
		fl = open(os.path.join(self.__madFilePath, fileName))
		
		str_local = ""
		print str_local
		for str in fl.readlines():
			#print str
			#take off the ";" at end of each line
			str0 = str
			if str0.rfind(";") > 0:
				str_local = ""
				for i in xrange(str0.rfind(";")):
					str_local = "".join([str_local,str0[i]])
			str_local.strip()
			#check if the line is a comment
			if str.find("!") == 0:
				str_local = ""
				continue
			#check the empty line
			if str == "":
				str_local = ""
				continue
			#Parse element or sequence definition. 
			if re.search(r'[\w]* *:.*',str_local):
				if(str_local.rfind("sequence") >=0):
					tokens = str_local.split(":")
					subtokens = tokens[1].split()
					self._sequencename = tokens[0]
					self._sequencelength = subtokens[3]
				else:
					elem = self.parseElem(str_local)
					self._accElemDict[elem.getName()] = elem
			#Add elements to the sequence list.
			if(str_local.rfind("at") >= 0):
				if(self._sequencename==""):
					print "Warning, adding elements to sequence with no name."
				tokens = str_local.split(",")
				subtokens = tokens[1].split()
				elem_name = tokens[0]
				position = subtokens[2]
				latt_elem = self._accElemDict[elem_name]
				latt_elem.addParameter("position", position)
				latt_drift = self.makeDrift(latt_elem)
				self._sequencelist.append(latt_drift)
				self._sequencelist.append(latt_elem)

		#If the last element is not at the end of the lattice, make a drift
	
		nelems = len(self._sequencelist)
		print 'number of elems', nelems
		if(nelems == 0):
			print "Warning: Creating empty lattice."
			sys.exit(1)
		endlength = float(self._sequencelength) - (float(self._sequencelist[nelems-1].getParameter("position")) + 0.5*float(self._sequencelist[nelems-1].getParameter("l")))
		if(endlength < -1e-10):
			print "Warning: Lattice parsing resulted in a lattice with length longer than specified by sequence command."
		if(endlength > 0):
			latt_drift = MADX_LattElement("end drift", "drift")
			latt_drift.addParameter("l", endlength)
			self._sequencelist.append(latt_drift)
		
						
		
	def makeDrift(self, downstreamelem):
	
		# Now we have to creat a drift between elements because MADX
		# sequence does not include drifts.
		
		seqlength = len(self._sequencelist)
		upstreamelemlength = 0
		downstreamelemlength = 0
		
		if(seqlength == 0):
			startpos = 0
		else:
			upstreamelem = self._sequencelist[seqlength-1]
			upstreamelemlength = float(upstreamelem.getParameter("l"))
			downstreamelemlength = float(downstreamelem.getParameter("l"))
			startpos = float(upstreamelem.getParameter("position")) + 0.5*upstreamelemlength

		endpos = float(downstreamelem.getParameter("position")) - 0.5*downstreamelemlength

		driftlength = endpos - startpos
	
		name = "Drift" + "_" + str(seqlength)
		type = "drift"
		length = 0.0
		strength = 0.0
		
		if(driftlength < -1e-10):
			print "Warning: Drift between, ', upstreamelem.getName(), ' and ', downstreamelem.getName(), ' has negative length."
			print "Setting length to zero."
			lattElem = MADX_LattElement(name, type)
		else:
			lattElem = MADX_LattElement(name, type)
			lattElem.addParameter("l", driftlength)

		return lattElem

			
	def parseElem(self,line_init):
		"""
			Method. It does the first parsing of the initial string.
		"""
		name = ""
		type = ""
		length = 0.0
		strength = 0.0
		
		tokens = line_init.split(",")
		nvalues = len(tokens)
		
		if nvalues >= 1:
			subtokens = tokens[0].split(":")
			name = subtokens[0]
			type = subtokens[1].strip()
			lattElem = MADX_LattElement(name, type)
			for i in range(1, nvalues):
				subtokens = tokens[i].split(":=")
				lattElem.addParameter(subtokens[0],float(subtokens[1]))
		else:
			lattElem = MADX_LattElement(name, type)
			print "Warning: Returning empty lattice element type."
		
		#Make sure there is always a length parameter.
		if(lattElem.hasParameter("l")==0):
			lattElem.addParameter("l",0.0)
				
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

#STOP parsing a MAD file if there is a start of mad commands

