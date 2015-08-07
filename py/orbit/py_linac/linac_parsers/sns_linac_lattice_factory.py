"""
The SNS Linac Lattice Factory generates the Linac Accelerator Lattice from the information
inside of the XML input file. This structure of this file is specific for the SNS.
Users from other facilities can use the same XML files with the same structure, but if they 
need something else they can create their own Factory for different structure.
Here we use XmlDataAdaptor to parse the XML file.
The SNS Linac Lattice Factory uses a predefined set of Linac Acc Elements.
"""

import os
import sys
import math

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor

from orbit.py_linac.lattice import LinacAccLattice
from orbit.py_linac.lattice import LinacAccNodes

from orbit.py_linac.lattice import BaseLinacNode, LinacNode, LinacMagnetNode, MarkerLinacNode, Drift, Quad, AbstractRF_Gap, Bend
from orbit.py_linac.lattice import DCorrectorH, DCorrectorV
from orbit.py_linac.lattice import RF_Cavity, Sequence
from orbit.py_linac.lattice import BaseRF_Gap

# import general accelerator elements
from orbit.lattice import AccNode

# import pyORBIT Python utilities classes for objects with names, types, and dictionary parameters
from orbit.utils import orbitFinalize

class SNS_LinacLatticeFactory():
	""" 
	The SNS Linac Lattice Factory generates the Linac Accelerator Lattice 
	from the XML file of the specific structure. 
	"""
	def __init__(self):
		#We need to compare positions, lengths etc. This is our delta in [m]
		self.zeroDistance = 0.00001
		#The maximal length of the drift. It will be devided if it is more than that.
		self.maxDriftLength = 1.
		
	def setMaxDriftLength(self, maxDriftLength = 1.0):
		"""
		Sets the maximal drift length that is used for 
		the purpose of the space charge calculations and diagnostics.
		"""
		self.maxDriftLength = maxDriftLength
		
	def getMaxDriftLength(self):
		"""
		Returns the maximal drift length that is used for the purpose 
		of the space charge calculations and diagnostics.
		"""
		return self.maxDriftLength
	
	def getLinacAccLattice(self,names,xml_file_name):
		"""
		Returns the linac accelerator lattice for specified sequence names and for a specified XML file.
		"""	
		if(len(names) < 1):
			msg = "The SNS_LinacLatticeFactory method getLinacAccLattice(names,xml_file_name): you have to specify the names array!"
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)
		#----- let's parse the XML file
		acc_da = XmlDataAdaptor.adaptorForFile(xml_file_name)
		return self.getLinacAccLatticeFromDA(names,acc_da)
	
	def getLinacAccLatticeFromDA(self,names,acc_da):
		"""
		Returns the linac accelerator lattice for specified sequence names.
		"""
		if(len(names) < 1):
			msg = "The SNS_LinacLatticeFactory method getLinacAccLatticeFromDA(names,): you have to specify the names array!"
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)
		#----- let's parse the XML DataAdaptor
		accSeq_da_arr = acc_da.childAdaptors()
		#let's check that the names in good order  ==start==
		seqencesLocal = accSeq_da_arr
		seqencesLocalNames = []
		for seq_da in seqencesLocal:
			seqencesLocalNames.append(seq_da.getName())
		ind_old = -1
		count = 0
		for name in names:
			ind = seqencesLocalNames.index(name)
			if(ind < 0 or (count > 0 and ind != (ind_old + 1))):
				msg = "The LinacLatticeFactory method getLinacAccLattice(names): sequence names array is wrong!"
				msg = msg + os.linesep
				msg = msg + "existing names=" + str(seqencesLocalNames)
				msg = msg + os.linesep
				msg = msg + "sequence names="+str(names)
				orbitFinalize(msg)
			ind_old = ind
			count += 1
		#	let's check that the names in good order  ==stop==			
		ind_start = seqencesLocalNames.index(names[0])
		accSeq_da_arr = accSeq_da_arr[ind_start:ind_start+len(names)]
		#----make linac lattice
		linacAccLattice = LinacAccLattice(acc_da.getName())
		#There are the folowing possible types of elements in the linac tree:
		#QUAD - quadrupole
		#RFGAP - RF Gap
		#DCH - horizontal dipole corrector
		#DCV - vertical dipole corrector
		#Marker - anything else with the length equals to 0
		#Before putting enerything into the linacAccLattice we will create sequences 
		# with all nodes.
		#----------------------------------------------------------------------
		# The DRIFTS will be generated additionally and put into right places
		#----------------------------------------------------------------------
		def positionComp(node1_da,node2_da):
			if(node1_da.getParam("pos") > node2_da.getParam("pos")):
				return 1
			else:
				if(node1_da.getParam("pos") == node2_da.getParam("pos")):
					return 0
			return -1
		accSeqs = []
		accRF_Cavs = [] 
		seqPosition = 0.
		for seq_da in accSeq_da_arr:
			#print "debug === seq=",seq_da.getName()
			accSeq = Sequence(seq_da.getName())
			accSeq.setLinacAccLattice(linacAccLattice)
			accSeq.setLength(seq_da.doubleValue("length"))
			accSeq.setPosition(seqPosition)
			seqPosition = seqPosition + accSeq.getLength()
			accSeqs.append(accSeq)
			#---- BPM frequnecy for this sequence
			bpmFrequency = seq_da.doubleValue("bpmFrequency")
			#---- create RF Cavities
			if(len(seq_da.childAdaptors("Cavities")) == 1):
				cavs_da = seq_da.childAdaptors("Cavities")[0]	
				cav_da_arr = cavs_da.childAdaptors("Cavity")
				for cav_da in cav_da_arr:
					frequency = cav_da.doubleValue("frequency")
					cav_amp = cav_da.doubleValue("ampl")
					cav_name = cav_da.stringValue("name")
					cav_pos = cav_da.doubleValue("pos")
					cav = RF_Cavity(cav_name)
					cav.setAmp(cav_amp)
					cav.setFrequency(frequency)
					cav.setPosition(cav_pos)
					accSeq.addRF_Cavity(cav)
			#----------------------------
			#node_da_arr - array of nodes. These nodes are not AccNodes. They are XmlDataAdaptor class instances
			node_da_arr = seq_da.childAdaptors("accElement")
			#put nodes in order according to the position in the sequence
			for node_da in node_da_arr:
				node_da.setParam("pos",node_da.doubleValue("pos"))
			node_da_arr.sort(positionComp)	
			#thinNodes - array of accNode nodes with zero length
			#They can be positioned inside the thick nodes, and this will be done at the end
			#of this method
			thinNodes = []
			for node_da in node_da_arr:
				params_da = node_da.childAdaptors("parameters")[0]
				node_type = node_da.stringValue("type")
				node_length = node_da.doubleValue("length")
				node_pos = node_da.getParam("pos")
				#------------QUAD-----------------
				if(node_type == "QUAD"):
					accNode = Quad(node_da.stringValue("name"))				
					accNode.setParam("dB/dr",params_da.doubleValue("field"))
					accNode.setParam("field",params_da.doubleValue("field"))
					accNode.setLength(params_da.doubleValue("effLength"))
					if(params_da.hasAttribute("poles")):
						accNode.setParam("poles",[int(x) for x in eval(params_da.stringValue("poles"))])
					if(params_da.hasAttribute("kls")):
						accNode.setParam("kls", [x for x in eval(params_da.stringValue("kls"))])
					if(params_da.hasAttribute("skews")):
						accNode.setParam("skews",[int(x) for x in eval(params_da.stringValue("skews"))])					
					if(0.5*accNode.getLength() > self.maxDriftLength):
						accNode.setnParts(2*int(0.5*accNode.getLength()/self.maxDriftLength  + 1.5 - 1.0e-12))
					accNode.setParam("pos",node_pos)
					accSeq.addNode(accNode)
				#------------BEND-----------------
				elif(node_type == "BEND"):
					accNode = Bend(node_da.stringValue("name"))                                                                                					
					if(params_da.hasAttribute("poles")):
						accNode.setParam("poles",[int(x) for x in eval(params_da.stringValue("poles"))])
					if(params_da.hasAttribute("kls")):
						accNode.setParam("kls", [x for x in eval(params_da.stringValue("kls"))])
					if(params_da.hasAttribute("skews")):
						accNode.setParam("skews",[int(x) for x in eval(params_da.stringValue("skews"))])
					accNode.setParam("ea1",params_da.doubleValue("ea1"))
					accNode.setParam("ea2",params_da.doubleValue("ea2"))
					accNode.setParam("theta",params_da.doubleValue("theta"))
					accNode.setLength(params_da.doubleValue("effLength"))
					if(0.5*accNode.getLength() > self.maxDriftLength):
						accNode.setnParts(2*int(0.5*accNode.getLength()/self.maxDriftLength  + 1.5 - 1.0e-12))
					accNode.setParam("pos",node_pos)
					accSeq.addNode(accNode)
				#------------RF_Gap-----------------	
				elif(node_type == "RFGAP"):
					accNode = BaseRF_Gap(node_da.stringValue("name"))					
					accNode.setLength(0.)
					accNode.setParam("E0TL",params_da.doubleValue("E0TL"))
					accNode.setParam("E0L",params_da.doubleValue("E0L"))
					accNode.setParam("mode",params_da.doubleValue("mode"))
					accNode.setParam("gap_phase",params_da.doubleValue("phase")*math.pi/180.)
					accNode.setParam("EzFile",params_da.stringValue("EzFile"))
					cav_name = params_da.stringValue("cavity")
					cav = accSeq.getRF_Cavity(cav_name)
					cav.addRF_GapNode(accNode)
					if(accNode.isFirstRFGap()):
						cav.setPhase(accNode.getParam("gap_phase"))
					#---- TTFs parameters
					ttfs_da = node_da.childAdaptors("TTFs")[0]
					accNode.setParam("beta_min",ttfs_da.doubleValue("beta_min"))
					accNode.setParam("beta_max",ttfs_da.doubleValue("beta_max"))					
					(polyT,polyS,polyTp,polySp) = accNode.getTTF_Polynimials()
					polyT_da = ttfs_da.childAdaptors("polyT")[0]
					polyS_da = ttfs_da.childAdaptors("polyS")[0]
					polyTp_da = ttfs_da.childAdaptors("polyTP")[0]
					polySp_da = ttfs_da.childAdaptors("polySP")[0]
					polyT.order(polyT_da.intValue("order"))
					polyS.order(polyS_da.intValue("order"))
					polyTp.order(polyTp_da.intValue("order"))
					polySp.order(polySp_da.intValue("order"))
					coef_arr = polyT_da.doubleArrayValue("pcoefs")
					for coef_ind in range(len(coef_arr)):
						polyT.coefficient(coef_ind,coef_arr[coef_ind])
					coef_arr = polyS_da.doubleArrayValue("pcoefs")
					for coef_ind in range(len(coef_arr)):
						polyS.coefficient(coef_ind,coef_arr[coef_ind])
					coef_arr = polyTp_da.doubleArrayValue("pcoefs")
					for coef_ind in range(len(coef_arr)):
						polyTp.coefficient(coef_ind,coef_arr[coef_ind])
					coef_arr = polySp_da.doubleArrayValue("pcoefs")
					for coef_ind in range(len(coef_arr)):
						polySp.coefficient(coef_ind,coef_arr[coef_ind])
					accNode.setParam("pos",node_pos)
					accSeq.addNode(accNode)
				else:
					if(node_length != 0.):
						msg = "The LinacLatticeFactory method getLinacAccLattice(names): there is a strange element!"
						msg = msg + os.linesep
						msg = msg + "name=" + node_da.stringValue("name")
						msg = msg + os.linesep
						msg = msg + "type="+node_type
						msg = msg + os.linesep
						msg = msg + "length(should be 0.)="+str(node_length)
						orbitFinalize(msg)
					#------ thin nodes analysis 
					accNode = None
					if(node_type == "DCV" or node_type == "DCH"):
						if(node_type == "DCV"): accNode = DCorrectorV(node_da.stringValue("name"))
						if(node_type == "DCH"): accNode = DCorrectorH(node_da.stringValue("name"))
						accNode.setParam("effLength",params_da.doubleValue("effLength"))
					else:
						accNode = MarkerLinacNode(node_da.stringValue("name"))
					accNode.setParam("pos",node_pos)
					thinNodes.append(accNode)
			#----- assign the thin nodes that are inside the thick nodes
			unusedThinNodes = []
			for thinNode in thinNodes:
				thinNode_pos = thinNode.getParam("pos")
				isInside = False
				for accNode in accSeq.getNodes():
					length = accNode.getLength()
					if(length > 0.):
						pos = accNode.getParam("pos")
						if(thinNode_pos >= (pos-length/2) and thinNode_pos <= (pos+length/2)):
							isInside = True
							delta_pos = thinNode_pos - (pos-length/2)
							s_path = 0.
							part_ind_in = -1
							for part_ind in range(accNode.getnParts()):
								part_ind_in = part_ind
								s_path += accNode.getLength(part_ind)
								if(delta_pos <= s_path + self.zeroDistance): 
									break
							accNode.addChildNode(thinNode, place = AccNode.BODY, part_index = part_ind_in , place_in_part = AccNode.AFTER)
							thinNode.setParam("pos",(pos-length/2)+s_path)
				if(not isInside):
					unusedThinNodes.append(thinNode)
			thinNodes = unusedThinNodes
			newAccNodes = accSeq.getNodes()[:] + thinNodes
			newAccNodes.sort(positionComp)
			accSeq.setNodes(newAccNodes)
			#insert the drifts ======================start ===========================
			#-----now check the integrety quads and rf_gaps should not overlap
			#-----and create drifts			
			copyAccNodes = accSeq.getNodes()[:]		
			firstNode = copyAccNodes[0]
			lastNode = copyAccNodes[len(copyAccNodes)-1]
			driftNodes_before = []
			driftNodes_after = []
			#insert the drift before the first element if its half length less than its position
			if(math.fabs(firstNode.getLength()/2.0 - firstNode.getParam("pos")) > self.zeroDistance):
				if(firstNode.getLength()/2.0 > firstNode.getParam("pos")):
					msg = "The LinacLatticeFactory method getLinacAccLattice(names): the first node is too long!"
					msg = msg + os.linesep
					msg = msg + "name=" + firstNode.getName()
					msg = msg + os.linesep
					msg = msg + "type=" + firstNode.getType()
					msg = msg + os.linesep
					msg = msg + "length=" + str(firstNode.getLength())
					msg = msg + os.linesep
					msg = msg + "pos=" + str(firstNode.getParam("pos"))						
					orbitFinalize(msg)
				else:
					driftNodes = []
					driftLength = firstNode.getParam("pos") - firstNode.getLength()/2.0
					nDrifts = int(driftLength/self.maxDriftLength) + 1
					driftLength = driftLength/nDrifts
					for idrift in range(nDrifts):
						drift = Drift(accSeq.getName()+":"+firstNode.getName()+":"+str(idrift+1)+":drift")
						drift.setLength(driftLength)
						drift.setParam("pos",0.+drift.getLength()*(idrift+0.5))
						driftNodes.append(drift)
					driftNodes_before = driftNodes
			#insert the drift after the last element if its half length less + position is less then the sequence length
			if(math.fabs(lastNode.getLength()/2.0 + lastNode.getParam("pos") - accSeq.getLength()) > self.zeroDistance):
				if(lastNode.getLength()/2.0 + lastNode.getParam("pos") > accSeq.getLength()):
					msg = "The LinacLatticeFactory method getLinacAccLattice(names): the last node is too long!"
					msg = msg + os.linesep
					msg = msg + "name=" + lastNode.getName()
					msg = msg + os.linesep
					msg = msg + "type=" + lastNode.getType()
					msg = msg + os.linesep
					msg = msg + "length=" + str(lastNode.getLength())
					msg = msg + os.linesep
					msg = msg + "pos=" + str(lastNode.getParam("pos"))					
					msg = msg + os.linesep
					msg = msg + "sequence name=" + accSeq.getName()				
					msg = msg + os.linesep
					msg = msg + "sequence length=" + str(accSeq.getLength())			
					orbitFinalize(msg)
				else:
					driftNodes = []
					driftLength = accSeq.getLength() - (lastNode.getParam("pos") + lastNode.getLength()/2.0)
					nDrifts = int(driftLength/self.maxDriftLength) + 1
					driftLength = driftLength/nDrifts
					for idrift in range(nDrifts):
						drift = Drift(accSeq.getName()+":"+lastNode.getName()+":"+str(idrift+1)+":drift")
						drift.setLength(driftLength)
						drift.setParam("pos",lastNode.getParam("pos")+lastNode.getLength()/2.0 + drift.getLength()*(idrift+0.5))
						driftNodes.append(drift)
					driftNodes_after = driftNodes	
			#now move on and generate drifts between (i,i+1) nodes from copyAccNodes				
			newAccNodes = driftNodes_before
			for node_ind in range(len(copyAccNodes)-1):
				accNode0 = copyAccNodes[node_ind]
				newAccNodes.append(accNode0)
				accNode1 = copyAccNodes[node_ind+1]
				dist = accNode1.getParam("pos") - accNode1.getLength()/2 - (accNode0.getParam("pos") + accNode0.getLength()/2)
				if(dist < 0.):
					msg = "The LinacLatticeFactory method getLinacAccLattice(names): two nodes are overlapping!"
					msg = msg + os.linesep
					msg = msg + "sequence name=" + accSeq.getName()		
					msg = msg + os.linesep					
					msg = msg + "node 0 name=" + accNode0.getName() + " pos="+ str(accNode0.getParam("pos")) + " L="+str(accNode0.getLength())
					msg = msg + os.linesep
					msg = msg + "node 1 name=" + accNode1.getName() + " pos="+ str(accNode1.getParam("pos")) + " L="+str(accNode1.getLength())			
					msg = msg + os.linesep
					orbitFinalize(msg)				
				elif(dist > self.zeroDistance):					
					driftNodes = []					
					nDrifts = int(dist/self.maxDriftLength) + 1
					driftLength = dist/nDrifts				
					for idrift in range(nDrifts):
						drift = Drift(accSeq.getName()+":"+accNode0.getName()+":"+str(idrift+1)+":drift")
						drift.setLength(driftLength)
						drift.setParam("pos",accNode0.getParam("pos")+accNode0.getLength()*0.5+drift.getLength()*(idrift+0.5))
						driftNodes.append(drift)
					newAccNodes += driftNodes
				else:
					pass				
			newAccNodes.append(lastNode)
			newAccNodes += driftNodes_after
			accSeq.setNodes(newAccNodes)
			#insert the drifts ======================stop ===========================
			#add all AccNodes to the linac lattice
			for accNode in accSeq.getNodes():
				linacAccLattice.addNode(accNode)
		#------- finalize the lattice construction
		linacAccLattice.initialize()
		return linacAccLattice


