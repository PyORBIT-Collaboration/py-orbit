"""
Module. Includes classes for all PTC elements.
Built on the back of TEAPOT by J. Holmes.
"""

import sys
import os
import math

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

# import the function that creates multidimensional arrays
from orbit.utils import orbitFinalize

# import some constants
from orbit.utils import consts

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode,\
    AccActionsContainer, AccNodeBunchTracker

# import the teapot classes
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import BaseTEAPOT

# import the interface to PTC
from libptc_orbit import *


class PTC_Lattice(TEAPOT_Lattice):
	"""
	PTC Subclass of the AccLattice class.
	Inherits from the TEAPOT lattice.
	"""
	def __init__(self, name = "no name"):
		TEAPOT_Lattice.__init__(self, name)

        def readPTC(self, PTC_File):
                """
                Reads the PTC file input and initializes all structures.
                Input PTC_File is the flat PTC file.
                """
		self.setName(PTC_File)
                length_of_name = len(PTC_File)
                ptc_init_(PTC_File, length_of_name - 1)
                (betax, betay, alphax, alphay, etax, etapx) =\
                    ptc_get_twiss_init_()
		self.betax0  = betax
		self.betay0  = betay
		self.alphax0 = alphax
		self.alphay0 = alphay
		self.etax0   = etax
		self.etapx0  = etapx
                (nNodes, nHarm, lRing, gammaT) = ptc_get_ini_params_()
		self.nNodes = nNodes
		self.nHarm  = nHarm
		self.lRing  = lRing
		self.gammaT = gammaT
                for node_index in range(nNodes):
                        (length, betax, betay,\
                                 alphax, alphay, etax, etapx) =\
                                 ptc_get_twiss_for_node_(node_index)
                        elem = PTC_Node("PTC_Node")
			elem.setparams(node_index, length,\
					       betax, betay, alphax, alphay,\
					       etax, etapx)
                        self.addNode(elem)
                self.initialize()


class PTC_Node(BaseTEAPOT):
	"""
	PTC element.
	"""
	def __init__(self, name = "ptc_node"):
		"""
		Constructor. Creates a PTC element.
		"""
		BaseTEAPOT.__init__(self, name)
		self.setType("ptc_node")

	def setparams(self, orbit_ptc_node_index, length,\
			      betax, betay, alphax, alphay,\
			      etax, etapx):
		"""
		Sets element parameters.
		"""
		self.addParam("node_index", orbit_ptc_node_index)
		self.setLength(length)
		self.addParam("betax" , betax)
		self.addParam("betay" , betay)
		self.addParam("alphax", alphax)
		self.addParam("alphay", alphay)
		self.addParam("etax"  , etax)
		self.addParam("etapx" , etapx)

	def track(self, paramsDict):
		"""
		The PTC class implementation of the
		AccNodeBunchTracker class track(probe) method.
		"""
		bunch = paramsDict["bunch"]
		PhaseLength = paramsDict["length"]
		orbit_ptc_node_index = self.getParam("node_index")
		action_type = -1
		#ptc_get_task_type_(orbit_ptc_node_index, action_type)
		
		if(action_type == 1):
			print "==============================="
			print "PTC_Node.track."
			print "Energy chnage actions have not been taken."
			print "STOP."
			sys.exit(1)
		ptc_trackBunch(bunch, PhaseLength, orbit_ptc_node_index)


def setBunchParamsPTC(bunch):
	"""
	Sets the synchronous particle parameters of the bunch.
	"""
	(mass, charge, kin_energy) = ptc_get_syncpart_()
	mass = mass * consts.mass_proton
	syncPart = bunch.getSyncParticle()
	syncPart.kinEnergy(kin_energy)
	bunch.charge(charge)
	bunch.mass(mass)


def readAccelTablePTC(Acc_File):
	"""
	Gets the information for acceleration.
	"""
	length_of_name = len(Acc_File)
	ptc_read_accel_table_(Acc_File, length_of_name - 1)


def readScriptPTC(Script_File):
	"""
	Reads a PTC Script file.
	"""
	length_of_name = len(Script_File)
	ptc_script_(Script_File, length_of_name - 1)


def updateParamsPTC(lattice, bunch):
	"""
	Updates Twiss parameters of lattice.
	Updates element parameters.
	Updates synchronous particle parameters of the bunch.
	"""
	(betax, betay, alphax, alphay, etax, etapx) =\
	    ptc_get_twiss_init_()
	lattice.betax0  = betax
	lattice.betay0  = betay
	lattice.alphax0 = alphax
	lattice.alphay0 = alphay
	lattice.etax0   = etax
	lattice.etapx0  = etapx
	(nNodes, nHarm, lRing, gammaT) = ptc_get_ini_params_()
	lattice.nNodes = nNodes
	lattice.nHarm  = nHarm
	lattice.lRing  = lRing
	lattice.gammaT = gammaT
	for node in lattice.getNodes():
		node_index = node.getParam("node_index")
		length     = node.getLength()
		ptc_get_twiss_for_node_(node_index)
		node.setparams(node_index, length,\
				       betax, betay, alphax, alphay,\
				       etax, etapx)
	setBunchParamsPTC(bunch)


def synchronousSetPTC(ival):
	"""
	Calls ptc_synchronous_set_.
	"""
	if(ival >= 0):
		print "==============================="
		print "synchronousSetPTC requires ival < 0"
		print "STOP."
		sys.exit(1)
	ptc_synchronous_set_(ival)


def synchronousAfterPTC(ival):
	"""
	Calls ptc_synchronous_set_.
	"""
	if(ival >= 0):
		print "==============================="
		print "synchronousAfterPTC requires ival < 0"
		print "STOP."
		sys.exit(1)
	ptc_synchronous_after_(ival)


def trackBunchThroughLatticePTC(lattice, bunch, PhaseLength):
	"""
	Tracks a bunch through the whole lattice.
	"""
	paramsDict = {}
	paramsDict["bunch"]= bunch
	paramsDict["length"]=PhaseLength
	for node in lattice.getNodes():
		node.track(paramsDict)


def trackBunchInRangePTC(lattice, bunch, PhaseLength, indexi, indexf):
	"""
	Tracks a bunch from indexi through indexf, inclusive.
	"""
	paramsDict = {}
	paramsDict["bunch"]= bunch
	paramsDict["length"]=PhaseLength
	for node in lattice.getNodes():
		node_index = node.getParam("node_index")
		if((node_index >= indexi) and (node_index <= indexf)):
			node.track(paramsDict)
