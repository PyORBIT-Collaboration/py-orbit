"""
The JPARC Linac Lattice Factory generates the Linac Accelerator Lattice from the information
inside of the XML input file. This structure of this file is very similar to SNS.
The difference is the order of accelerator sequences. In JPARC file they are in  arbitrary order.
The user will combine them on the application script level. The reason for this: JPARC linac
has several different dump sequences, so the order check that exists in the SNS lattice factory 
could fail for JPARC. It is responsibility of the user to check the order and names of the 
sequences in the lattice.
"""

import os
import sys
import math

from sns_linac_lattice_factory import SNS_LinacLatticeFactory

class JPARC_LinacLatticeFactory(SNS_LinacLatticeFactory):
	""" 
	The JPARC Linac Lattice Factory generates the Linac Accelerator Lattice 
	from the XML file of the specific structure. 
	"""
	def __init__(self):
		SNS_LinacLatticeFactory.__init__(self)
		
		
	def filterSequences_and_OptionalCheck(self,accSeq_da_arr,names):
		"""
		This method is overridden here.
		This method will filter the sequences according to names list.
		It returns the filtered array with data adapters with the names
		that are in the names array.
		"""
		accSeq_da_fltr_arr = []
		for name in names:
			for accSeq_da in accSeq_da_arr:
				if(accSeq_da.getName() == name):
					accSeq_da_fltr_arr.append(accSeq_da)
		return accSeq_da_fltr_arr
				