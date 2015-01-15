## \namespace orbit::impedances
## \brief The classes and functions for impedances
##
## Classes:
##   ImpedanceNode - impedance node for TEAPOT lattices
##
## Functions:
##   addImpedanceNode - function to add one impedance
##                      node to the lattice

from orbit.impedances.ImpedanceNode import LImpedance_Node
from orbit.impedances.ImpedanceNode import FreqDep_LImpedance_Node
from orbit.impedances.ImpedanceNode import BetFreqDep_LImpedance_Node

from orbit.impedances.ImpedanceNode import TImpedance_Node
from orbit.impedances.ImpedanceNode import FreqDep_TImpedance_Node
from orbit.impedances.ImpedanceNode import BetFreqDep_TImpedance_Node

from orbit.impedances.ImpedanceLatticeModifications import addImpedanceNode
from orbit.impedances.ImpedanceLatticeModifications import addImpedanceNodeAsChild

__all__ = []
__all__.append("ImpedanceNode")
__all__.append("addImpedanceNode")
__all__.append("addImpedanceNodeAsChild")

