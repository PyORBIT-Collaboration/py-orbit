## \namespace orbit::sc1d
## \brief The classes and functions for foils
##
## Classes:
##  - sc1DNode -longitudinal space charge node for the TEAPOT lattices
## 
## addLongitudinalSpaceChargeNode- function to add one longitudinal space charge node to the lattice

from orbit.space_charge.sc1d.sc1DNode import SC1D_AccNode
from orbit.space_charge.sc1d.scLatticeModifications import addLongitudinalSpaceChargeNode

__all__ = []
__all__.append("sc1DNode")
__all__.append("addLongitudinalSpaceChargeNode")

