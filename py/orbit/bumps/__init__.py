## \namespace orbit::injection
## \brief These classes are for displacement bumps
##
## Classes:
## - simpleBump    - Class for simple coordinate transverse bump.
## - TDsimpleBump  - Class for time-dependent simple coordinate 
##   		     transverse bump.
## - addTeapotBumpNode - Adds a teapot bump node to a teapot lattice
## - TeapotSimpleBumpNode - Creates a teapot instance of a simple bump nodes

from bumps import simpleBump, close_orbit_bumps, TDsimpleBump
from BumpLatticeModifications import addTeapotBumpNode
from TeapotBumpNode import TeapotSimpleBumpNode, TDTeapotSimpleBumpNode

__all__ = []
__all__.append("addTeapotBumpNode")
__all__.append("TeapotSimpleBumpNode")
__all__.append("TDTeapotSimpleBumpNode")
__all__.append("simpleBump")
__all__.append("close_orbit_bumps")
__all__.append("TDsimpleBump")
