## \namespace orbit::errors
## \brief The classes and functions for errors
##
## Classes:
##   ErrorNode - error node for TEAPOT lattices
##
## Functions:
##   addErrorNode - function to add one error
##                  node to the lattice

from orbit.errors.ErrorNode import coorddisplacement
from orbit.errors.ErrorNode import longdisplacement
from orbit.errors.ErrorNode import straightrotationxy
from orbit.errors.ErrorNode import straightrotationxsi
from orbit.errors.ErrorNode import straightrotationxsf
from orbit.errors.ErrorNode import straightrotationysi
from orbit.errors.ErrorNode import straightrotationysf
from orbit.errors.ErrorNode import bendfieldi
from orbit.errors.ErrorNode import bendfieldf
from orbit.errors.ErrorNode import benddisplacementxi
from orbit.errors.ErrorNode import benddisplacementxf
from orbit.errors.ErrorNode import benddisplacementyi
from orbit.errors.ErrorNode import benddisplacementyf
from orbit.errors.ErrorNode import benddisplacementli
from orbit.errors.ErrorNode import benddisplacementlf
from orbit.errors.ErrorNode import rotationi
from orbit.errors.ErrorNode import rotationf
from orbit.errors.ErrorNode import dipolekicker
from orbit.errors.ErrorNode import dipolekickerosc
from orbit.errors.ErrorNode import quadkicker
from orbit.errors.ErrorNode import quadkickerosc
from orbit.errors.ErrorNode import AddErrorNode
from orbit.errors.ErrorNode import AddErrorSet

from orbit.errors.ErrorLatticeModifications import addErrorNode
from orbit.errors.ErrorLatticeModifications import addErrorNodeAsChild
from orbit.errors.ErrorLatticeModifications import addErrorNodeAsChild_I
from orbit.errors.ErrorLatticeModifications import addErrorNodeAsChild_F

__all__ = []
__all__.append("")
__all__.append("addErrorNode")
__all__.append("addErrorNodeAsChild")
__all__.append("addErrorNodeAsChild_I")
__all__.append("addErrorNodeAsChild_F")
