## \namespace orbit::utils
## \brief Utility classes.
##
## Classes:
## - multiDimArray - function for multi-dimension array generation
## - orbitFinalize - function to finalize ORBIT script execution
## - NamedObject - class represents an object with name
## - TypedObject - class represents an object with type
## - ParamsDictObject - class represents an object that has a parameters dictionary

from orbit.utils.multiDimArray import multiDimArray
from orbit.utils.orbitFinalize import orbitFinalize
from orbit.utils.NamedObject import NamedObject
from orbit.utils.TypedObject import TypedObject
from orbit.utils.ParamsDictObject import ParamsDictObject

__all__ = []
__all__.append("multiDimArray")
__all__.append("orbitFinalize")
__all__.append("NamedObject")
__all__.append("TypedObject")
__all__.append("ParamsDictObject")
