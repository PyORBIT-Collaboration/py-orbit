## \namespace orbit::utils
## \brief Utility classes.
##
## Classes:
## - multiDimArray    - Method. Generates multi-dimensional array.
## - orbitFinalize    - Method. Finalizes ORBIT script execution.
## - NamedObject      - Class. Represents an object with a name.
## - TypedObject      - Class. Represents an object with a type.
## - ParamsDictObject - Class. Represents an object that has a parameters dictionary.

from orbit.utils.multiDimArray    import multiDimArray
from orbit.utils.orbitFinalize    import orbitFinalize
from orbit.utils.NamedObject      import NamedObject
from orbit.utils.TypedObject      import TypedObject
from orbit.utils.ParamsDictObject import ParamsDictObject

__all__ = []
__all__.append("multiDimArray")
__all__.append("orbitFinalize")
__all__.append("NamedObject")
__all__.append("TypedObject")
__all__.append("ParamsDictObject")
