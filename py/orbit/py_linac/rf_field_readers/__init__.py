## \namespace orbit::sns_linac::rf_field_readers
## \brief The classes to read and analyze the SuperFish output files 
##        with the electric and magnetic fields of the RF cavities.
##        The field are assumed to have axial symmetry.
##
## Classes:
## - RF_AxisFieldAnalysis         - Class. The analyzer for the RF electric filed 
##                       on the axis of RF cavity. It will calculate TTF T,T',S,S'
## - SuperFish_3D_RF_FiledReader  - Class. Will read the SuperFish file and create 
##           3D filed. It generates the field on the axis for RF_AxisFieldAnalysis
#-------------------------------------------------------------------------

from RF_AxisFieldAnalysis import RF_AxisFieldAnalysis
from SuperFish_3D_RF_FieldReader import SuperFish_3D_RF_FieldReader

__all__ = []


# Classes
__all__.append("RF_AxisFieldAnalysis")
__all__.append("SuperFish_3D_RF_FiledReader")

