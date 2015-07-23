"""
Module. Includes classes for error accelerator nodes.
"""

import sys
import os
import math

# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import physical constants
from orbit.utils import consts

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode,\
     AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

#import error packages
from error_base import drifti
from error_base import CoordDisplacement
from error_base import QuadKicker
from error_base import QuadKickerOsc
from error_base import DipoleKickerOsc
from error_base import LongDisplacement
from error_base import StraightRotationXY
from error_base import StraightRotationXSI
from error_base import StraightRotationXSF
from error_base import StraightRotationYSI
from error_base import StraightRotationYSF
from error_base import BendFieldI
from error_base import BendFieldF
from error_base import BendDisplacementXI
from error_base import BendDisplacementXF
from error_base import BendDisplacementYI
from error_base import BendDisplacementYF
from error_base import BendDisplacementLI
from error_base import BendDisplacementLF
from error_base import RotationI
from error_base import RotationF
from error_base import derf
from error_base import root_normal
from error_base import getGauss


# Displace the coordinates of a bunch
class CoordDisplacement(DriftTEAPOT):

    def __init__(self, dx, dxp, dy, dyp, dz, dE,\
	    name = "Coordinate Displacement"):
        """
            Constructor. Creates CoordDisplacement element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("coordinate displacement node")
        self.setLength(0.0)
        self.dx  = dx
        self.dxp = dxp
        self.dy  = dy
        self.dyp = dyp
        self.dz  = dz
        self.dE  = dE

    def trackBunch(self, bunch):
        """
            The CoordDisplacement-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        dx  = self.dx
        dxp = self.dxp
        dy  = self.dy
        dyp = self.dyp
        dz  = self.dz
        dE  = self.dE
        CoordDisplacement(bunch, dx, dxp, dy, dyp, dz, dE)

    def track(self, paramsDict):
        """ 
            The CoordDisplacement-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        dx  = self.dx
        dxp = self.dxp
        dy  = self.dy
        dyp = self.dyp
        dz  = self.dz
        dE  = self.dE
        CoordDisplacement(bunch, dx, dxp, dy, dyp, dz, dE)


# Quadrupole kick a bunch
class QuadKicker(DriftTEAPOT):

    def __init__(self, k,\
	    name = "Quad Kicker"):
        """
            Constructor. Creates QuadKicker element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("quadrupole kicker node")
        self.setLength(0.0)
        self.k = k

    def trackBunch(self, bunch):
        """
            The QuadKicker-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        k = self.k
        QuadKicker(bunch, k)

    def track(self, paramsDict):
        """ 
            The QuadKicker-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        k = self.k
        QuadKicker(bunch, k)


# Oscillating quadrupole kick a bunch
class QuadKickerOsc(DriftTEAPOT):

    def __init__(self, k, phaselength, phase\
	    name = "Quad Kicker Osc"):
        """
            Constructor. Creates QuadKickerOsc element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("oscillating quadrupole kicker node")
        self.setLength(0.0)
        self.k = k
        self.phaselength = phaselength
        self.phase = phase

    def trackBunch(self, bunch):
        """
            The QuadKickerOsc-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        k = self.k
        phaselength = self.phaselength
        phase = self.phase
        QuadKickerOsc(bunch, k, phaselength, phase)

    def track(self, paramsDict):
        """ 
            The QuadKickerOsc-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        k = self.k
        phaselength = self.phaselength
        phase = self.phase
        QuadKickerOsc(bunch, k, phaselength, phase)


# Oscillating dipole kick a bunch
class DipoleKickerOsc(DriftTEAPOT):

    def __init__(self, k, phaselength, phase\
	    name = "Dipole Kicker Osc"):
        """
            Constructor. Creates DipoleKickerOsc element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("oscillating dipole kicker node")
        self.setLength(0.0)
        self.k = k
        self.phaselength = phaselength
        self.phase = phase

    def trackBunch(self, bunch):
        """
            The DipoleKickerOsc-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        k = self.k
        phaselength = self.phaselength
        phase = self.phase
        DipoleKickerOsc(bunch, k, phaselength, phase)

    def track(self, paramsDict):
        """ 
            The DipoleKickerOsc-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        k = self.k
        phaselength = self.phaselength
        phase = self.phase
        DipoleKickerOsc(bunch, k, phaselength, phase)
