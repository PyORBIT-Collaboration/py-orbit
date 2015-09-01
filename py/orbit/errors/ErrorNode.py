"""
Module. Includes classes for error accelerator nodes.
"""

import sys
import os
import math

# import random number generators
from random import random
from random import gauss

# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import physical constants
from orbit.utils import consts

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode,\
     AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT
from orbit.teapot import KickTEAPOT
from orbit.teapot import SolenoidTEAPOT
from orbit.teapot import MultipoleTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import BendTEAPOT

# import error packages
from error_base import CoordDisplacement
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
from error_base import DipoleKickerOsc
from error_base import QuadKicker
from error_base import QuadKickerOsc

# import node adders
from orbit.errors.ErrorLatticeModifications import addErrorNode
from orbit.errors.ErrorLatticeModifications import addErrorNodeAsChild
from orbit.errors.ErrorLatticeModifications import addErrorNodeAsChild_I
from orbit.errors.ErrorLatticeModifications import addErrorNodeAsChild_F

# Displace the coordinates of a bunch
class coorddisplacement(DriftTEAPOT):

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


# Longitudinally displace a bunch
class longdisplacement(DriftTEAPOT):

    def __init__(self, ds,\
	    name = "Longitudinal Displacement"):
        """
            Constructor. Creates LongDisplacement element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("longitudinal displacement node")
        self.setLength(0.0)
        self.ds = ds

    def trackBunch(self, bunch):
        """
            The LongDisplacement-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        ds = self.ds
        LongDisplacement(bunch, ds)

    def track(self, paramsDict):
        """ 
            The LongDisplacement-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        ds = self.ds
        LongDisplacement(bunch, ds)


# XY rotate a bunch
class straightrotationxy(DriftTEAPOT):

    def __init__(self, angle,\
	    name = "XY Rotation"):
        """
            Constructor. Creates StraightRotationXY element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("xy rotation node")
        self.setLength(0.0)
        self.angle = angle

    def trackBunch(self, bunch):
        """
            The StraightRotationXY-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        angle = self.angle
        StraightRotationXY(bunch, angle)

    def track(self, paramsDict):
        """ 
            The StraightRotationXY-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        angle = self.angle
        StraightRotationXY(bunch, angle)


# XS rotate a bunch entering element
class straightrotationxsi(DriftTEAPOT):

    def __init__(self, angle, lengthelt,\
	    name = "XSI Rotation"):
        """
            Constructor. Creates StraightRotationXSI element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("xsi rotation node")
        self.setLength(0.0)
        self.angle = angle
        self.lengthelt = lengthelt

    def trackBunch(self, bunch):
        """
            The StraightRotationXSI-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        angle = self.angle
        lengthelt = self.lengthelt
        StraightRotationXSI(bunch, angle, lengthelt)

    def track(self, paramsDict):
        """ 
            The StraightRotationXSI-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        angle = self.angle
        lengthelt = self.lengthelt
        StraightRotationXSI(bunch, angle, lengthelt)


# XS rotate a bunch leaving element
class straightrotationxsf(DriftTEAPOT):

    def __init__(self, angle, lengthelt,\
	    name = "XSF Rotation"):
        """
            Constructor. Creates StraightRotationXSF element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("xsf rotation node")
        self.setLength(0.0)
        self.angle = angle
        self.lengthelt = lengthelt

    def trackBunch(self, bunch):
        """
            The StraightRotationXSF-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        angle = self.angle
        lengthelt = self.lengthelt
        StraightRotationXSF(bunch, angle, lengthelt)

    def track(self, paramsDict):
        """ 
            The StraightRotationXSF-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        angle = self.angle
        lengthelt = self.lengthelt
        StraightRotationXSF(bunch, angle, lengthelt)


# YS rotate a bunch entering element
class straightrotationysi(DriftTEAPOT):

    def __init__(self, angle, lengthelt,\
	    name = "YSI Rotation"):
        """
            Constructor. Creates StraightRotationYSI element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("ysi rotation node")
        self.setLength(0.0)
        self.angle = angle
        self.lengthelt = lengthelt

    def trackBunch(self, bunch):
        """
            The StraightRotationYSI-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        angle = self.angle
        lengthelt = self.lengthelt
        StraightRotationYSI(bunch, angle, lengthelt)

    def track(self, paramsDict):
        """ 
            The StraightRotationYSI-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        angle = self.angle
        lengthelt = self.lengthelt
        StraightRotationYSI(bunch, angle, lengthelt)


# YS rotate a bunch leaving element
class straightrotationysf(DriftTEAPOT):

    def __init__(self, angle, lengthelt,\
	    name = "YSF Rotation"):
        """
            Constructor. Creates StraightRotationYSF element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("ysf rotation node")
        self.setLength(0.0)
        self.angle = angle
        self.lengthelt = lengthelt

    def trackBunch(self, bunch):
        """
            The StraightRotationYSF-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        angle = self.angle
        lengthelt = self.lengthelt
        StraightRotationYSF(bunch, angle, lengthelt)

    def track(self, paramsDict):
        """ 
            The StraightRotationYSF-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        angle = self.angle
        lengthelt = self.lengthelt
        StraightRotationYSF(bunch, angle, lengthelt)


# Bend field strength error to a bunch entering element
class bendfieldi(DriftTEAPOT):

    def __init__(self, drho,\
	    name = "BendI Field"):
        """
            Constructor. Creates BendFieldI element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("bendi field node")
        self.setLength(0.0)
        self.drho = drho

    def trackBunch(self, bunch):
        """
            The BendFieldI-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        drho = self.drho
        BendFieldI(bunch, drho)

    def track(self, paramsDict):
        """ 
            The BendFieldI-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        drho = self.drho
        BendFieldI(bunch, drho)


# Bend field strength error to a bunch leaving element
class bendfieldf(DriftTEAPOT):

    def __init__(self, drho,\
	    name = "BendF Field"):
        """
            Constructor. Creates BendFieldF element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("bendf field node")
        self.setLength(0.0)
        self.drho = drho

    def trackBunch(self, bunch):
        """
            The BendFieldF-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        drho = self.drho
        BendFieldF(bunch, drho)

    def track(self, paramsDict):
        """ 
            The BendFieldF-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        drho = self.drho
        BendFieldF(bunch, drho)


# X displacement error to a bunch entering bend
class benddisplacementxi(DriftTEAPOT):

    def __init__(self, theta, disp,\
	    name = "BendXI Displacement"):
        """
            Constructor. Creates BendDisplacementXI element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("xi bend displacement node")
        self.setLength(0.0)
        self.theta = theta
        self.disp = disp

    def trackBunch(self, bunch):
        """
            The BendDisplacementXI-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        theta = self.theta
        disp = self.disp
        BendDisplacementXI(bunch, theta, disp)

    def track(self, paramsDict):
        """ 
            The BendDisplacementXI-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        theta = self.theta
        disp = self.disp
        BendDisplacementXI(bunch, theta, disp)


# X displacement error to a bunch leaving bend
class benddisplacementxf(DriftTEAPOT):

    def __init__(self, theta, disp,\
	    name = "BendXF Displacement"):
        """
            Constructor. Creates BendDisplacementXF element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("xf bend displacement node")
        self.setLength(0.0)
        self.theta = theta
        self.disp = disp

    def trackBunch(self, bunch):
        """
            The BendDisplacementXF-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        theta = self.theta
        disp = self.disp
        BendDisplacementXF(bunch, theta, disp)

    def track(self, paramsDict):
        """ 
            The BendDisplacementXF-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        theta = self.theta
        disp = self.disp
        BendDisplacementXF(bunch, theta, disp)


# Y displacement error to a bunch entering bend
class benddisplacementyi(DriftTEAPOT):

    def __init__(self, disp,\
	    name = "BendYI Displacement"):
        """
            Constructor. Creates BendDisplacementYI element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("yi bend displacement node")
        self.setLength(0.0)
        self.disp = disp

    def trackBunch(self, bunch):
        """
            The BendDisplacementYI-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        disp = self.disp
        BendDisplacementYI(bunch, disp)

    def track(self, paramsDict):
        """ 
            The BendDisplacementYI-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        disp = self.disp
        BendDisplacementYI(bunch, disp)


# Y displacement error to a bunch leaving bend
class benddisplacementyf(DriftTEAPOT):

    def __init__(self, disp,\
	    name = "BendYF Displacement"):
        """
            Constructor. Creates BendDisplacementYF element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("yf bend displacement node")
        self.setLength(0.0)
        self.disp = disp

    def trackBunch(self, bunch):
        """
            The BendDisplacementYF-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        disp = self.disp
        BendDisplacementYF(bunch, disp)

    def track(self, paramsDict):
        """ 
            The BendDisplacementYF-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        disp = self.disp
        BendDisplacementYF(bunch, disp)


# L displacement error to a bunch entering bend
class benddisplacementli(DriftTEAPOT):

    def __init__(self, theta, disp,\
	    name = "BendLI Displacement"):
        """
            Constructor. Creates BendDisplacementLI element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("li bend displacement node")
        self.setLength(0.0)
        self.theta = theta
        self.disp = disp

    def trackBunch(self, bunch):
        """
            The BendDisplacementLI-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        theta = self.theta
        disp = self.disp
        BendDisplacementLI(bunch, theta, disp)

    def track(self, paramsDict):
        """ 
            The BendDisplacementLI-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        theta = self.theta
        disp = self.disp
        BendDisplacementLI(bunch, theta, disp)


# L displacement error to a bunch leaving bend
class benddisplacementlf(DriftTEAPOT):

    def __init__(self, theta, disp,\
	    name = "BendLF Displacement"):
        """
            Constructor. Creates BendDisplacementLF element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("lf bend displacement node")
        self.setLength(0.0)
        self.theta = theta
        self.disp = disp

    def trackBunch(self, bunch):
        """
            The BendDisplacementLF-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        theta = self.theta
        disp = self.disp
        BendDisplacementLF(bunch, theta, disp)

    def track(self, paramsDict):
        """ 
            The BendDisplacementLF-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        theta = self.theta
        disp = self.disp
        BendDisplacementLF(bunch, theta, disp)


# General rotation error to a bunch entering element
class rotationi(DriftTEAPOT):

    def __init__(self, angle, rhoi, theta, lengthelt, et, rotype,\
	    name = "RotationI General"):
        """
            Constructor. Creates RotationI element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("generali rotation node")
        self.setLength(0.0)
        self.angle = angle
        self.rhoi = rhoi
        self.theta = theta
        self.lengthelt = lengthelt
        self.et = et
        self.rotype = rotype

    def trackBunch(self, bunch):
        """
            The RotationI-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        angle = self.angle
        rhoi = self.rhoi
        theta = self.theta
        lengthelt = self.lengthelt
        et = self.et
        rotype = self.rotype
        RotationI(bunch, angle, rhoi, theta, lengthelt, et, rotype)

    def track(self, paramsDict):
        """ 
            The RotationI-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        angle = self.angle
        rhoi = self.rhoi
        theta = self.theta
        lengthelt = self.lengthelt
        et = self.et
        rotype = self.rotype
        RotationI(bunch, angle, rhoi, theta, lengthelt, et, rotype)


# General rotation error to a bunch leaving element
class rotationf(DriftTEAPOT):

    def __init__(self, angle, rhoi, theta, lengthelt, et, rotype,\
	    name = "RotationF General"):
        """
            Constructor. Creates RotationF element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("generalf rotation node")
        self.setLength(0.0)
        self.angle = angle
        self.rhoi = rhoi
        self.theta = theta
        self.lengthelt = lengthelt
        self.et = et
        self.rotype = rotype

    def trackBunch(self, bunch):
        """
            The RotationF-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        angle = self.angle
        rhoi = self.rhoi
        theta = self.theta
        lengthelt = self.lengthelt
        et = self.et
        rotype = self.rotype
        RotationF(bunch, angle, rhoi, theta, lengthelt, et, rotype)

    def track(self, paramsDict):
        """ 
            The RotationF-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        angle = self.angle
        rhoi = self.rhoi
        theta = self.theta
        lengthelt = self.lengthelt
        et = self.et
        rotype = self.rotype
        RotationF(bunch, angle, rhoi, theta, lengthelt, et, rotype)


# Dipole kick a bunch
class dipolekicker(DriftTEAPOT):

    def __init__(self, dxp, dyp,\
	    name = "Dipole Kicker"):
        """
            Constructor. Creates DipoleKicker element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("dipole kicker node")
        self.setLength(0.0)
        self.dxp = dxp
        self.dyp = dyp

    def setkick(self, dxp, dyp):
        """
            Sets or changes kick values.
        """
        self.dxp = dxp
        self.dyp = dyp

    def trackBunch(self, bunch):
        """
            The DipoleKicker-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        dxp = self.dxp
        dyp = self.dyp
        CoordDisplacement(bunch, 0.0, dxp, 0.0, dyp, 0.0, 0.0)

    def track(self, paramsDict):
        """ 
            The DipoleKicker-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        dxp = self.dxp
        dyp = self.dyp
        CoordDisplacement(bunch, 0.0, dxp, 0.0, dyp, 0.0, 0.0)


# Oscillating dipole kick a bunch
class dipolekickerosc(DriftTEAPOT):

    def __init__(self, k, phaselength, phase,\
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


# Quadrupole kick a bunch
class quadkicker(DriftTEAPOT):

    def __init__(self, k,\
	    name = "Quad Kicker"):
        """
            Constructor. Creates QuadKicker element.
        """
        DriftTEAPOT.__init__(self, name)
        self.setType("quadrupole kicker node")
        self.setLength(0.0)
        self.k = k

    def setkick(self, k):
        """
            Sets or changes kick values.
        """
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
class quadkickerosc(DriftTEAPOT):

    def __init__(self, k, phaselength, phase,\
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


######################################################################


# Add an error node
class AddErrorNode():

    def __init__(self, lattice, positioni, positionf, paramsDict,\
	    name = "Error Node"):
        """
            Constructor. Adds the node.
        """
        self.lattice = lattice
        self.positioni = positioni
        self.positionf = positionf
        (self.nodeIndexi, self.zi, f) = FindNode(lattice, positioni)
        (self.nodeIndexf, i, self.zf) = FindNode(lattice, positionf)
        self.localDict = paramsDict
        self.errtype = self.localDict["errtype"]
        if(self.errtype == "StraightError"):
        	self.AddStraightError()
        if(self.errtype == "FieldError"):
        	self.AddFieldError()
        if(self.errtype == "BendDisplacementError"):
        	self.AddBendDisplacementError()
        if(self.errtype == "RotationError"):
        	self.AddRotationError()

    def AddStraightError(self):
    	"""
    		Adds a StraightError to nodes
    		between nodeIndexi and nodeIndexf.
    	"""
    	lattice = self.lattice
    	lattice.initialize()
    	nodeIndexi = self.nodeIndexi
    	nodeIndexf = self.nodeIndexf
    	Indexf     = nodeIndexf + 1
    	if(nodeIndexi > nodeIndexf):
    		Indexf = len(lattice.getNodes())
    		for node in lattice.getNodes()[0 : nodeIndexf + 1]:
    			if(isinstance(node, BendTEAPOT)):
    				print "Bend node = ", node.getName(), " type = ",\
    				node.getType(), " L = ", node.getLength()
    				orbitFinalize("Can't add StraightError around a Bend! Stop!")
    	for node in lattice.getNodes()[nodeIndexi : Indexf]:
    		if(isinstance(node, BendTEAPOT)):
    			print "Bend node = ", node.getName(), " type = ",\
    			node.getType(), " L = ", node.getLength()
    			orbitFinalize("Can't add StraightError around a Bend! Stop!")
    	nodei = lattice.getNodes()[nodeIndexi]
    	nodef = lattice.getNodes()[nodeIndexf]
    	sample = self.localDict["sample"]
    	if(self.localDict["subtype"] == "TransDisp"):
    		multdx = 1.0
    		multdy = 1.0
    		if(sample == "Uniform"):
    			minimum = self.localDict["minimum"]
    			maximum = self.localDict["maximum"]
    			multdx  = minimum + (maximum - minimum) * random()
    			multdy  = minimum + (maximum - minimum) * random()
    		if(sample == "Gaussian"):
    			mean   = self.localDict["mean"]
    			sigma  = self.localDict["sigma"]
    			multdx = gauss(mean, sigma)
    			multdy = gauss(mean, sigma)
    		dx = self.localDict["dx"] * multdx
    		dy = self.localDict["dy"] * multdy
    		errori = coorddisplacement( dx, 0.0,  dy, 0.0, 0.0, 0.0)
    		errorf = coorddisplacement(-dx, 0.0, -dy, 0.0, 0.0, 0.0)
    	if(self.localDict["subtype"] == "LongDisp"):
    		multds = 1.0
    		if(sample == "Uniform"):
    			minimum = self.localDict["minimum"]
    			maximum = self.localDict["maximum"]
    			multds  = minimum + (maximum - minimum) * random()
    		if(sample == "Gaussian"):
    			mean   = self.localDict["mean"]
    			sigma  = self.localDict["sigma"]
    			multds = gauss(mean, sigma)
    		ds = self.localDict["ds"] * multds
    		errori = longdisplacement(ds)
    		errorf = longdisplacement(-ds)
    	if(self.localDict["subtype"] == "XYRot"):
    		multangle = 1.0
    		if(sample == "Uniform"):
    			minimum = self.localDict["minimum"]
    			maximum = self.localDict["maximum"]
    			multangle  = minimum + (maximum - minimum) * random()
    		if(sample == "Gaussian"):
    			mean   = self.localDict["mean"]
    			sigma  = self.localDict["sigma"]
    			multangle = gauss(mean, sigma)
    		angle = self.localDict["angle"] * multangle
    		errori = straightrotationxy( angle)
    		errorf = straightrotationxy(-angle)
    	if(self.localDict["subtype"] == "XSRot"):
    		multangle = 1.0
    		if(sample == "Uniform"):
    			minimum = self.localDict["minimum"]
    			maximum = self.localDict["maximum"]
    			multangle  = minimum + (maximum - minimum) * random()
    		if(sample == "Gaussian"):
    			mean   = self.localDict["mean"]
    			sigma  = self.localDict["sigma"]
    			multangle = gauss(mean, sigma)
    		angle = self.localDict["angle"] * multangle
    		if(self.zf >= self.zi):
    			lengtherr = self.zf - self.zi
    		else:
    			lengtherr = lattice.getLength() - self.zf + self.zi
    		errori = straightrotationxsi(angle, lengtherr)
    		errorf = straightrotationxsf(angle, lengtherr)
    	if(self.localDict["subtype"] == "YSRot"):
    		multangle = 1.0
    		if(sample == "Uniform"):
    			minimum = self.localDict["minimum"]
    			maximum = self.localDict["maximum"]
    			multangle  = minimum + (maximum - minimum) * random()
    		if(sample == "Gaussian"):
    			mean   = self.localDict["mean"]
    			sigma  = self.localDict["sigma"]
    			multangle = gauss(mean, sigma)
    		angle = self.localDict["angle"] * multangle
    		if(self.zf >= self.zi):
    			lengtherr = self.zf - self.zi
    		else:
    			lengtherr = lattice.getLength() - self.zf + self.zi
    		errori = straightrotationysi(angle, lengtherr)
    		errorf = straightrotationysf(angle, lengtherr)
    	addErrorNodeAsChild_I(lattice, nodei, errori)
    	addErrorNodeAsChild_F(lattice, nodef, errorf)

    def AddFieldError(self):
    	"""
    		Adds a FieldError to nodes at nodeIndexi to nodeIndexf
    	"""
    	lattice = self.lattice
    	lattice.initialize()
    	nodeIndexi = self.nodeIndexi
    	nodeIndexf = self.nodeIndexf
    	Indexf     = nodeIndexf + 1
    	if(nodeIndexi > nodeIndexf):
    		Indexf = len(lattice.getNodes())
    	sample = self.localDict["sample"]
    	multfrac = 1.0
    	if(sample == "Uniform"):
    		minimum = self.localDict["minimum"]
    		maximum = self.localDict["maximum"]
    		multfrac  = minimum + (maximum - minimum) * random()
    	if(sample == "Gaussian"):
    		mean   = self.localDict["mean"]
    		sigma  = self.localDict["sigma"]
    		multfrac = gauss(mean, sigma)
    	fracerr = self.localDict["fracerr"] * multfrac
    	if(self.localDict["subtype"] == "KickField"):
    		if(nodeIndexi > nodeIndexf):
    			for node in lattice.getNodes()[0 : nodeIndexf + 1]:
    				if(not isinstance(node, KickTEAPOT)):
    					print "node = ", node.getName(), " type = ",\
    					node.getType(), " L = ", node.getLength()
    					orbitFinalize("Field Error: Wanted a Kick node! Stop!")
    				kx = node.getParam("kx") * (1.0 + fracerr)
    				ky = node.getParam("ky") * (1.0 + fracerr)
    				node.setParam("kx", kx)
    				node.setParam("ky", ky)
    		for node in lattice.getNodes()[nodeIndexi : Indexf]:
    			if(not isinstance(node, KickTEAPOT)):
    				print "node = ", node.getName(), " type = ",\
    				node.getType(), " L = ", node.getLength()
    				orbitFinalize("Field Error: Wanted a Kick node! Stop!")
    			kx = node.getParam("kx") * (1.0 + fracerr)
    			ky = node.getParam("ky") * (1.0 + fracerr)
    			node.setParam("kx", kx)
    			node.setParam("ky", ky)
    	if(self.localDict["subtype"] == "SolenoidField"):
    		if(nodeIndexi > nodeIndexf):
    			for node in lattice.getNodes()[0 : nodeIndexf + 1]:
    				if(not isinstance(node, SolenoidTEAPOT)):
    					print "node = ", node.getName(), " type = ",\
    					node.getType(), " L = ", node.getLength()
    					orbitFinalize("Field Error: Wanted a Solenoid node! Stop!")
    				B = node.getParam("B") * (1.0 + fracerr)
    				node.setParam("B", B)
    		for node in lattice.getNodes()[nodeIndexi : Indexf]:
    			if(not isinstance(node, SolenoidTEAPOT)):
    				print "node = ", node.getName(), " type = ",\
    				node.getType(), " L = ", node.getLength()
    				orbitFinalize("Field Error: Wanted a Solenoid node! Stop!")
    			B = node.getParam("B") * (1.0 + fracerr)
    			node.setParam("B", B)
    	if(self.localDict["subtype"] == "MultipoleField"):
    		if(nodeIndexi > nodeIndexf):
    			for node in lattice.getNodes()[0 : nodeIndexf + 1]:
    				if(not isinstance(node, MultipoleTEAPOT)):
    					print "node = ", node.getName(), " type = ",\
    					node.getType(), " L = ", node.getLength()
    					orbitFinalize("Field Error: Wanted a Multipole node! Stop!")
    				klArr = node.getParam("kls")
    				for i in xrange(len(klArr)):
    					klArr[i] *= (1.0 + fracerr)
    				node.setParam("kls", klArr)
    		for node in lattice.getNodes()[nodeIndexi : Indexf]:
    			if(not isinstance(node, MultipoleTEAPOT)):
    				print "node = ", node.getName(), " type = ",\
    				node.getType(), " L = ", node.getLength()
    				orbitFinalize("Field Error: Wanted a Multipole node! Stop!")
    			klArr = node.getParam("kls")
    			for i in xrange(len(klArr)):
    				klArr[i] *= (1.0 + fracerr)
    			node.setParam("kls", klArr)
    	if(self.localDict["subtype"] == "QuadField"):
    		if(nodeIndexi > nodeIndexf):
    			for node in lattice.getNodes()[0 : nodeIndexf + 1]:
    				if(not isinstance(node, QuadTEAPOT)):
    					print "node = ", node.getName(), " type = ",\
    					node.getType(), " L = ", node.getLength()
    					orbitFinalize("Field Error: Wanted a Quadrupole node! Stop!")
    				kq = node.getParam("kq") * (1.0 + fracerr)
    				node.setParam("kq", kq)
    				klArr = node.getParam("kls")
    				for i in xrange(len(klArr)):
    					klArr[i] *= (1.0 + fracerr)
    				node.setParam("kls", klArr)
    		for node in lattice.getNodes()[nodeIndexi : Indexf]:
    			if(not isinstance(node, QuadTEAPOT)):
    				print "node = ", node.getName(), " type = ",\
    				node.getType(), " L = ", node.getLength()
    				orbitFinalize("Field Error: Wanted a Quadrupole node! Stop!")
    			kq = node.getParam("kq") * (1.0 + fracerr)
    			node.setParam("kq", kq)
    			klArr = node.getParam("kls")
    			for i in xrange(len(klArr)):
    				klArr[i] *= (1.0 + fracerr)
    			node.setParam("kls", klArr)
    	if(self.localDict["subtype"] == "BendField"):
    		if(nodeIndexi > nodeIndexf):
    			for node in lattice.getNodes()[0 : nodeIndexf + 1]:
    				if(not isinstance(node, BendTEAPOT)):
    					print "node = ", node.getName(), " type = ",\
    					node.getType(), " L = ", node.getLength()
    					orbitFinalize("Field Error: Wanted a Bend node! Stop!")
    				klArr = node.getParam("kls")
    				for i in xrange(len(klArr)):
    					klArr[i] *= (1.0 + fracerr)
    				node.setParam("kls", klArr)
    		for node in lattice.getNodes()[nodeIndexi : Indexf]:
    			if(not isinstance(node, BendTEAPOT)):
    				print "node = ", node.getName(), " type = ",\
    				node.getType(), " L = ", node.getLength()
    				orbitFinalize("Field Error: Wanted a Bend node! Stop!")
    			klArr = node.getParam("kls")
    			for i in xrange(len(klArr)):
    				klArr[i] *= (1.0 + fracerr)
    			node.setParam("kls", klArr)
    		nodei = lattice.getNodes()[nodeIndexi]
    		nodef = lattice.getNodes()[nodeIndexf]
    		drhoi = -nodei.getParam("rho") * fracerr / (1.0 + fracerr)
    		drhof = -nodef.getParam("rho") * fracerr / (1.0 + fracerr)
    		errori = bendfieldi(drhoi)
    		errorf = bendfieldf(drhof)
    		addErrorNodeAsChild_I(lattice, nodei, errori)
    		addErrorNodeAsChild_F(lattice, nodef, errorf)

    def AddBendDisplacementError(self):
    	"""
    		Adds a BendDisplacementError to nodes
    		between nodeIndexi and nodeIndexf.
    	"""
    	lattice = self.lattice
    	lattice.initialize()
    	nodeIndexi = self.nodeIndexi
    	nodeIndexf = self.nodeIndexf
    	Indexf = nodeIndexf + 1
    	sample = self.localDict["sample"]
    	multfrac = 1.0
    	if(sample == "Uniform"):
    		minimum = self.localDict["minimum"]
    		maximum = self.localDict["maximum"]
    		multfrac  = minimum + (maximum - minimum) * random()
    	if(sample == "Gaussian"):
    		mean   = self.localDict["mean"]
    		sigma  = self.localDict["sigma"]
    		multfrac = gauss(mean, sigma)
    	disp = self.localDict["disp"] * multfrac
    	theta = 0.0
    	if(nodeIndexi > nodeIndexf):
    		Indexf = len(lattice.getNodes())
    		for node in lattice.getNodes()[0 : nodeIndexf + 1]:
    			if(not isinstance(node, BendTEAPOT)):
    				print "Node = ", node.getName(), " type = ",\
    				node.getType(), " L = ", node.getLength()
    				orbitFinalize("Can't add BendError around a Straight! Stop!")
    			else:
    				theta += node.getParam("theta")
    	for node in lattice.getNodes()[nodeIndexi : Indexf]:
    		if(not isinstance(node, BendTEAPOT)):
    			print "Node = ", node.getName(), " type = ",\
    			node.getType(), " L = ", node.getLength()
    			orbitFinalize("Can't add BendError around a Straight! Stop!")
    		else:
    			theta += node.getParam("theta")
    	nodei = lattice.getNodes()[nodeIndexi]
    	nodef = lattice.getNodes()[nodeIndexf]
    	if(self.localDict["subtype"] == "XDisp"):
    		errori = benddisplacementxi(theta, disp)
    		errorf = benddisplacementxf(theta, disp)
    	if(self.localDict["subtype"] == "YDisp"):
    		errori = benddisplacementyi(disp)
    		errorf = benddisplacementyf(disp)
    	if(self.localDict["subtype"] == "LongDisp"):
    		errori = benddisplacementli(theta, disp)
    		errorf = benddisplacementlf(theta, disp)
    	addErrorNodeAsChild_I(lattice, nodei, errori)
    	addErrorNodeAsChild_F(lattice, nodef, errorf)

    def AddRotationError(self):
    	"""
    		Adds a RotationError to nodes
    		between nodeIndexi and nodeIndexf.
    	"""
    	lattice = self.lattice
    	lattice.initialize()
    	nodeIndexi = self.nodeIndexi
    	nodeIndexf = self.nodeIndexf
    	Indexf = nodeIndexf + 1
    	sample = self.localDict["sample"]
    	multfrac = 1.0
    	if(sample == "Uniform"):
    		minimum = self.localDict["minimum"]
    		maximum = self.localDict["maximum"]
    		multfrac  = minimum + (maximum - minimum) * random()
    	if(sample == "Gaussian"):
    		mean   = self.localDict["mean"]
    		sigma  = self.localDict["sigma"]
    		multfrac = gauss(mean, sigma)
    	angle = self.localDict["angle"] * multfrac
    	print "multfrac, angle = ", multfrac, angle
    	rhoi = 0.0
    	theta = 0.0
    	if(self.zf >= self.zi):
    		lengtherr = self.zf - self.zi
    	else:
    		lengtherr = lattice.getLength() - self.zf + self.zi
    	et = self.localDict["elementtype"]
    	rotype = self.localDict["subtype"]
    	if((et == "SBEN") or (et == "sbend") or (et == "rbend")):
    		if(nodeIndexi > nodeIndexf):
    			Indexf = len(lattice.getNodes())
    			for node in lattice.getNodes()[0 : nodeIndexf + 1]:
    				if(not isinstance(node, BendTEAPOT)):
    					print "Node = ", node.getName(), " type = ",\
    					node.getType(), " L = ", node.getLength()
    					orbitFinalize("Can't add BendError around a Straight! Stop!")
    				else:
    					rhoi = 1.0 / node.getParam("rho")
    					theta += node.getParam("theta")
    		for node in lattice.getNodes()[nodeIndexi : Indexf]:
    			if(not isinstance(node, BendTEAPOT)):
    				print "Node = ", node.getName(), " type = ",\
    				node.getType(), " L = ", node.getLength()
    				orbitFinalize("Can't add BendError around a Straight! Stop!")
    			else:
    				rhoi = 1.0 / node.getParam("rho")
    				theta += node.getParam("theta")
    	else:
    		if(nodeIndexi > nodeIndexf):
    			Indexf = len(lattice.getNodes())
    			for node in lattice.getNodes()[0 : nodeIndexf + 1]:
    				if(isinstance(node, BendTEAPOT)):
    					print "Node = ", node.getName(), " type = ",\
    					node.getType(), " L = ", node.getLength()
    					orbitFinalize("Can't add StraigtError around a Bend! Stop!")
    		for node in lattice.getNodes()[nodeIndexi : Indexf]:
    			if(isinstance(node, BendTEAPOT)):
    				print "Node = ", node.getName(), " type = ",\
    				node.getType(), " L = ", node.getLength()
    				orbitFinalize("Can't add StraigtError around a Bend! Stop!")
    	nodei = lattice.getNodes()[nodeIndexi]
    	nodef = lattice.getNodes()[nodeIndexf]
    	errori = rotationi(angle, rhoi, theta, lengtherr, et, rotype)
    	errorf = rotationf(angle, rhoi, theta, lengtherr, et, rotype)
    	addErrorNodeAsChild_I(lattice, nodei, errori)
    	addErrorNodeAsChild_F(lattice, nodef, errorf)


######################################################################


# Add a set of error nodes
class AddErrorSet():

    def __init__(self, lattice, positioni, positionf, setDict, paramsDict,\
	    name = "Error Set"):
        """
            Constructor. Adds the nodes.
        """
        et   = setDict["elementtype"]
        ringline = setDict["ringline"]
        errnodecandidates = []
        index = -1
        zf = 0.0
        for node in lattice.getNodes():
        	index += 1
        	zi = zf
        	zf += node.getLength()
        	status = []
        	status.append(index)
        	status.append(zi)
        	status.append(zf)
        	onoff = 0
        	if((et == "drift") and isinstance(node, DriftTEAPOT)):
        		onoff = 1
        	if((et == "kick")  and isinstance(node, KickTEAPOT)):
        		onoff = 1
        	if((et == "soln")  and isinstance(node, SolenoidTEAPOT)):
        		onoff = 1
        	if((et == "mult")  and isinstance(node, MultipoleTEAPOT)):
        		onoff = 1
        	if((et == "quad")  and isinstance(node, QuadTEAPOT)):
        		onoff = 1
        	if((et == "SBEN") or (et == "sbend") or (et == "rbend")):
        		if(isinstance(node, BendTEAPOT)): onoff = 1
        	status.append(onoff)
        	errnodecandidates.append(status)
        nodelist = []
        onoff = 0
        for index in range(0, len(errnodecandidates)):
        	if(errnodecandidates[index][3]) == 0:
        		if(onoff == 1):
        			status = []
        			status.append(istart)
        			status.append(istop)
        			status.append(zi)
        			status.append(zf)
        			nodelist.append(status)
        		onoff = 0
        	else:
        		if(onoff == 0):
        			istart = index
        			istop  = index
        			zi     = errnodecandidates[index][1]
        			zf     = errnodecandidates[index][2]
        		else:
        			istop  = index
        			zf     = errnodecandidates[index][2]
        		onoff = 1
        if(onoff == 1):
        	if(ringline == "line"):
        		status = []
        		status.append(istart)
        		status.append(istop)
        		status.append(zi)
        		status.append(zf)
        		nodelist.append(status)
        	elif(nodelist[0][0] == 0):
        		nodelist[0][0] = istart
        		nodelist[0][2] = zi
        tiny = 1.0e-07
        for index in range(0, len(nodelist)):
        	zi = nodelist[index][2] + tiny
        	zf = nodelist[index][3] - tiny
        	if(((positioni <= zf) and (zf <=  positionf)) or \
        		((positioni <= zi) and (zi <=  positionf)) or \
        		((zi <= positioni) and (positionf <= zf))):
        		AddErrorNode(lattice, zi, zf, paramsDict)
        		print "zi, zf = ", zi, zf

######################################################################


def FindNode(lattice, position):
	"""
	Finds node at position in lattice.
	"""
	lattice.initialize()
	nodeIndex = 0
	zf = 0.0
	yes = 1
	for node in lattice.getNodes():
		if(yes == 1):
			zi = zf
			zf += node.getLength()
		else:
			return (nodeIndex, zi, zf)
		if(position >= zf):
			nodeIndex += 1
		else:
			yes = 0
		if(position >= lattice.getLength()):
			nodeIndex -= 1
	return (nodeIndex, zi, zf)
