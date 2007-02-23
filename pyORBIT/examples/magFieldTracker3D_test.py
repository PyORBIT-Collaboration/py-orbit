import sys
from fieldtracker import MagneticFieldTracker3D
from fieldtracker import MultipoleExpansion3D
from bunch import Bunch

print "pyScript -> Begin"

# Initialize the MagneticFieldTracker3D class instance
test = MagneticFieldTracker3D()

# Initialize the Bunch to the proper initial conditions
bunch = Bunch()
bunch.addPartAttr("masscharge")
bunch.addParticle(0.0,0.0,0.0,0.0,0.0,0.0)
bunch.partAttrValue("masscharge",0,1,-1.0)     # this sets the charge

syncPart = bunch.getSyncParticle()
syncPart.x(0.0);
syncPart.y(0.0);
syncPart.z(0.0);
syncPart.momentum(1.0e3);  # sets px = py = 0.0 and pz to the argument

# Initialize the field
field = MultipoleExpansion3D("Chicane_combined_4c.dat")


test.track(0.0,2.0,0.0,0.0,0.0,0.0,0.0,bunch,field)
del bunch
del field
del test

print "pyScript -> Done!"




