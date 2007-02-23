# AnalyticTest10.py
#
# This is run using the initial parameters that Jeff provided and a bunch that
#    he has also provided

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


######### REF PARTICLE (at 1GeV) ########
syncPart = bunch.getSyncParticle()
syncPart.x(0.0);
syncPart.y(0.0);
syncPart.z(0.0);
syncPart.momentum(1696.037);  # sets px = py = 0.0 and pz to the argument in (MeV)


######### READ IN BUNCH ########
inputFile = open('Matt_Bunch.dat','r')    # open a read-only file
numParticles = 10000                      # number of particles to read in

# read particles in 
for i in range(numParticles):
   temp = inputFile.readline()
   lst = temp.split()
   bunch.addParticle(float(lst[0]),float(lst[1])*1000.0,float(lst[2]),float(lst[3])*1000.0,
                     float(lst[4]),float(lst[5])*syncPart.momentum())
   bunch.partAttrValue("masscharge",i,1,+1.0)

# Initialize the field
field = MultipoleExpansion3D("Chicane_combined_4c.dat")

# track the particles from foil #1 to foil #2
x_0 = 40.0
y_0 = 23.0
alpha = 0.0
beta = 0.0   
gamma = 0.0
angles = test.track(0.309193,2.938333,x_0,y_0,alpha,beta,gamma,bunch,field,0)

# destructor
del bunch
del field
del test

print "pyScript -> Done!"
