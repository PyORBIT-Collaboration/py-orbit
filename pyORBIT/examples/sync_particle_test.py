import sys
from bunch import Bunch
from bunch import SyncParticle
#-----------------------------------------------------
#Memory leak test - creates and destroys empty bunches
#-----------------------------------------------------

print "Start."

b = Bunch()

syncPart = b.getSyncParticle()
m = syncPart.mass()
print "m=",m
syncPart.x(0.1)
syncPart.y(0.2)
syncPart.z(0.3)
print "x=",syncPart.x()
print "y=",syncPart.y()
print "z=",syncPart.z()
syncPart.px(1.1)
syncPart.py(1.2)
syncPart.pz(1.3)
print "px=",syncPart.px()
print "py=",syncPart.py()
print "pz=",syncPart.pz()

p = syncPart.momentum()
print "p=",p

#memory test
#the amount of memory should be the same
nIter = 100000000
for i in xrange(nIter):
	m = syncPart.mass()
	p = syncPart.momentum()
	syncPart.px(syncPart.px())
	syncPart.py(syncPart.py())
	syncPart.pz(syncPart.pz())
	syncPart.x(syncPart.x())
	syncPart.y(syncPart.y())
	syncPart.z(syncPart.z())
	beta = syncPart.beta()
	g = syncPart.gamma()
	e = syncPart.kinEnergy()
	p = syncPart.energyToMomentum(e)
	e = syncPart.momentumToEnergy(p)
	f = syncPart.rfFrequency()
	syncPart.rfFrequency(f)
	if i % 1000 == 0:
		print "i=",i

print "Stop."

