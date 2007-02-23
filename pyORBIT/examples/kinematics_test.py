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

step = 0.00000001*m
nStep = int(0.01*m/step)

diff_max = 0.
p_max = 0.
e_max = 0.

for i in xrange(nStep):
	p0 = (i+1)*step + 10.0
	e = syncPart.momentumToEnergy(p0)
	p = syncPart.energyToMomentum(e)
	e1 = syncPart.momentumToEnergy(p)
	diff = abs(p0-p)/p
	if(diff > diff_max):
		diff_max = diff
		p_max = p
	diff = abs(e1-e)/e
	if(diff > diff_max):
		diff_max = diff
		e_max = e
	if(i%10000 == 0):
		print "i=",i," p=",p0

print "diff. max=",diff_max
print "p_max =",p_max
print "e_max =",e_max	

