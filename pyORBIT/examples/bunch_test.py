import sys
from bunch import Bunch

#-----------------------------------------------------
#Memory leak test - creates and destroys empty bunches
#-----------------------------------------------------

print "Start."

b = Bunch()

nParts = 100000
for i in xrange(nParts):
	b.addParticle(0.1+i,0.2+i,0.3+i,0.4+i,0.5+i,0.6+i)
	if(i%5000 == 0):
		print "i=",i," size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()

print "bunch size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()

for i in xrange(nParts/2):
	b.deleteParticleFast(i);

print "bunch size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()
b.compress()
print "bunch size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()

n = 0

while(1 > 0):
	for i in xrange(10000):
		b.addParticle(0.1+i,0.2+i,0.3+i,0.4+i,0.5+i,0.6+i)
	b.compress()
	for i in xrange(10000):
		b.deleteParticleFast(i)
	if(n%10 == 0):
		print "n=",n," bunch size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()
	b.compress()
	n = n + 1
	if(n%10 == 0):
		print "n=",n," bunch size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()


del b

print "Stop."
