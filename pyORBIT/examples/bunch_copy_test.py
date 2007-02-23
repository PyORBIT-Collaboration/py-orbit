import sys
from bunch import Bunch

#-----------------------------------------------------
#Copy bunch test
#-----------------------------------------------------

print "Start."

b = Bunch()

nParts = 20
for i in xrange(nParts):
	b.addParticle(0.1+i,0.2+i,0.3+i,0.4+i,0.5+i,0.6+i)
print "Initial bunch size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()

b.mass(1000.1111)
b.classicalRadius(222.222)
b.charge(10.5555)
b.macroSize(6666.666)
b.addPartAttr("macrosize")

syncPart = b.getSyncParticle()
syncPart.pz(100.0)
syncPart.rfFrequency(1.0123456789123456e+12)
syncPart.time(333.333)

b1 = Bunch()
#Copy b->b1
b.copyBunchTo(b1)

print "Old bunch ================================"
b.dumpBunch()


print "New bunch ================================"
b1.dumpBunch()

nParts = b1.getSize()
print "New bunch size=",b1.getSize()," total=",b1.getTotalCount()," capacity=",b1.getCapacity()
print "New bunch part attrs =",b1.getPartAttrNames()


sys.exit(1)

#memory leak test
count = 0
while(1 < 2):
	count += 1
	b2 = Bunch()
	b.copyBunchTo(b2)
	if((count % 1000) == 0):
		print "i=",count


print "Stop."
