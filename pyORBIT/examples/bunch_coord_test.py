import sys
from bunch import Bunch

#-----------------------------------------------------
#Coordinates function test
#-----------------------------------------------------

print "Start."

b = Bunch()

nParts = 20
for i in xrange(nParts):
	b.addParticle(0.1+i,0.2+i,0.3+i,0.4+i,0.5+i,0.6+i)
print "======before compress==========="
print "bunch size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()
print "======after compress==========="
b.compress()
print "bunch size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()
print "======coords test==========="
for i in xrange(nParts):
	print "i=",i," x=",b.x(i)," y=",b.y(i)," z=",b.z(i)," phi=",b.phi(i),\
	" px=",b.px(i)," py=",b.py(i)," pz=",b.pz(i)," dE=",b.dE(i)," flag=",b.flag(i)

nParts = b.getSize()

print "======coords test==========="
for i in xrange(nParts):
	b.x(i,2*b.x(i))
	b.y(i,2*b.y(i))
	b.z(i,2*b.z(i))
	b.px(i,2*b.px(i))
	b.py(i,2*b.py(i))
	b.pz(i,2*b.pz(i))
	#b.flag(i,0)
	print "i=",i," x=",b.x(i)," y=",b.y(i)," z=",b.z(i)," phi=",b.phi(i),\
	" px=",b.px(i)," py=",b.py(i)," pz=",b.pz(i)," dE=",b.dE(i)," flag=",b.flag(i)
print "======after deleteAllParticles()==========="
b.deleteAllParticles()
print "bunch size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()
print "======after compress==========="
b.compress()
print "bunch size=",b.getSize()," total=",b.getTotalCount()," capacity=",b.getCapacity()
b.dumpBunch("bunch_test.dat")

#populate bunch again
for i in xrange(nParts):
	b.addParticle(0.1+i,0.2+i,0.3+i,0.4+i,0.5+i,0.6+i)
b.compress()

print "============================================="
print "mass=",b.mass()
b.mass(1000.1111)
print "mass=",b.mass()

print "R=",b.classicalRadius()
b.classicalRadius(222.222)
print "R=",b.classicalRadius()

print "charge=",b.charge()
b.charge(10.5555)
print "charge=",b.charge()

print "macroSize=",b.macroSize()
b.macroSize(6666.666)
print "macroSize=",b.macroSize()

b.deleteAllParticles()
b.dumpBunch()
b.readBunch("bunch_test.dat",5)
print "============================================="
b.dumpBunch()

b.bunchAttrDouble("mass",222.333)
print "bunch attr name=mass value=",b.bunchAttrDouble("mass")

b.bunchAttrDouble("aaa",222)
b.bunchAttrDouble("bbb",333)
b.bunchAttrInt("aaa",111)
b.bunchAttrInt("bbb",444)
print "bunch attr dbl name=aaa value=",b.bunchAttrDouble("aaa")
print "bunch attr dbl name=bbb value=",b.bunchAttrDouble("bbb")
print "bunch attr int name=aaa value=",b.bunchAttrInt("aaa")
print "bunch attr int name=bbb value=",b.bunchAttrInt("bbb")

print "bunch double attrs=",b.bunchAttrDoubleNames()
print "bunch int attrs=",b.bunchAttrIntNames()

b.dumpBunch("bunch_test.dat")

b1 = Bunch()
print "bunch b1 double attrs=",b1.bunchAttrDoubleNames()
print "bunch b1 int attrs=",b1.bunchAttrIntNames()
b1.initBunchAttr("bunch_test.dat")
print "bunch b1 double attrs=",b1.bunchAttrDoubleNames()
print "bunch b1 int attrs=",b1.bunchAttrIntNames()

print "bunch b1 ccc int attrs has=",b1.hasBunchAttrInt("ccc")
print "bunch b1 aaa int attrs has=",b1.hasBunchAttrInt("aaa")
print "bunch b1 ccc double attrs has=",b1.hasBunchAttrDouble("ccc")
print "bunch b1 aaa double attrs has=",b1.hasBunchAttrDouble("aaa")


print "bunch possible part attrs =",b.getPossiblePartAttrNames()
print "bunch part attrs =",b.getPartAttrNames()
b.addPartAttr("macrosize")
print "bunch part attrs =",b.getPartAttrNames()
b.dumpBunch()
b.removePartAttr("macrosize")
print "bunch part attrs =",b.getPartAttrNames()
b.dumpBunch()
b.addPartAttr("macrosize")
b.removePartAttr("macrosize")
for i in xrange(100):
	b.addParticle(0.1+i,0.2+i,0.3+i,0.4+i,0.5+i,0.6+i)
nParts = b.getSize()
print "number of macroparticles=",nParts
count = 0
while(1<2):
	b.addPartAttr("macrosize")
	for i in xrange(nParts):
		b.partAttrValue("macrosize",i,0,i+0.1)
	b.removePartAttr("macrosize")
	count = count + 1
	if(count%1000 == 0):
		print "count =",count
b.dumpBunch()
print "====before delete===="
del b
del b1
print "====after delete===="

print "Stop."
