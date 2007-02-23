import sys
from fieldtracker import MultipoleExpansion3D

multExp3D = MultipoleExpansion3D("Chicane_combined_4c.dat")

fieldTuple = multExp3D.getMagneticField(0.,1.,2.,3.)

print "fieldTuple=",fieldTuple

print "Done!"
sys.exit(1)

i = 0

while(True):
	multExp3D = MultipoleExpansion3D()
	fieldTuple = multExp3D.getMagneticField(0.,1.,2.,3.)
	#print "fieldTuple=",fieldTuple
	del multExp3D
	i = i + 1
	if(i % 1000 == 0):
		print "i=",i


print "Done!"


