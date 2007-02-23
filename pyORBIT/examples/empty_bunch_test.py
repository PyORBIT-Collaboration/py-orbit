import sys
from bunch import Bunch
from bunch import SyncParticle
#-----------------------------------------------------
#Memory leak test - creates and destroys empty bunches
#-----------------------------------------------------

print "Start."

b = Bunch()
print "bunch->cpp_ptr =",b.cpp_ptr
del b

#memory test
#the amount of memory should be the same
nIter = 1000000000
for i in xrange(nIter):
    b = Bunch()
    if i % 1000 == 0:
        print "i=",i
    del b

print "Stop."

