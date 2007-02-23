import teapot
from bunch import Bunch

from lattice import AccLattice, AccLine, AccElement, AccActionsConatainer

print "Start."

b = Bunch()
b.addParticle(1.0e-3,0.0,0.0,0.0,0.0,0.0)
b.addParticle(0.0,1.0e-3,0.0,0.0,0.0,0.0)
b.addParticle(0.0,0.0,1.0e-3,0.0,0.0,0.0)
b.addParticle(0.0,0.0,0.0,1.0e-3,0.0,0.0)
b.addParticle(0.0,0.0,0.0,0.0,1.0,0.0)
b.addParticle(0.0,0.0,0.0,0.0,0.0,1.0)
b.compress()

syncPart = b.getSyncParticle()
energy = 1.0e+3                          #energy in MeV
#p = syncPart.energyToMomentum(energy)
#syncPart.pz(p)
syncPart.kinEnergy(energy)

latt = AccLattice("test_lattice")

elem1 = teapot.DriftTEAPOT("drift1")
elem2 = teapot.QuadTEAPOT("quad1")
elem3 = teapot.QuadTEAPOT("quad2")
elem4 = teapot.BendTEAPOT("bend1")
elem5 = teapot.BendTEAPOT("bend2")
elem6 = teapot.MultipoleTEAPOT("sextupole")

latt.addChildNode(elem1)
latt.addChildNode(elem2)
latt.addChildNode(elem3)
latt.addChildNode(elem4)
latt.addChildNode(elem5)
latt.addChildNode(elem6)

#-----------------------------
# Set TEAPOT nodes parameters
#-----------------------------
elem1.setLength(0.2)
elem2.setLength(0.3)
elem3.setLength(0.4)
elem4.setLength(0.5)
elem5.setLength(0.6)
elem6.setLength(0.7)

elem2.setnParts(5)
elem2.addParam("kq",-0.7)

elem3.setnParts(5)
elem3.addParam("kq",+0.7)

elem4.setnParts(11)
elem4.addParam("theta",+0.1)
elem4.addParam("ea1",+0.01)
elem4.addParam("ea2",+0.02)

elem5.setnParts(11)
elem5.addParam("theta",+0.2)
elem5.addParam("ea1",+0.01)
elem5.addParam("ea2",+0.02)

elem6.setnParts(5)
elem6.getParam("poles").append(2)
elem6.getParam("kls").append(0.28)
elem6.getParam("skews").append(0)

latt.initialize()

#set phase scale
ringLength = 2.7
c = 2.99792458e+8
syncPart.rfFrequency((syncPart.beta()*c)/ringLength)

print "==============BEFORE============================"
#b.dumpBunch()
print "=========================================="


#=====STOP example ============
def stopAction(paramsDict):
    node = paramsDict["node"]
    actions = paramsDict["actions"]
    if(node == elem4):
        actions.setShouldStop(True)
    
accContainer = AccActionsConatainer()
accContainer.addEntranceAction(stopAction)

latt.trackBunch(b)
#latt.trackBunch(b,accContainer)

print "=============AFTER============================="
b.dumpBunch()
print "=========================================="

print "lattice length=",latt.getLength()
print "beta=",b.getSyncParticle().beta()
print "TEAPOT time[sec]=",b.getSyncParticle().time()
print "SIMPLE time[sec]=",latt.getLength()/(b.getSyncParticle().beta()*2.99792458e+8)
print "Stop."

#========================================================
#         TEST OUTPUT
#========================================================
#==============BEFORE============================
#==========================================
#=============AFTER=============================
#% PARTICLE_ATTRIBUTES_CONTROLLERS_NAMES
#% BUNCH_ATTRIBUTE_DOUBLE charge   1
#% BUNCH_ATTRIBUTE_DOUBLE classical_radius   1.5347e-18
#% BUNCH_ATTRIBUTE_DOUBLE macro_size   0
#% BUNCH_ATTRIBUTE_DOUBLE mass   938.272
#%  SYNC_PART_COORDS 0 0 0  x, y, z positions in [m]
#%  SYNC_PART_MOMENTUM 0 0 1696.037912  px, py, pz momentum component in MeV/c
#%  info only: energy of the synchronous particle [MeV] = 1000
#%  info only: momentum of the synchronous particle [MeV/c] = 1696.037912
#%  info only: beta=v/c of the synchronous particle = 0.8750256155
#%  info only: gamma=1/sqrt(1-(v/c)**2) of the synchronous particle = 2.065788684
#%  SYNC_PART_TIME 0  time in [sec]
#%  SYNC_PART_RF_FREQUENCY 97157807.44  rf frequency in [Hz]
#% x[m] px[rad] y[m] py[rad] (z_or_phi or pz_or_dE)
#0.0008284196 -0.00015325622 0 0 0.00069582754 -2.1657666e-13
#0.002363176 0.00076945303 0 0 0.0010010072 -2.1665711e-13
#2.3859816e-08 1.3257048e-07 0.0010236879 4.1172715e-05 1.0454023e-07 -2.1657666e-13
#-4.7883175e-07 6.8646526e-07 0.0028841539 0.0010928603 3.6930219e-06 -2.1673757e-13
#-2.7891822e-17 -2.6327478e-17 0 0 1 -2.1657666e-13
#0.00023639448 0.00019942888 0 0 -0.00096634215 1

