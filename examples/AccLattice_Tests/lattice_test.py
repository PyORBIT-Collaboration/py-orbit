import sys
import math
import posix

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

lattice = AccLattice("test_lattice")

elem1 = AccNode("el-1")
elem2 = AccNode("el-2")
elem3 = AccNode("el-3")

elem1.setLength(1.1)
elem2.setLength(2.1)
elem3.setLength(3.1)

lattice.addNode(elem1)
lattice.addNode(elem2)
lattice.addNode(elem3)

elem1_1 = AccNode("el-1-1")
elem1_1.setnParts(2)
elem1_1_1 = AccNode("el-1-1-1")
elem1_1_2 = AccNode("el-1-1-2")
elem1_1_3 = AccNode("el-1-1-3")
elem1_1_4 = AccNode("el-1-1-4")

elem1.addChildNode(elem1_1,AccNode.ENTRANCE)
elem1_1.addChildNode(elem1_1_1,AccNode.ENTRANCE)
elem1_1.addChildNode(elem1_1_2,AccNode.BODY,0)
elem1_1.addChildNode(elem1_1_3,AccNode.BODY,1)
elem1_1.addChildNode(elem1_1_4,AccNode.EXIT)


elem1_2 = AccNode("el-1-2")
elem2.addChildNode(elem1_2,AccNode.EXIT)

acts = AccActionsContainer()

def Blanks(n):
    s = ""
    for i in xrange(n):
        s += " "
    return s

nLevel = [0]
nElems = [0]

def funcEntrance(paramsDict):
    nLevel[0] += 1
    node = paramsDict["node"]
    if(paramsDict.has_key("print") and paramsDict["print"] == True):
			print Blanks(nLevel[0]),"ENTER level=",nLevel[0]," node=",node.getName()
			nElems[0] += 1

def funcExit(paramsDict):
    node = paramsDict["node"]
    if(paramsDict.has_key("print") and paramsDict["print"] == True):
        print Blanks(nLevel[0]),"EXIT  level=",nLevel[0]," node=",node.getName()
    nLevel[0] -= 1

def funcTrack(paramsDict):
    node = paramsDict["node"]
    if(paramsDict.has_key("print") and paramsDict["print"] == True):
        print Blanks(nLevel[0]),"BODY TRACK through node =",node.getName()," level=",nLevel[0]

acts.addAction(funcEntrance,AccActionsContainer.ENTRANCE)
acts.addAction(funcTrack,AccActionsContainer.BODY)
acts.addAction(funcExit,AccActionsContainer.EXIT)

lattice.initialize()

print "Total length=",lattice.getLength()

nodes = lattice.getNodes()
for node in nodes:
	print "node=",node.getName()," s start,stop = %4.3f %4.3f "%lattice.getNodePositionsDict()[node]


d = {"print":True}

lattice.trackActions(acts,d)

print "Total number of nodes=",nElems[0]
#========Speed test==========================
count = 1
while(count <100000):
    #lattice.initialize()
    lattice.trackActions(acts)
    if( count % 10000 == 0):
        print "i=",count, " time= %9.8f "%(posix.times()[0]/count)
    count += 1

print "====STOP==="

sys.exit(0)


