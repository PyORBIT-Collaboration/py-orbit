import sys
import math
import posix

from lattice import AccLattice, AccLine, AccElement, AccActionsConatainer

lat = AccLattice("test_lattice")

line1 = AccLine("line1")
line2 = AccLine("line2")
line3 = AccLine("line3")
line4 = AccLine("line4")

elem1 = AccElement("el-1")
elem2 = AccElement("el-2")
elem3 = AccElement("el-3")
elem4 = AccElement("el-4")

elem5 = AccElement("el-5")

elem1.setLength(1.1)
elem2.setLength(2.1)
elem3.setLength(3.1)
elem4.setLength(4.1)

line1.addChildNode(elem1)
line2.addChildNode(elem2)
line3.addChildNode(elem3)
line4.addChildNode(elem4)

lat.addChildNode(line1)
lat.addChildNode(line2)
lat.addChildNode(line3)
lat.addChildNode(line4)
lat.addChildNode(elem5)

elem1_1 = AccElement("el-1-1")
elem1_1_1 = AccElement("el-1-1-1")
elem1_1_2 = AccElement("el-1-1-2")
elem1_1_3 = AccElement("el-1-1-3")
elem1_1_4 = AccElement("el-1-1-4")

elem1.addChildNode(elem1_1)
elem1_1.addChildNodeAtStart(elem1_1_1)
elem1_1.addChildNodeAtFinish(elem1_1_4)
elem1_1.addChildNode(elem1_1_2)
elem1_1.addChildNode(elem1_1_3)


elem1_2 = AccElement("el-1-2")
elem2.addChildNode(elem1_2)

acts = AccActionsConatainer()

def Blanks(n):
    s = ""
    for i in xrange(n):
        s += " "
    return s
    
nLevel = 0

def funcEntrance(paramsDict):
    global nLevel
    nLevel += 1
    node = paramsDict["node"]
    if(paramsDict.has_key("print") and paramsDict["print"] == True):
        print Blanks(nLevel),"ENTER level=",nLevel," node=",node.getName()

def funcExit(paramsDict):
    global nLevel
    node = paramsDict["node"]
    if(paramsDict.has_key("print") and paramsDict["print"] == True):
        print Blanks(nLevel),"EXIT  level=",nLevel," node=",node.getName()
    nLevel -= 1
    
def funcTrack(paramsDict):
    global nLevel
    node = paramsDict["node"]
    node.track(paramsDict)
    if(paramsDict.has_key("print") and paramsDict["print"] == True):
        print Blanks(nLevel),"TRACK through node =",node.getName()," level=",nLevel

acts.addEntranceAction(funcEntrance)
acts.addAction(funcTrack)
acts.addExitAction(funcExit)

lat.initialize()

print "Total length=",lat.getLength()

d = {"print":True}

lat.trackActions(acts,d)

#========Speed test==========================
count = 0
while(True):
    #lat.initialize()
    lat.trackActions(acts)
    if( count % 10000 == 0):
        print "i=",count, " time=",posix.times()[0]
    count += 1

print "====STOP==="

sys.exit(1)


