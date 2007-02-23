from bunch import Bunch

#-----------------------------------------------------
#Documentation test for Bunch class
#-----------------------------------------------------

print "Start."

m_list = dir(Bunch)

print "============Bunch class documentation==START=========="

for m in m_list:
    print "method=",str(m)," DocStr=",getattr(Bunch,m).__doc__

print "============Bunch class documentation===STOP ========="

print "Stop."
