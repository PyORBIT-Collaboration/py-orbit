"""
Module. Includes classes for impedance accelerator nodes.
"""

import sys
import os
import math

# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import physical constants
from orbit.utils import consts

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode,\
     AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

#import impedance packages
from impedances import LImpedance
from impedances import TImpedance

#-----------------------------------------------------------------------------
# Node for LImpedance as function of node number
#-----------------------------------------------------------------------------

class LImpedance_Node(DriftTEAPOT):

    def __init__(self, phaseLength, nMacrosMin, nBins,\
	    name = "LImpedance node"):
        """
            Constructor. Creates LImpedance element.
        """
        DriftTEAPOT.__init__(self, name)
        self.limpedance = LImpedance(phaseLength, nMacrosMin, nBins)
        self.setType("limpedance node")
        self.setLength(0.0)

    def trackBunch(self, bunch):
        """
            The LImpedance-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        self.limpedance.trackBunch(bunch)	#track method goes here

    def track(self, paramsDict):
        """ 
            The LImpedance-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.limpedance.trackBunch(bunch)	#track method goes here

    def assignImpedance(self, py_cmplx_arr):
        self.limpedance.assignImpedance(py_cmplx_arr)
        print "Assigning the impedance array"

#-----------------------------------------------------------------------------
# Node for LImpedance as function of frequency
#-----------------------------------------------------------------------------

class FreqDep_LImpedance_Node(DriftTEAPOT):

    def __init__(self, phaseLength, nMacrosMin, nBins,\
	    bunch, impeDict,\
	    name = "freq. dep. LImpedance node"):
        """
            Constructor. Creates the FreqDep_LImpedance-teapot element.
        """
        DriftTEAPOT.__init__(self, name)
        self.limpedance = LImpedance(phaseLength, nMacrosMin, nBins)
        self.setType("freq. dep. limpedance node")
        self.setLength(0.0)
        self.phaseLength = phaseLength
        self.nBins = nBins
        self.localDict = impeDict
        self.freq_tuple = self.localDict["freqs"]
        self.freq_range = (len(self.freq_tuple) - 1)
        self.z_tuple = self.localDict["z_imp"]
        self.c = consts.speed_of_light
        myBeta = bunch.getSyncParticle().beta()
        myFreq = (myBeta * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = myFreq * (n + 1)
            checkFrequency(n, freq_mode, self.freq_tuple)
            z_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.z_tuple)
            Z.append(z_mode)
        self.limpedance.assignImpedance(Z)

    def trackBunch(self, bunch):
        """
            The FreqDep_LImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        myBeta = bunch.getSyncParticle().beta()
        myFreq = (myBeta * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = myFreq * (n + 1)
            checkFrequency(n, freq_mode, self.freq_tuple)
            z_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.z_tuple)
            Z.append(z_mode)
        self.limpedance.assignImpedance(Z)
        self.limpedance.trackBunch(bunch)

    def track(self, paramsDict):
        """
            The FreqDep_LImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        myBeta = bunch.getSyncParticle().beta()
        myFreq = (myBeta * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = myFreq * (n + 1)
            checkFrequency(n, freq_mode, self.freq_tuple)
            z_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.z_tuple)
            Z.append(z_mode)
        self.limpedance.assignImpedance(Z)
        self.limpedance.trackBunch(bunch)

#-----------------------------------------------------------------------------
# Node for LImpedance as function of beta and frequency
#-----------------------------------------------------------------------------

class BetFreqDep_LImpedance_Node(DriftTEAPOT):

    def __init__(self, phaseLength, nMacrosMin, nBins,\
	    bunch, impeDict,\
	    name = "freq. dep long sc node"):
        """
            Constructor. Creates the BetFreqDep_LImpedance-teapot element.
        """
        DriftTEAPOT.__init__(self, name)
        self.limpedance = LImpedance(phaseLength, nMacrosMin, nBins)
        self.setType("beta-freq. dep. limpedance node")
        self.setLength(0.0)
        self.phaseLength = phaseLength
        self.nBins = nBins
        self.localDict = impeDict
        self.bet_tuple = self.localDict["betas"]
        self.bet_range = (len(self.bet_tuple) - 1)
        self.freq_tuple = self.localDict["freqs"]
        self.freq_range = (len(self.freq_tuple) - 1)
        self.z_bf = self.localDict["z_imp"]
        self.c = consts.speed_of_light
        myBeta = bunch.getSyncParticle().beta()
        checkBeta(myBeta, self.bet_tuple)
        myFreq = (myBeta * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = myFreq * (n + 1)
            checkFrequency(n, freq_mode, self.freq_tuple)
            z_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.z_bf)
            Z.append(z_mode)
        self.limpedance.assignImpedance(Z)

    def trackBunch(self, bunch):
        """
            The BetFreqDep_LImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        myBeta = bunch.getSyncParticle().beta()
        checkBeta(myBeta, self.bet_tuple)
        myFreq = (myBeta * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = myFreq * (n + 1)
            checkFrequency(n, freq_mode, self.freq_tuple)
            z_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.z_bf)
            Z.append(z_mode)
        self.limpedance.assignImpedance(Z)
        self.limpedance.trackBunch(bunch)

    def track(self, paramsDict):
        """
            The BetFreqDep_LImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        myBeta = bunch.getSyncParticle().beta()
        checkBeta(myBeta, self.bet_tuple)
        myFreq = (myBeta * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = myFreq * (n + 1)
            checkFrequency(n, freq_mode, self.freq_tuple)
            z_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.z_bf)
            Z.append(z_mode)
        self.limpedance.assignImpedance(Z)
        self.limpedance.trackBunch(bunch)

#-----------------------------------------------------------------------------
# Node for TImpedance as function of node number
#-----------------------------------------------------------------------------

class TImpedance_Node(DriftTEAPOT):

    def __init__(self, phaseLength, nMacrosMin, nBins,\
	    useX, useY,\
	    name = "TImpedance node"):
        """
            Constructor. Creates TImpedance element.
        """
        DriftTEAPOT.__init__(self, name)
        self.timpedance = TImpedance(phaseLength, nMacrosMin, nBins,\
		useX, useY)
        self.setType("timpedance node")
        self.setLength(0.0)

    def trackBunch(self, bunch):
        """
            The TImpedance-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        self.timpedance.trackBunch(bunch)	#track method goes here

    def track(self, paramsDict):
        """ 
            The TImpedance-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.timpedance.trackBunch(bunch)	#track method goes here

    def assignLatFuncs(self, qX, alphaX, betaX, qY, alphaY, betaY):
        self.timpedance.assignLatFuncs(qX, alphaX, betaX, qY, alphaY, betaY)
        print "Assigning lattice functions at timpedance location"

    def assignImpedance(self, XorY, py_cmplx_arrp, py_cmplx_arrm):
        self.timpedance.assignImpedance(XorY, py_cmplx_arrp, py_cmplx_arrm)
        print "Assigning the horizontal or vertical impedance arrays"

#-----------------------------------------------------------------------------
# Node for TImpedance as function of frequency
#-----------------------------------------------------------------------------

class FreqDep_TImpedance_Node(DriftTEAPOT):

    def __init__(self, phaseLength, nMacrosMin, nBins,\
	    useX, useY,\
	    bunch, impeDict,\
	    name = "freq. dep. TImpedance node"):
        """
            Constructor. Creates the FreqDep_TImpedance-teapot element.
        """
        DriftTEAPOT.__init__(self, name)
        self.timpedance = TImpedance(phaseLength, nMacrosMin, nBins,\
		useX, useY)
        self.setType("freq. dep. timpedance node")
        self.setLength(0.0)
        self.phaseLength = phaseLength
        self.nBins = nBins
        self.localDict = impeDict
        self.freq_tuple = self.localDict["freqs"]
        self.freq_range = (len(self.freq_tuple) - 1)
        self.c = consts.speed_of_light
        myBeta = bunch.getSyncParticle().beta()
        myFreq = (myBeta * self.c) / self.phaseLength
	if(useX != 0):
		self.zp_tuple = self.localDict["zxp_imp"]
		self.zm_tuple = self.localDict["zxm_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zp_tuple)
			zm_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zm_tuple)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("X", Zp, Zm)
	if(useY != 0):
		self.zp_tuple = self.localDict["zyp_imp"]
		self.zm_tuple = self.localDict["zym_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zp_tuple)
			zm_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zm_tuple)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("Y", Zp, Zm)

    def trackBunch(self, bunch):
        """
            The FreqDep_TImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        myBeta = bunch.getSyncParticle().beta()
        myFreq = (myBeta * self.c) / self.phaseLength
	if(useX != 0):
		self.zp_tuple = self.localDict["zxp_imp"]
		self.zm_tuple = self.localDict["zxm_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zp_tuple)
			zm_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zm_tuple)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("X", Zp, Zm)
	if(useY != 0):
		self.zp_tuple = self.localDict["zyp_imp"]
		self.zm_tuple = self.localDict["zym_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zp_tuple)
			zm_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zm_tuple)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("Y", Zp, Zm)
        self.timpedance.trackBunch(bunch)

    def track(self, paramsDict):
        """
            The FreqDep_TImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        myBeta = bunch.getSyncParticle().beta()
	if(useX != 0):
		self.zp_tuple = self.localDict["zxp_imp"]
		self.zm_tuple = self.localDict["zxm_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zp_tuple)
			zm_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zm_tuple)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("X", Zp, Zm)
	if(useY != 0):
		self.zp_tuple = self.localDict["zyp_imp"]
		self.zm_tuple = self.localDict["zym_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zp_tuple)
			zm_mode = interp(freq_mode, self.freq_range,\
                            self.freq_tuple, self.zm_tuple)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("Y", Zp, Zm)
        self.timpedance.trackBunch(bunch)

#-----------------------------------------------------------------------------
# Node for TImpedance as function of beta and frequency
#-----------------------------------------------------------------------------

class BetFreqDep_TImpedance_Node(DriftTEAPOT):

    def __init__(self, phaseLength, nMacrosMin, nBins,\
	    useX, useY,\
	    bunch, impeDict,\
	    name = "freq. dep long sc node"):
        """
            Constructor. Creates the BetFreqDep_TImpedance-teapot element.
        """
        DriftTEAPOT.__init__(self, name)
        self.timpedance = TImpedance(phaseLength, nMacrosMin, nBins,\
		useX, useY)
        self.setType("beta-freq. dep. timpedance node")
        self.setLength(0.0)
        self.phaseLength = phaseLength
        self.nBins = nBins
        self.localDict = impeDict
        self.bet_tuple = self.localDict["betas"]
        self.bet_range = (len(self.bet_tuple) - 1)
        self.freq_tuple = self.localDict["freqs"]
        self.freq_range = (len(self.freq_tuple) - 1)
        self.c = consts.speed_of_light
        myBeta = bunch.getSyncParticle().beta()
        checkBeta(myBeta, self.bet_tuple)
        myFreq = (myBeta * self.c) / self.phaseLength
	if(useX != 0):
		self.zp_bf = self.localDict["zxp_imp"]
		self.zm_bf = self.localDict["zxm_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zp_bf)
			zm_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zm_bf)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("X", Zp, Zm)
	if(useY != 0):
		self.zp_bf = self.localDict["zyp_imp"]
		self.zm_bf = self.localDict["zym_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zp_bf)
			zm_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zm_bf)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("Y", Zp, Zm)

    def trackBunch(self, bunch):
        """
            The BetFreqDep_TImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        myBeta = bunch.getSyncParticle().beta()
        checkBeta(myBeta, self.bet_tuple)
        myFreq = (myBeta * self.c) / self.phaseLength
	if(useX != 0):
		self.zp_bf = self.localDict["zxp_imp"]
		self.zm_bf = self.localDict["zxm_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zp_bf)
			zm_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zm_bf)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("X", Zp, Zm)
	if(useY != 0):
		self.zp_bf = self.localDict["zyp_imp"]
		self.zm_bf = self.localDict["zym_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zp_bf)
			zm_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zm_bf)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("Y", Zp, Zm)
        self.timpedance.trackBunch(bunch)

    def track(self, paramsDict):
        """
            The BetFreqDep_TImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        myBeta = bunch.getSyncParticle().beta()
        checkBeta(myBeta, self.bet_tuple)
        myFreq = (myBeta * self.c) / self.phaseLength
	if(useX != 0):
		self.zp_bf = self.localDict["zxp_imp"]
		self.zm_bf = self.localDict["zxm_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zp_bf)
			zm_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zm_bf)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("X", Zp, Zm)
	if(useY != 0):
		self.zp_bf = self.localDict["zyp_imp"]
		self.zm_bf = self.localDict["zym_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2 - 1):
			freq_mode = myFreq * (n + 1)
			checkFrequency(n, freq_mode, self.freq_tuple)
			zp_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zp_bf)
			zm_mode = bilinterp(myBeta, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.zm_bf)
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("Y", Zp, Zm)
        self.timpedance.trackBunch(bunch)

#-----------------------------------------------------------------------------
# Methods used by LImpedance and TImpedance classes
#-----------------------------------------------------------------------------

def interp(x, n_tuple, x_tuple, y_tuple):
    """
        Linear interpolation: Given n-tuple + 1 points,
        x_tuple and y_tuple, routine finds y = y_tuple
        at x in x_tuple. Assumes x_tuple is increasing array.
    """
    if x < x_tuple[0]:
        print "*******************--Warning--*********************"
        print "interp: x < Min(x_tuple)"
        y = y_tuple[0]
        return y
    if x > x_tuple[n_tuple]:
        print "*******************--Warning--*********************"
        print "interp: x > Max(x_tuple)"
        y = y_tuple[n_tuple]
        return y
    dxp = x - x_tuple[0]
    for n in range(n_tuple):
        dxm = dxp
        dxp = x - x_tuple[n + 1]
        dxmp = dxm * dxp
        if dxmp <= 0:
            break
    y = (-dxp * y_tuple[n] + dxm * y_tuple[n + 1]) /\
		(dxm - dxp)
    return y


def bilinterp(x, y, nx_tuple, ny_tuple, x_tuple, y_tuple, fxy):
    """
        Bilinear interpolation: Given nx-tuple + 1 x-points,
        ny-tuple + 1 y-points, x_tuple and y_tuple,
        routine finds f(x, y) = fxy at (x, y) in (x_tuple, y_tuple).
        Assumes x_tuple and y_tuple are increasing arrays.
    """
    f_tuple = []
    if x < x_tuple[0]:
        print "*******************--Warning--*********************"
        print "bilinterp: x < Min(x_tuple)"
        for ny in range(ny_tuple + 1):
            vf = fxy[0][ny]
            f_tuple.append(vf)
    elif x > x_tuple[nx_tuple]:
        print "*******************--Warning--*********************"
        print "bilinterp: x > Max(x_tuple)"
        for ny in range(ny_tuple + 1):
            vf = fxy[x_tuple][ny]
            f_tuple.append(vf)
    else:
        dxp = x - x_tuple[0]
        for nx in range(nx_tuple):
            dxm = dxp
            dxp = x - x_tuple[nx + 1]
            dxmp = dxm * dxp
            if dxmp <= 0:
                break
	for ny in range(ny_tuple + 1):
		vf = (-dxp * fxy[nx][ny] + dxm * fxy[nx + 1][ny]) /\
		(dxm - dxp)
		f_tuple.append(vf)
    f = interp(y, ny_tuple, y_tuple, f_tuple)
    return f

def checkFrequency(n, myFreq, freq_tuple):
    """
        This function returns a warning that if the frequency is
        outside the range of the provided values, the corresponding
        impedance will be set appropriately to that at the maximum
        or minimum frequency.
    """
    if myFreq < min(freq_tuple):
        print "******************--Warning--*********************"
        print "minError in bin %d: \n Your frequency is %f. \n This is below the minimum provided frequency of %f." % (n, myFreq, min(freq_tuple))
        print " -The interpolation function will use the minimum provided frequency."
        print "**************************************************"
    elif myFreq > max(freq_tuple):
        print "******************--Warning--*********************"
        print "maxerror in bin %d: \n Your frequency is %f. \n This is above the maximum provided frequency of %f." % (n, myFreq, max(freq_tuple))
        print " -The interpolation function will use the maximum provided frequency."
        print "**************************************************"

def checkBeta(myBeta, bet_tuple):
    """
        This function returns a warning that if the beta is
        outside the range of the provided values, the corresponding
        impedance will be set appropriately to that at the maximum
        or minimum beta.
    """
    if myBeta < min(bet_tuple):
        print "******************--Warning--*********************"
        print "minError: \n Your beta is %f. \n This is below the minimum provided beta of %f." % (myBeta, min(bet_tuple))
        print " -The interpolation function will use the minimum provided beta."
        print "***************************************************"
    elif myBeta > max(bet_tuple):
        print "******************--Warning--**********************"
        print "maxerror: \n Your beta is %f. \n This is above the maximum provided beta of %f." % (n, myBeta, max(bet_tuple))
        print " -The interpolation function will use the maximum provided beta."
        print "***************************************************"
