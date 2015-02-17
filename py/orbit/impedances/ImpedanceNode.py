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
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
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
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
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
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
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
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = bilinterp(BetaRel, freq_mode,\
		    self.bet_range, self.freq_range,\
		    self.bet_tuple, self.freq_tuple,\
		    self.z_bf)
            Z.append(z_mode)
        self.limpedance.assignImpedance(Z)

    def trackBunch(self, bunch):
        """
            The BetFreqDep_LImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = bilinterp(BetaRel, freq_mode,\
		    self.bet_range, self.freq_range,\
		    self.bet_tuple, self.freq_tuple,\
		    self.z_bf)
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
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = bilinterp(BetaRel, freq_mode,\
		    self.bet_range, self.freq_range,\
		    self.bet_tuple, self.freq_tuple,\
		    self.z_bf)
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

    def assignImpedance(self, XorY, py_cmplx_arrp, py_cmplx_arrm):
        self.timpedance.assignImpedance(XorY, py_cmplx_arrp, py_cmplx_arrm)

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
	self.useX = useX
	self.useY = useY
        self.localDict = impeDict
        self.freq_tuple = self.localDict["freqs"]
        self.freq_range = (len(self.freq_tuple) - 1)
        self.c = consts.speed_of_light
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
	self.calcImpedance(Freq0, self.qX, self.qY)

    def trackBunch(self, bunch):
        """
            The FreqDep_TImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
	self.calcImpedance(Freq0, self.qX, self.qY)
        self.timpedance.trackBunch(bunch)

    def track(self, paramsDict):
        """
            The FreqDep_TImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        bunch = paramsDict["bunch"]
        length = self.getLength(self.getActivePartIndex())
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
	self.calcImpedance(Freq0, self.qX, self.qY)
        self.timpedance.trackBunch(bunch)

    def assignLatFuncs(self, qX, alphaX, betaX, qY, alphaY, betaY):
	self.qX = qX
	self.qY = qY
        self.timpedance.assignLatFuncs(qX, alphaX, betaX, qY, alphaY, betaY)

    def calcImpedance(self, Freq0, qX, qY):
	if(self.useX != 0):
		z_tuple = self.localDict["zx_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2):
			Freq_p = Freq0 * (n + qX)
			Freq_m = Freq0 * (n - qX)
			sign_m = 1.0
			if(Freq_m < 0):
				Freq_m = -Freq_m
				sign_m = -1.0
			zp_mode = interp(Freq_p, self.freq_range,\
                            self.freq_tuple, z_tuple)
			zm_mode = interp(Freq_m, self.freq_range,\
                            self.freq_tuple, z_tuple)
			zm_mode.real = sign_m * zm_mode.real
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		zp_mode = 0.0 + 0.0j
		zm_mode = 0.0 + 0.0j
		Zp.append(zp_mode)
		Zm.append(zm_mode)
		for n in range(1, self.nBins / 2):
			zp_mode.real = -Zm[self.nBins / 2 - n].real
			zp_mode.imag =  Zm[self.nBins / 2 - n].imag
			zm_mode.real = -Zp[self.nBins / 2 - n].real
			zp_mode.imag =  Zp[self.nBins / 2 - n].imag
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("X", Zp, Zm)
	if(self.useY != 0):
		z_tuple = self.localDict["zy_imp"]
		Zp = []
		Zm = []
		for n in range(self.nBins / 2):
			Freq_p = Freq0 * (n + qY)
			Freq_m = Freq0 * (n - qY)
			sign_m = 1.0
			if(Freq_m < 0):
				Freq_m = -Freq_m
				sign_m = -1.0
			zp_mode = interp(Freq_p, self.freq_range,\
                            self.freq_tuple, z_tuple)
			zm_mode = interp(Freq_m, self.freq_range,\
                            self.freq_tuple, z_tuple)
			zm_mode.real = sign_m * zm_mode.real
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		zp_mode = 0.0 + 0.0j
		zm_mode = 0.0 + 0.0j
		Zp.append(zp_mode)
		Zm.append(zm_mode)
		for n in range(1, self.nBins / 2):
			zp_mode.real = -Zm[self.nBins / 2 - n].real
			zp_mode.imag =  Zm[self.nBins / 2 - n].imag
			zm_mode.real = -Zp[self.nBins / 2 - n].real
			zp_mode.imag =  Zp[self.nBins / 2 - n].imag
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("Y", Zp, Zm)

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
	self.useX = useX
	self.useY = useY
        self.localDict = impeDict
        self.bet_tuple = self.localDict["betas"]
        self.bet_range = (len(self.bet_tuple) - 1)
        self.freq_tuple = self.localDict["freqs"]
        self.freq_range = (len(self.freq_tuple) - 1)
	self.zx_bf = self.localDict["zx_imp"]
	self.zy_bf = self.localDict["zy_imp"]
        self.c = consts.speed_of_light
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
	self.calcImpedance(BetaRel, Freq0, self.qX, self.qY)

    def trackBunch(self, bunch):
        """
            The BetFreqDep_TImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
	self.calcImpedance(BetaRel, Freq0, self.qX, self.qY)
        self.timpedance.trackBunch(bunch)

    def track(self, paramsDict):
        """
            The BetFreqDep_TImpedance-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        bunch = paramsDict["bunch"]
        length = self.getLength(self.getActivePartIndex())
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
	self.calcImpedance(BetaRel, Freq0, self.qX, self.qY)
        self.timpedance.trackBunch(bunch)

    def assignLatFuncs(self, qX, alphaX, betaX, qY, alphaY, betaY):
	self.qX = qX
	self.qY = qY
        self.timpedance.assignLatFuncs(qX, alphaX, betaX, qY, alphaY, betaY)

    def calcImpedance(self, BetaRel, Freq0, qX, qY):
	if(self.useX != 0):
		Zp = []
		Zm = []
		for n in range(self.nBins / 2):
			Freq_p = Freq0 * (n + qX)
			Freq_m = Freq0 * (n - qX)
			sign_m = 1.0
			if(Freq_m < 0):
				Freq_m = -Freq_m
				sign_m = -1.0
			zp_mode = bilinterp(BetaRel, Freq_p,\
				self.bet_range, self.freq_range,\
				self.bet_tuple, self.freq_tuple,\
				self.zx_bf)
			zm_mode = bilinterp(BetaRel, Freq_m,\
				self.bet_range, self.freq_range,\
				self.bet_tuple, self.freq_tuple,\
				self.zx_bf)
			zm_mode.real = sign_m * zm_mode.real
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		zp_mode = 0.0 + 0.0j
		zm_mode = 0.0 + 0.0j
		Zp.append(zp_mode)
		Zm.append(zm_mode)
		for n in range(1, self.nBins / 2):
			zp_mode.real = -Zm[self.nBins / 2 - n].real
			zp_mode.imag =  Zm[self.nBins / 2 - n].imag
			zm_mode.real = -Zp[self.nBins / 2 - n].real
			zp_mode.imag =  Zp[self.nBins / 2 - n].imag
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("X", Zp, Zm)
	if(self.useY != 0):
		Zp = []
		Zm = []
		for n in range(self.nBins / 2):
			Freq_p = Freq0 * (n + qY)
			Freq_m = Freq0 * (n - qY)
			sign_m = 1.0
			if(Freq_m < 0):
				Freq_m = -Freq_m
				sign_m = -1.0
			zp_mode = bilinterp(BetaRel, Freq_p,\
				self.bet_range, self.freq_range,\
				self.bet_tuple, self.freq_tuple,\
				self.zy_bf)
			zm_mode = bilinterp(BetaRel, Freq_m,\
				self.bet_range, self.freq_range,\
				self.bet_tuple, self.freq_tuple,\
				self.zy_bf)
			zm_mode.real = sign_m * zm_mode.real
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		zp_mode = 0.0 + 0.0j
		zm_mode = 0.0 + 0.0j
		Zp.append(zp_mode)
		Zm.append(zm_mode)
		for n in range(1, self.nBins / 2):
			zp_mode.real = -Zm[self.nBins / 2 - n].real
			zp_mode.imag =  Zm[self.nBins / 2 - n].imag
			zm_mode.real = -Zp[self.nBins / 2 - n].real
			zp_mode.imag =  Zp[self.nBins / 2 - n].imag
			Zp.append(zp_mode)
			Zm.append(zm_mode)
		self.timpedance.assignImpedance("Y", Zp, Zm)

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
        y = y_tuple[0]
        return y
    if x > x_tuple[n_tuple]:
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
        for ny in range(ny_tuple + 1):
            vf = fxy[0][ny]
            f_tuple.append(vf)
    elif x > x_tuple[nx_tuple]:
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

