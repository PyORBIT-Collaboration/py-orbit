"""
Module. Includes classes for 1D longidutinal space charge accelerator nodes.
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

#import longitudinal space charge package
from spacecharge import LSpaceChargeCalc

#-----------------------------------------------------------------------------
# Node for impedance as function of node number
#-----------------------------------------------------------------------------

class SC1D_AccNode(DriftTEAPOT):

    def __init__(self, b_a, phaseLength, nMacrosMin, useSpaceCharge,\
                 nBins, name = "long sc node"):
        """
            Constructor. Creates the SC1D-teapot element.
        """
        DriftTEAPOT.__init__(self, name)
        self.lspacecharge = LSpaceChargeCalc(b_a, phaseLength, nMacrosMin,\
                                             useSpaceCharge, nBins)
        self.setType("long sc node")
        self.setLength(0.0)

    def trackBunch(self, bunch):
        """
            The SC1D-teapot class implementation of the
            AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        self.lspacecharge.trackBunch(bunch)	#track method goes here

    def track(self, paramsDict):
        """ 
            The SC1D-teapot class implementation of the
            AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.lspacecharge.trackBunch(bunch)	#track method goes here

    def assignImpedance(self, py_cmplx_arr):
        self.lspacecharge.assignImpedance(py_cmplx_arr)

#-----------------------------------------------------------------------------
# Node for impedance as function of frequency
#-----------------------------------------------------------------------------

class FreqDep_SC1D_AccNode(DriftTEAPOT):

    def __init__(self, b_a, phaseLength, nMacrosMin, useSpaceCharge,\
                 nBins, bunch, impeDict,\
		 name = "freq. dep. long sc node"):
        """
            Constructor. Creates the FreqDep_SC1D-teapot element.
        """
        DriftTEAPOT.__init__(self, name)
        self.lspacecharge = LSpaceChargeCalc(b_a, phaseLength, nMacrosMin,\
                                             useSpaceCharge, nBins)
        self.setType("freq. dep. long sc node")
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
        self.lspacecharge.assignImpedance(Z)

    def trackBunch(self, bunch):
        """
            The FreqDep_SC1D-teapot class implementation of
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
        self.lspacecharge.assignImpedance(Z)
        self.lspacecharge.trackBunch(bunch)

    def track(self, paramsDict):
        """
            The FreqDep_SC1D-teapot class implementation of
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
        self.lspacecharge.assignImpedance(Z)
        self.lspacecharge.trackBunch(bunch)

#-----------------------------------------------------------------------------
# Node for impedance as function of beta and frequency
#-----------------------------------------------------------------------------

class BetFreqDep_SC1D_AccNode(DriftTEAPOT):

    def __init__(self, b_a, phaseLength, nMacrosMin, useSpaceCharge,\
                 nBins, bunch, impeDict,\
		 name = "freq. dep. long sc node"):
        """
            Constructor. Creates the BetFreqDep_SC1D-teapot element.
        """
        DriftTEAPOT.__init__(self, name)
        self.lspacecharge = LSpaceChargeCalc(b_a, phaseLength, nMacrosMin,\
                                             useSpaceCharge, nBins)
        self.setType("beta-freq. dep. long sc node")
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
            z_mode = bilinterp(BetaRel, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.z_bf)
            Z.append(z_mode)
        self.lspacecharge.assignImpedance(Z)

    def trackBunch(self, bunch):
        """
            The BetFreqDep_SC1D-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = bilinterp(BetaRel, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.z_bf)
            Z.append(z_mode)
        self.lspacecharge.assignImpedance(Z)
        self.lspacecharge.trackBunch(bunch)

    def track(self, paramsDict):
        """
            The BetFreqDep_SC1D-teapot class implementation of
            the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phaseLength
        Z = []
        for n in range(self.nBins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = bilinterp(BetaRel, freq_mode, self.bet_range,\
                               self.freq_range, self.bet_tuple,\
                               self.freq_tuple, self.z_bf)
            Z.append(z_mode)
        self.lspacecharge.assignImpedance(Z)
        self.lspacecharge.trackBunch(bunch)


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

