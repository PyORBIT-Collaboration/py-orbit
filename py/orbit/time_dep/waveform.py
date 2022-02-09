import sys
import os
import math

from orbit.utils import NamedObject, TypedObject

class Waveform(NamedObject, TypedObject):
	"""
	The base abstract class of waveforms hierarchy.
	"""
	def __init__(self, name = "no name"):
		NamedObject.__init__(self, name)
		TypedObject.__init__(self, "base waveform")

class KickerWaveform(Waveform):
	"""
	Subclass of Waveform class. The abstract class of kicker waveforms.
	"""
	def __init__(self, name = "no name"):
		Waveform.__init__(self, name)
		self.setType("kicker waveform")
		self.kx = None
		self.ky = None

	def setKx(self,kx):
		self.kx = kx

	def getKx(self):
		return self.kx

	def setKy(self,ky):
		self.ky = ky

	def getKy(self):
		return self.ky

class ConstantKickerWaveform(KickerWaveform):
	"""
	Kicker waveform of constant strength.
	"""
	def __init__(self, name = "no name"):
		KickerWaveform.__init__(self, name)

	def initialize(self, kx, ky):
		self.setKx(kx)
		self.setKy(ky)

	def calc(self, time):
		pass

class SquareRootKickerWaveform(KickerWaveform):
	"""
	Square Root Waveform of Kicker.
	"""
	def __init__(self, name = "no name"):
		KickerWaveform.__init__(self, name)
		self.__initial = []

	def initialize(self, t1, t2, kx, ky, dx, dy):
		self.__initial = (t1, t2, dx, dy)
		self.setKx(kx)
		self.setKy(ky)

	def calc(self, time):
		t = time
		(t1, t2, dx, dy) = self.__initial
		dt = sqrt((t - t1) / (t2 - t1))
		self.setKx(self.getKx() * (1 - dt) + dx)
		self.setKy(self.getKy() * (1 - dt) + dy)

class MagnetWaveform(Waveform):
	"""
	Subclass of Waveform class. The abstract class of magnet waveforms.
	"""
	def __init__(self, name = "no name"):
		Waveform.__init__(self, name)
		self.setType("magnet waveform")
		self.strength = None

	def setStrength(self, strength):
		self.strength = strength

	def getStrength(self):
		return self.strength

class ConstantMagnetWaveform(MagnetWaveform):
	"""
	Magnet waveform of constant strength.
	"""
	def __init__(self, name = "no name"):
		MagnetWaveform.__init__(self, name)

	def initialize(self, strength):
		self.setStrength(strength)

	def calc(self, time):
		pass

class LinearMagnetWaveform(MagnetWaveform):
	"""
	Linear lattice strength variation between t1 and t2
	"""
	def __init__(self, name = "no name"):
		MagnetWaveform.__init__(self, name)
		self.__initial = []

	def initialize(self, t1, t2, s1, s2):
		self.__initial = (t1,t2,s1,s2)

	def calc(self, t):
		(t1, t2, s1, s2) = self.__initial
		if   t < t1: self.strength = s1
		elif t > t2: self.strength = s2
		else:
		  dt = (t - t1) / (t2 - t1)
		  self.setStrength(s1 + dt * (s2 - s1))

class ConstantWaveform:
	"""
	Waveform of constant strength.
	"""
	def __init__(self, syncpart, lattlength, strength):
		self.name = "constant waveform"
		self.syncpart = syncpart
                self.lattlength = lattlength
                self.strength = strength

	def getStrength(self):
		return self.strength


class SquareRootWaveform:
	"""
	Square root waveform.
	"""
	def __init__(self, syncpart, lattlength, ti, tf, si, sf):
		self.name = "square root waveform"
		self.syncpart = syncpart
		self.lattlength = lattlength
                self.ti = ti
                self.tf = tf
                self.si = si
                self.sf = sf

	def getStrength(self):
		time = self.syncpart.time()
                if(time < self.ti):
                        strength = self.si
                elif(time > self.tf):
                        strength = self.sf
                else:
                        dt = math.sqrt((time - self.ti) / (self.tf - self.ti))
                        strength = ((self.si - self.sf) * (1.0 - dt) \
                                         + self.sf)
		return strength


class LinearWaveform:
	"""
	Linear strength variation between ti and tf
	"""
	def __init__(self, syncpart, lattlength, ti, tf, si, sf):
		self.name = "linear waveform"
		self.syncpart = syncpart
		self.lattlength = lattlength
                self.ti = ti
                self.tf = tf
                self.si = si
                self.sf = sf

	def getStrength(self):
		time = self.syncpart.time()
                if(time < self.ti):
                        strength = self.si
                elif(time > self.tf):
                        strength = self.sf
		else:
                        dt = (time - self.ti) / (self.tf - self.ti)
                        strength = (self.sf - self.si) * dt + self.si
		return strength


class OneTimeKick:
	"""
	Deliver a single kick on a specific turn
	"""
	def __init__(self, syncpart, lattlength, kickturn, strength):
		self.name = "linear waveform"
		self.syncpart   = syncpart
		self.lattlength = lattlength
                self.kickturn   = kickturn
                self.strength   = strength
                self.turn = 0

        def setturn(self, turn):
                self.turn = turn

	def getStrength(self):
                strength = 0.0
                if(self.turn == self.kickturn): strength = self.strength
		return strength


class JPARC_08:
	"""
	Generates noise signal QFL
	"""
	def __init__(self, syncpart, lattlength):
		self.name = "JPARC_08"
		self.syncpart   = syncpart
		self.lattlength = lattlength

	def getStrength(self):
                # time in msec, since this taken from old ORBIT
                time = 1000.0 * self.syncpart.time()
                QFLBase = 0.156944849901
                QFLNoise = time * (-0.0002569    + \
                        time * (0.00054243       + \
                        time * (-0.00012826      + \
                        time * (0.000016387      + \
                        time * (-0.0000012873    + \
                        time * (0.000000061252   + \
                        time * (-0.0000000016099 + \
                        time * (0.000000000017885))))))))
                strength = (QFLBase + QFLNoise) / QFLBase
		return strength


class JPARC_09:
	"""
	Generates noise signal QDL
	"""
	def __init__(self, syncpart, lattlength):
		self.name = "JPARC_09"
		self.syncpart   = syncpart
		self.lattlength = lattlength

	def getStrength(self):
                # time in msec, since this taken from old ORBIT
                time = 1000.0 * self.syncpart.time()
                QDLBase = -0.167375195725
                QDLNoise = time * (0.000031569   + \
                        time * (-0.00022683      + \
                        time * (0.000075352      + \
                        time * (-0.000012569     + \
                        time * (0.0000011936     + \
                        time * (-0.000000064555  + \
                        time * (0.0000000018493  + \
                        time * (-0.00000000002179))))))))
                strength = (QDLBase + QDLNoise) / QDLBase
		return strength


class JPARC_15:
	"""
	Generates noise signal QFM
	"""
	def __init__(self, syncpart, lattlength):
		self.name = "JPARC_15"
		self.syncpart   = syncpart
		self.lattlength = lattlength

	def getStrength(self):
                # time in msec, since this taken from old ORBIT
                time = 1000.0 * self.syncpart.time()
                QFMBase = 0.150484620592
                QFMNoise = time * (-0.00039206   + \
                        time * (0.00072671       + \
                        time * (-0.00015812      + \
                        time * (0.000018344      + \
                        time * (-0.0000013108    + \
                        time * (0.000000057474   + \
                        time * (-0.0000000014143 + \
                        time * (0.000000000014925))))))))
                strength = (QFMBase + QFMNoise) / QFMBase
		return strength


class JPARC_8to9and15:
	"""
	Generates noise signal
	"""
	def __init__(self, syncpart, lattlength, QBase, Coeffs):
		self.name = "JPARC_8to9and15"
		self.syncpart   = syncpart
		self.lattlength = lattlength
                self.QBase      = QBase
                self.Coeffs     = Coeffs

	def getStrength(self):
                # time in msec, since this taken from old ORBIT
                time = 1000.0 * self.syncpart.time()
                QNoise = time * (self.Coeffs[0] + \
                        time * (self.Coeffs[1] + \
                        time * (self.Coeffs[2] + \
                        time * (self.Coeffs[3] + \
                        time * (self.Coeffs[4] + \
                        time * (self.Coeffs[5] + \
                        time * (self.Coeffs[6] + \
                        time * (self.Coeffs[7]))))))))
                strength = (self.QBase + QNoise) / self.QBase 
		return strength


class JPARC_12:
	"""
	Pranab 2011/07/27: copied from rcsdel05.
        For time-dependent Sextupoles
	"""
	def __init__(self, syncpart, lattlength):
		self.name = "JPARC_12"
		self.syncpart = syncpart
		self.lattlength = lattlength

	def getStrength(self):
                Bsync0 = 0.1808950855 # JPARC 181 MeV
                gamma0 = (self.syncpart.mass + 0.181) / self.syncpart.mass
                mom0   = math.sqrt((gamma0 + 1.0) * (gamma0 - 1.0))
                mom    = math.sqrt((self.syncpart.gamma + 1.0) * \
                                   (self.syncpart.gamma - 1.0))
                strength = mom0 / mom
		return strength


class JPARC_14:
	"""
	Pranab 2011/07/27: copied from rcsdel05.
        For time-dependent Sextupoles
	"""
	def __init__(self, syncpart, lattlength):
		self.name = "JPARC_14"
		self.syncpart = syncpart
		self.lattlength = lattlength

	def getStrength(self):
                Bsync0 = 0.2828661203 # JPARC 400 MeV
                gamma0 = (self.syncpart.mass + 0.400) / self.syncpart.mass
                mom0   = math.sqrt((gamma0 + 1.0) * (gamma0 - 1.0))
                mom    = math.sqrt((self.syncpart.gamma + 1.0) * \
                                   (self.syncpart.gamma - 1.0))
                strength = mom0 / mom
		return strength


class JPARC_12and14:
	"""
	Pranab 2011/07/27: copied from rcsdel05.
        For time-dependent Sextupoles
	"""
	def __init__(self, syncpart, lattlength, ke0):
		self.name = "JPARC_12and14"
		self.syncpart = syncpart
		self.lattlength = lattlength
                self.ke0 = ke0

	def getStrength(self):
                gamma0 = (self.syncpart.mass + ke0) / self.syncpart.mass
                mom0   = math.sqrt((gamma0 + 1.0) * (gamma0 - 1.0))
                mom    = math.sqrt((self.syncpart.gamma + 1.0) * \
                                   (self.syncpart.gamma - 1.0))
                strength = mom0 / mom
		return strength


class JPARC_16to19:
	"""
        Pranab: BM 1kZ ripple WF16-19
        tinit time is adjusted to match the measured beam
        oscillation(turn-by-turn)
	"""
	def __init__(self, syncpart, lattlength, kebot, \
                     tinit, tdt, BMripple):
		self.name = "JPARC_16to19"
		self.syncpart = syncpart
		self.lattlength = lattlength
                self.kebot = kebot
                self.tinit = tinit
                self.tdt = tdt
                self.BMripple = BMripple

	def getStrength(self):
                time = self.syncpart.time()
                BcalNew1 = math.sin(math.pi * \
                                    (time - self.tinit) / self.tdt)
                strength = BcalNew1 * self.kebot * self.BMripple / \
                           self.syncpart.kinEnergy()
		return strength


class JPARC_20:
	"""
	Pranab: DC leakage field with strength factor
	"""
	def __init__(self, syncpart, lattlength, ke0, LeakFld):
		self.name = "JPARC_20"
		self.syncpart = syncpart
		self.lattlength = lattlength
                self.ke0 = ke0
                self.LeakFld = LeakFld

	def getStrength(self):
                gamma0 = (self.syncpart.mass() + self.ke0) / \
                         self.syncpart.mass()
                mom0   = math.sqrt((gamma0 + 1.0) * (gamma0 - 1.0))
                mom    = math.sqrt((self.syncpart.gamma() + 1.0) * \
                                   (self.syncpart.gamma() - 1.0))
                strength = mom0 * self.LeakFld / mom
		return strength


class JPARC_21to27:
	"""
	Pranab: DC leakage field with strength factor
	"""
	def __init__(self, syncpart, lattlength, QNoise):
		self.name = "JPARC_21to27"
		self.syncpart = syncpart
		self.lattlength = lattlength
                self.QNoise = QNoise
                self.time = []
                self.dQ = []
                Noise_in = open(self.QNoise, "r")
                for line in Noise_in.readlines():
                        splitline = line.split()
                        value = map(float, splitline)
                        self.time.append(value[0])
                        self.dQ.append(value[1])
                self.turn = 0

	def getStrength(self):
                dQVal = self.dQ[self.turn]
                self.turn = self.turn + 1
                strength = 1.0 + dQVal
		return strength


class JPARC_28:
	"""
	Pranab 2016/03/14: Sextupole user-defined pattern.
        Linear between ti -> tf, constant outside.
        Set to J-PARC parameters.
	"""
	def __init__(self, syncpart, lattlength, ti, tf, si, sf):
		self.name = "JPARC_28"
		self.syncpart = syncpart
		self.lattlength = lattlength
                self.ti = ti
                self.tf = tf
                self.si = si
                self.sf = sf

	def getStrength(self):
                Bsync0 = 0.2828661203 # JPARC 400 MeV
                gamma0 = (self.syncpart.mass() + 0.400) / \
                         self.syncpart.mass()
                mom0   = math.sqrt((gamma0 + 1.0) * (gamma0 - 1.0))
                mom    = math.sqrt((self.syncpart.gamma() + 1.0) * \
                                   (self.syncpart.gamma() - 1.0))
                Bsync  = Bsync0 * mom / mom0
                stri   = Bsync0 * self.si / Bsync
                strf   = self.sf
		time = self.syncpart.time()
                if(time <= self.ti):
                        strength = stri
                elif(time > self.tf):
                        strength = strf
		else:
                        dt = (time - self.ti) / \
                             (self.tf - self.ti)
                        strength = (strf - stri) * dt + stri
		return strength


class JPARC_31to34:
	"""
	Pranab 2016/04/12: Square Root Waveform for multipole hkickers.
	"""
	def __init__(self, syncpart, lattlength, ti, tf, tflin, tsclin, \
                     tcutoff, six, sfx, siy, sfy):
		self.name = "JPARC_31to34"
		self.syncpart = syncpart
		self.lattlength = lattlength
                self.ti = ti
                self.tf = tf
                self.tflin = tflin
                self.tsclin = tsclin
                self.tcutoff = tcutoff
                self.six = six
                self.sfx = sfx
                self.siy = siy
                self.sfy = sfy

	def getStrength(self):
		time = self.syncpart.time()
                dt = 0.0
                if(time <= self.tflin):
                        dt = (time - self.ti) / (self.tsclin - self.ti)
                        if(dt < 0.0): dt = 0.0
                        if(dt > 1.0): dt = 1.0
                        sx = (self.sfx - self.six) * dt + self.six
                        sy = (self.sfy - self.siy) * dt + self.siy
                elif((time > self.tflin) and (time  <= self.tf)):
                        dt = math.sqrt((time - self.ti) / \
                                       (self.tf - self.ti))
                        if(dt < 0.0): dt = 0.0
                        if(dt > 1.0): dt = 1.0
                        sx = (self.sfx - self.six) * dt + self.six
                        sy = (self.sfy - self.siy) * dt + self.siy
		elif((time > self.tf) and (time  <= self.tcutoff)):
                        dt = (time - self.tf) / \
                             (self.tcutoff - self.tf)
                        if(dt < 0.0): dt = 0.0
                        if(dt > 1.0): dt = 1.0
                        sx = self.sfx * (1.0 - dt)
                        sy = self.sfy * (1.0 - dt)
                else:
                        sx = 0.0
                        sy = 0.0
                strength = sx
		return strength
