#!/usr/bin/env python

"""
This is not a parallel version! 
"""


import math
import random
import sys

class UniformLongDist:
	""" 
	This class generates uniform longitudinal distribution coordinates
	"""
	def __init__(self, zmin, zmax, sp, eoffset, deltaEfrac):
		self.name = "JohoLongitudinal"
		self.zmin = zmin
		self.zmax = zmax
		self.sp = sp
		self.ekinetic = sp.kinEnergy()
		self.eoffset = eoffset
		self.deltaEfrac = deltaEfrac

	def getCoordinates(self):
		
		zinj = (self.zmin + (self.zmax - self.zmin) * random.random())
		dEinj = self.eoffset + (self.ekinetic * -self.deltaEfrac * (1 - 2*random.random()))

		return (zinj,dEinj)


class UniformLongDistPaint:
    '''
        This class generates time dependent uniform longitudinal distribution coordinates
        according to user defined (mathematical) functions for zmin and zmax
        '''
    def __init__(self, zminFunc, zmaxFunc, sp, eoffset, deltaEfrac):
        
        
        #checks for correct input arguement types (to a certain extent)
        if not (type(zminFunc) is list and type(zmaxFunc) is list):
            sys.exit("ERROR from class 'UniformLongDistPaint': input arguements zminFunc and zmaxFunc must be a list of pairs")
        else:
            if len(zminFunc) < 2 or len(zmaxFunc) < 2:
                sys.exit("ERROR from class 'UniformLongDistPaint': input arguements zminFunc and zmaxFunc must have at least two elements")
            else:
                if not (type(zminFunc[0]) is list and type(zmaxFunc[0]) is list):
                    sys.exit("ERROR from class 'UniformLongDistPaint': input arguements zminFunc and zmaxFunc must be a list of pairs")
            
            
        #sorts functions for quick interpolation calculations
        zminFunc = sorted(zminFunc, key= lambda x: x[0])
        zmaxFunc = sorted(zmaxFunc, key= lambda x: x[0])
            
        self.name = "JohoLongitudinalPaint"
        self.zminFunc = zminFunc
        self.zmaxFunc = zmaxFunc
        self.sp = sp
        self.ekinetic = sp.kinEnergy()
        self.eoffset = eoffset
        self.deltaEfrac = deltaEfrac
		# These help vary the number of macro particles according to the PW
        self.frac_change = 1.0
        self.last_length = -1
        
        
    def getCoordinates(self):
        #time in seconds
		zminNow = interpolate(self.zminFunc,self.sp.time())
		zmaxNow = interpolate(self.zmaxFunc,self.sp.time())
		if zminNow >= zmaxNow: print "Warning from getCoordinates call: zmin >= zmax"
		
		length = zmaxNow-zminNow
		if self.last_length != -1: self.frac_change = length/self.last_length
		self.last_length = length
		
		zinj = (zminNow + (zmaxNow - zminNow) * random.random())
		dEinj = self.eoffset + (self.ekinetic * -self.deltaEfrac * (1 - 2*random.random()))

		return (zinj,dEinj)




class GULongDist():
	"""
	This class generates random intial longitudinal coordinates for a
	distribution uniform in phi and gaussian in dE.
	"""

	def __init__(self, zmin, zmax, sp, emean, esigma, etrunc, emin, emax):
        
		self.name = "JohoLongitudinal"
		self.zmin = zmin
		self.zmax = zmax
		self.sp = sp
		self.ekinetic = sp.kinEnergy()
		self.emean = emean
		self.esigma = esigma
		self.etrunc = etrunc
		self.emin = emin
		self.emax = emax
	

	def getCoordinates(self):
	
		zinj = self.zmin + (self.zmax - self.zmin) * random.random()
		
		tol = 1e-6
		ymin = 0.
		ymax = 10.
		pmin = 0.
		pmax = 1.
		

		if(self.etrunc != 0):
			if(self.emin >= self.emean):
				pmin = 0.5 + 0.5 * erf( (self.emin - self.emean) / (math.sqrt(2.) * self.esigma) )
			else:
				pmin = 0.5 - 0.5 * erf( (self.emean - self.emin) / (math.sqrt(2.) * self.esigma) )

			if(self.emax >= self.emean):
				pmax = 0.5 + 0.5 * erf( (self.emax - self.emean) / (math.sqrt(2.) * self.esigma) )
			else:
				pmax = 0.5 - 0.5 * erf( (self.emean - self.emax) / (math.sqrt(2.) * self.esigma) )

		prand = pmin + (pmax - pmin) * random.random()
		
		while( (erf(ymax) - math.fabs(2.*prand - 1)) < 0. ):
			ymax *= 10.

		root = rootNorm(ymin, ymax, prand, tol)

		if (prand >= 0.5 ):
			einj = self.emean + (math.sqrt(2.) * self.esigma * root)
		else:
			einj = self.emean - (math.sqrt(2.) * self.esigma * root)

		dEinj = einj - self.ekinetic
		
		return (zinj, dEinj)

			

class SNSESpreadDist():
	"""
	Class for generating random initial particle coordinat4es for a uniform longitudinal
	distribution and a gaussian energy distribution, and then adding sinusoidal energy spread
	and random centroid jitter
	"""
	def __init__(self, lattlength, zmin, zmax, tailfraction, sp, emean, esigma, etrunc, emin, emax, ecparams, esparams):
		self.name = "JohoLongitudinal"
		self.lattlength = lattlength
		self.zmin = zmin
		self.zmax = zmax
		self.tailfraction = tailfraction
		self.sp = sp
		self.emean = emean
		self.ekinetic = sp.kinEnergy()
		self.esigma = esigma
		self.etrunc = etrunc
		self.emin = emin
		self.emax = emax
		self.ecparams = ecparams
		self.esparams = esparams
	
	
	def getCoordinates(self):
	
		if(random.random() > (self.tailfraction * self.lattlength / (self.lattlength - self.zmax + self.zmin))):
			#Put it in the main distribution
			zinj = self.zmin + (self.zmax - self.zmin) * random.random() 
		
		else:
			#Put it in an extended tail
			zinj = -self.lattlength/2.0 + self.lattlength * random.random()

		tol = 1e-6
		ymin = 0.
		ymax = 10.
		pmin = 0.
		pmax = 1.

		if(self.etrunc != 0):
			if(self.emin >= self.emean):
				pmin = 0.5 + 0.5 * erf( (self.emin - self.emean) / (math.sqrt(2.) * self.esigma) )
			else:
				pmin = 0.5 - 0.5 * erf( (self.emean - self.emin) / (math.sqrt(2.) * self.esigma) )
				
			if(self.emax >= self.emean):
				pmax = 0.5 + 0.5 * erf( (self.emax - self.emean) / (math.sqrt(2.) * self.esigma) )
			else:
				pmax = 0.5 - 0.5 * erf( (self.emean - self.emax) / (math.sqrt(2.) * self.esigma) )
			
		prand = pmin + (pmax - pmin) * random.random()
			
		while( (erf(ymax) - math.fabs(2.*prand - 1)) < 0. ):
			ymax *= 10.
			
		root = rootNorm(ymin, ymax, prand, tol)
			
		if (prand >= 0.5 ):
			einj = self.emean + (math.sqrt(2.) * self.esigma * root)
		else:
			einj = self.emean - (math.sqrt(2.) * self.esigma * root)	

		#If ecmin = 0  then do regular Gaussian distribution.
		#If ecmin >= 0, then this will indicate that a truncated Gaussian distribution is desired.

		(ecmean, ecsigma, ectrunc, ecmin, ecmax, ecdrifti, ecdriftf, drifttime) = self.ecparams
		
		
		pmin = 0.
		pmax = 1.
		ec = 0.
		
		ecdrift = ecdrifti + (ecdriftf - ecdrifti) * self.sp.time() / drifttime
				
		if(ectrunc != 0):
			if(ecmin >= ecmean):
				pmin = 0.5 + 0.5 * erf( (ecmin - ecmean)/(math.sqrt(2.)*ecsigma) )
			else:
				pmin = 0.5 - 0.5 * erf( (ecmean - ecmin)/(math.sqrt(2.)*ecsigma) )
			if(ecmax >= ecmean):
				pmax = 0.5 + 0.5 * erf( (ecmax - ecmean)/(math.sqrt(2.)*ecsigma) )
			else:
				pmax = 0.5 + 0.5 * erf( (ecmean - ecmax)/(math.sqrt(2.)*ecsigma) )

		prand = pmin + (pmax - pmin) * random.random()
	
		while( (erf(ymax) - math.fabs(2.*prand - 1)) < 0. ):
			ymax *= 10.
		
		root = rootNorm(ymin, ymax, prand, tol)
		
		if(prand >= 0.5):
			ec = ecmean + ecdrift + math.sqrt(2.)*ecsigma * root
		else:
			ec = ecmean + ecdrift - math.sqrt(2.)*ecsigma * root
		
		(esnu, esphase, esmax, nulltime) = self.esparams		
				
		iphase = esnu * self.sp.time()
		phasec = 2. * math.pi * (esnu * self.sp.time() - iphase)
		turntime = self.lattlength / (self.sp.beta() * 2.998e8)
		phasephifac = esnu * turntime
		phasephi = phasephifac * zinj*math.pi/(self.lattlength/2.)
		phase = phasec + phasephi + esphase
		es = esmax * math.sin(phase)
		tfac = 0.

		if(nulltime > 0.):
			if(self.sp.time() > nulltime):
				tfac = 0.
			else:
				tfac = math.pow( 1. - self.sp.time() / nulltime, 0.5)

		es *= tfac;
	
		dEinj = einj + ec + es - self.ekinetic

		return(zinj, dEinj)



class SNSESpreadDistPaint():
	"""
		This class generates time dependent SNSESpreadDistPaint distribution coordinates
        according to user defined (mathematical) functions for zmin and zmax
		"""
	def __init__(self, lattlength, zminFunc, zmaxFunc, tailfraction, sp, emean, esigma, etrunc, emin, emax, ecparams, esparams):
		
		#checks for correct input arguement types (to a certain extent)
		if not (type(zminFunc) is list and type(zmaxFunc) is list):
			sys.exit("ERROR from class 'UniformLongDistPaint': input arguements zminFunc and zmaxFunc must be a list of pairs")
		else:
			if len(zminFunc) < 2 or len(zmaxFunc) < 2:
				sys.exit("ERROR from class 'UniformLongDistPaint': input arguements zminFunc and zmaxFunc must have at least two elements")
			else:
				if not (type(zminFunc[0]) is list and type(zmaxFunc[0]) is list):
					sys.exit("ERROR from class 'UniformLongDistPaint': input arguements zminFunc and zmaxFunc must be a list of pairs")
		
		
        #sorts functions for quick interpolation calculations
		zminFunc = sorted(zminFunc, key= lambda x: x[0])
		zmaxFunc = sorted(zmaxFunc, key= lambda x: x[0])

		self.name = "JohoLongitudinalPaint"
		self.lattlength = lattlength
		self.zminFunc = zminFunc
		self.zmaxFunc = zmaxFunc
		self.tailfraction = tailfraction
		self.sp = sp
		self.emean = emean
		self.ekinetic = sp.kinEnergy()
		self.esigma = esigma
		self.etrunc = etrunc
		self.emin = emin
		self.emax = emax
		self.ecparams = ecparams
		self.esparams = esparams
		#These help vary the number of macro partices according to PW
		self.frac_change = 1.0
		self.last_length = -1




	def getCoordinates(self):
		
		zminNow = interpolate(self.zminFunc,self.sp.time())
		zmaxNow = interpolate(self.zmaxFunc,self.sp.time())
		
		length = zmaxNow-zminNow
		if self.last_length != -1: self.frac_change = length/self.last_length
		self.last_length = length

		
		if(random.random() > (self.tailfraction * self.lattlength / (self.lattlength - zmaxNow + zminNow))):
			#Put it in the main distribution
			zinj = zminNow + (zmaxNow - zminNow) * random.random()
		
		else:
			#Put it in an extended tail
			zinj = -self.lattlength/2.0 + self.lattlength * random.random()
		
		tol = 1e-6
		ymin = 0.
		ymax = 10.
		pmin = 0.
		pmax = 1.
		
		if(self.etrunc != 0):
			if(self.emin >= self.emean):
				pmin = 0.5 + 0.5 * erf( (self.emin - self.emean) / (math.sqrt(2.) * self.esigma) )
			else:
				pmin = 0.5 - 0.5 * erf( (self.emean - self.emin) / (math.sqrt(2.) * self.esigma) )
			
			if(self.emax >= self.emean):
				pmax = 0.5 + 0.5 * erf( (self.emax - self.emean) / (math.sqrt(2.) * self.esigma) )
			else:
				pmax = 0.5 - 0.5 * erf( (self.emean - self.emax) / (math.sqrt(2.) * self.esigma) )
		
		prand = pmin + (pmax - pmin) * random.random()
		
		while( (erf(ymax) - math.fabs(2.*prand - 1)) < 0. ):
			ymax *= 10.
		
		root = rootNorm(ymin, ymax, prand, tol)
		
		if (prand >= 0.5 ):
			einj = self.emean + (math.sqrt(2.) * self.esigma * root)
		else:
			einj = self.emean - (math.sqrt(2.) * self.esigma * root)
		
		#If ecmin = 0  then do regular Gaussian distribution.
		#If ecmin >= 0, then this will indicate that a truncated Gaussian distribution is desired.
		
		(ecmean, ecsigma, ectrunc, ecmin, ecmax, ecdrifti, ecdriftf, drifttime) = self.ecparams
		
		
		pmin = 0.
		pmax = 1.
		ec = 0.
		
		ecdrift = ecdrifti + (ecdriftf - ecdrifti) * self.sp.time() / drifttime
		
		if(ectrunc != 0):
			if(ecmin >= ecmean):
				pmin = 0.5 + 0.5 * erf( (ecmin - ecmean)/(math.sqrt(2.)*ecsigma) )
			else:
				pmin = 0.5 - 0.5 * erf( (ecmean - ecmin)/(math.sqrt(2.)*ecsigma) )
			if(ecmax >= ecmean):
				pmax = 0.5 + 0.5 * erf( (ecmax - ecmean)/(math.sqrt(2.)*ecsigma) )
			else:
				pmax = 0.5 + 0.5 * erf( (ecmean - ecmax)/(math.sqrt(2.)*ecsigma) )
		
		prand = pmin + (pmax - pmin) * random.random()
		
		while( (erf(ymax) - math.fabs(2.*prand - 1)) < 0. ):
			ymax *= 10.
		
		root = rootNorm(ymin, ymax, prand, tol)
		
		if(prand >= 0.5):
			ec = ecmean + ecdrift + math.sqrt(2.)*ecsigma * root
		else:
			ec = ecmean + ecdrift - math.sqrt(2.)*ecsigma * root
		
		(esnu, esphase, esmax, nulltime) = self.esparams
		
		iphase = esnu * self.sp.time()
		phasec = 2. * math.pi * (esnu * self.sp.time() - iphase)
		turntime = self.lattlength / (self.sp.beta() * 2.998e8)
		phasephifac = esnu * turntime
		phasephi = phasephifac * zinj*math.pi/(self.lattlength/2.)
		phase = phasec + phasephi + esphase
		es = esmax * math.sin(phase)
		tfac = 0.
		
		if(nulltime > 0.):
			if(self.sp.time() > nulltime):
				tfac = 0.
			else:
				tfac = math.pow( 1. - self.sp.time() / nulltime, 0.5)
		
		es *= tfac;
		
		dEinj = einj + ec + es - self.ekinetic
		
		return(zinj, dEinj)


def interpolate(List, time):
    
    ListMax = List[len(List)-1][0]
    ListMin = List[0][0]
    
    #boundry cases
    if time >= ListMax:
        p1 = List[len(List)-2]
        p2 = List[len(List)-1]
    elif time <= ListMin:
        p1 = List[1]
        p2 = List[0]
    else:
        #binary search that returns the closest x-coordinate
        e = len(List) -1 #end
        f = 0 # front
        i = int(math.floor((e+f)/2))
        closest = i
        while f<=e:
            if math.fabs(time-List[i][0]) <= math.fabs(time-List[closest][0]):
                closest = i
            if List[i][0] > time:
                e = i -1
                i = int(math.floor((e+f)/2))
            elif List[i][0] < time:
                f = i +1
                i = int(math.floor((e+f)/2))
            else:
                return List[i][1]
        
        if List[i][0] < time:
            p1 = List[i]
            p2 = List[i+1]
        else:
            p1 = List[i-1]
            p2 = List[i]
    
    m = (p2[1]-p1[1])/(p2[0]-p1[0])
    b = p1[1] - m*p1[0]
    
    return m*time + b


		   
def rootNorm(ymin, ymax, prand, tol):
	"""
	Finds the roots of the (Gauss?) function
	"""

	imax = 40
	rtbis = 0.
	dx = 0. 
	xmid = 0.
	fmid = 0.

	if( (erf(ymin) - math.fabs(2*prand -1 )) < 0.):
		rtbis = ymin
		dx = ymax - ymin
	else:
		rtbis = ymax
		dx = ymin - ymax

	for i in xrange(imax):
		dx = dx * 0.5
		xmid = rtbis + dx
		fmid = erf(xmid) - math.fabs(2 * prand - 1)
		if( fmid <= 0.):
			rtbis = xmid
		if( (math.fabs(dx) < tol) or (fmid == 0) ):
			return rtbis
	
	return rtbis

def erf(z):
	t = 1.0 / (1.0 + 0.5 * abs(z))
	# use Horner's method
	ans = 1 - t * math.exp( -z*z -  1.26551223 +
						   t * ( 1.00002368 +
								t * ( 0.37409196 + 
									 t * ( 0.09678418 + 
										  t * (-0.18628806 + 
											   t * ( 0.27886807 + 
													t * (-1.13520398 + 
														 t * ( 1.48851587 + 
															  t * (-0.82215223 + 
																   t * ( 0.17087277))))))))))
	if z >= 0.0:
		return ans
	else:
		return -ans


