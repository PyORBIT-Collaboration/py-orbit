#!/usr/bin/env python

"""
This is not a parallel version! ???
"""


import math
import random
import sys

class TwissContainer:
	""" 
	Keeps the twiss paremeters alpha, beta and the emittance.
	Calculates the normalized value (u**2+(alpha*u + beta*u')**2)/(beta*emittance),
	which is (gamma*u**2+2*alpha*u*u'+beta*u'**2)/(emittance).
	Translates the normalized values u and up to the non-normalized ones.
	"""
	
	def __init__(self, alpha,beta, emittance):
		self.alpha = alpha
		self.beta = beta
		self.gamma = (1.0+self.alpha**2)/self.beta
		self.emittance = emittance
		self.__initialize()
		
	def __initialize(self):
		self.u_max = math.sqrt(self.beta*self.emittance)
		self.up_coeff = math.sqrt(self.emittance/self.beta)
		self.gamma = (1.0+self.alpha**2)/self.beta
		self.up_max = math.sqrt(self.gamma*self.emittance)		
		
	def setEmittance(self, emittance):
		self.emittance = emittance
		self.__initialize()		
		
	def getNormalizedH(self,u,up):
		""" 
		Returns (u**2+(alpha*u + beta*u')**2)/(beta*emittance) 
		"""
		return (u**2+(self.alpha*u + self.beta*up)**2)/(self.beta*self.emittance)
			
	def getU_Max(self):
		"""
		Returns the maximal value of u.
		"""
		return self.u_max
			
	def getUP_Max(self):
		"""
		Returns the maximal value of u'.
		"""
		return self.up_max
		
	def getU_UP(self,u_norm,up_norm):
		"""
		Returns the coordinate and momentum for normilized ones
		u = sqrt(beta*emittance)*u_norm
		up = sqrt(emittance/beta)*(up_norm - alpha*u_norm)
		"""
		u = self.u_max*u_norm
		up = self.up_coeff*(up_norm - self.alpha*u_norm)
		return (u,up)
		
	def getAlphaBetaGammaEmitt(self):
		return (self.alpha,self.beta,self.gamma,self.emittance)
		
	def getAlphaBetaEmitt(self):
		return (self.alpha,self.beta,self.emittance)
		

class KVDist1D:
	"""
	Generates the 1D KV-distribution. The input emittance in the TwissConatainer
	is a rms emittance. The generated distribution will give the same value. Remember
	that 100% emittance is 2 times bigger for 1D KV distribution.
	"""
	def __init__(self, twiss = TwissContainer(0.,1.,1.)):
		""" Constructor """
		(alpha,beta,emittance) = twiss.getAlphaBetaEmitt()
		self.twiss = TwissContainer(alpha,beta,2*emittance)
		self.sign_choices = (-1.,1.)
		
	def getCoordinates(self):
		""" Return (u,up) distributed for the 1D KV-distribution. """
		x_norm = math.sin(2*math.pi*(random.random()-0.5))
		xp_norm = random.choice(self.sign_choices)*math.sqrt(1.0 - x_norm**2)
		return self.twiss.getU_UP(x_norm,xp_norm)
	
	def getTwissContainers(self):
		""" Returns the twiss container. """
		return (self.twiss,)
		
	
class WaterBagDist1D:
	""" 
	Generates the Water Bag 1D distribution.The input emittance in the TwissConatainer
	is a rms emittance. The generated distribution will give the same value. Remember
	that 100% emittance is 4 times bigger for 1D WaterBag distribution. 
	"""
	def __init__(self, twiss = TwissContainer(0.,1.,1.)):
		""" Constructor """
		(alpha,beta,emittance) = twiss.getAlphaBetaEmitt()
		twiss = TwissContainer(alpha,beta,2*emittance)
		self.kv_dist = KVDist1D(twiss)
		
	def getCoordinates(self):
		""" Return (u,up) distributed for the 1D WaterBag-distribution. """
		(u,up) = self.kv_dist.getCoordinates()
		g = math.sqrt(random.random())
		return (g*u,g*up)

	def getTwissContainers(self):
		""" Returns the twiss container. """
		return self.kv_dist.getTwissContainers()


class KVDist2D:
	"""
	Generates the 2D KV-distribution. The input emittance in the TwissConatainer
	is a rms emittance. The generated distribution will give the same value. Remember
	that 100% emittance is 4 times bigger for 2D KV distribution.
	"""
	def __init__(self, twissX = TwissContainer(0.,1.,1.), twissY = TwissContainer(0.,1.,1.) ):
		""" Constructor """
		(alpha_x,beta_x,emittance_x) = twissX.getAlphaBetaEmitt()
		(alpha_y,beta_y,emittance_y) = twissY.getAlphaBetaEmitt()
		self.twissX = TwissContainer(alpha_x,beta_x,4*emittance_x)
		self.twissY = TwissContainer(alpha_y,beta_y,4*emittance_y)
		
	def getCoordinates(self):
		""" Return (x,xp,y,yp) distributed for the 2D KV-distribution. """
		#x-y plane
		phi = 2*math.pi*(random.random()-0.5)
		rho = math.sqrt(random.random())
		x_norm = rho*math.cos(phi)
		y_norm = rho*math.sin(phi)
		#momentum
		p0 = math.sqrt(math.fabs(1. - rho**2))
		phi = 2*math.pi*(random.random()-0.5)
		xp_norm = p0*math.cos(phi)
		yp_norm = p0*math.sin(phi)
		(x,xp) = self.twissX.getU_UP(x_norm,xp_norm)
		(y,yp) = self.twissY.getU_UP(y_norm,yp_norm)
		return (x,xp,y,yp)

	def getTwissContainers(self):
		""" Returns the (twissX,twissY) containers. """
		return (self.twissX,self.twissY)
	
	
class WaterBagDist2D:
	""" 
	Generates the Water Bag 2D distribution. The input emittance in the TwissConatainer
	is a rms emittance. The generated distribution will give the same value. Remember
	that 100% emittance is 6 times bigger for 2D WaterBag distribution. 
	"""
	def __init__(self, twissX = TwissContainer(0.,1.,1.), twissY = TwissContainer(0.,1.,1.) ):
		""" Constructor """
		(alpha_x,beta_x,emittance_x) = twissX.getAlphaBetaEmitt()
		(alpha_y,beta_y,emittance_y) = twissY.getAlphaBetaEmitt()
		twissX = TwissContainer(alpha_x,beta_x,6.0*emittance_x/4.0)
		twissY = TwissContainer(alpha_y,beta_y,6.0*emittance_y/4.0)		
		self.kv_dist = KVDist2D(twissX,twissY)
		
	def getCoordinates(self):
		""" Return (x,xp,y,yp) distributed for the 2D WaterBag-distribution. """
		(x,xp,y,yp) = self.kv_dist.getCoordinates()
		g = math.sqrt(math.sqrt(random.random()))
		return (g*x,g*xp,g*y,g*yp)

	def getTwissContainers(self):
		""" Returns the (twissX,twissY) containers. """
		return self.kv_dist.getTwissContainers()

class KVDist3D:
	"""
	Generates the 3D KV-distribution.The input emittance in the TwissConatainer
	is a rms emittance. The generated distribution will give the same value. Remember
	that 100% emittance is 8 times bigger for 3D KV distribution.
	"""
	def __init__(self, twissX = TwissContainer(0.,1.,1.), twissY = TwissContainer(0.,1.,1.), twissZ = TwissContainer(0.,1.,1.) ):
		""" Constructor """
		(alpha_x,beta_x,emittance_x) = twissX.getAlphaBetaEmitt()
		(alpha_y,beta_y,emittance_y) = twissY.getAlphaBetaEmitt()		
		(alpha_z,beta_z,emittance_z) = twissZ.getAlphaBetaEmitt()		
		self.twissX = TwissContainer(alpha_x,beta_x,6*emittance_x)
		self.twissY = TwissContainer(alpha_y,beta_y,6*emittance_y)
		self.twissZ = TwissContainer(alpha_z,beta_z,6*emittance_z)
		
	def getCoordinates(self):
		""" Return (x,xp,y,yp,z,zp) distributed for the 3D KV-distribution. """
		#x-y-z-zp plane
		n_limit = 1000
		n_count = 0
		pxy2 = 0.
		x_norm = 0.
		y_norm = 0.
		z_norm = 0.
		zp_norm = 0.
		while(1 < 2):
			n_count = n_count + 1
			x_norm = 2*(random.random()-0.5)
			y_norm = 2*(random.random()-0.5)
			z_norm = 2*(random.random()-0.5)
			zp_norm = 2*(random.random()-0.5)
			pxy2 = 1.0 - x_norm**2 - y_norm**2 - z_norm**2 - zp_norm**2
			if(pxy2 > 0.):
				break
			if(n_count > n_limit):
				print "KVDist3D generator has a problem with Python random module!"
				print "Stop."
				sys.exit(1)
		#make xp-yp plane
		pxy = math.sqrt(pxy2)
		phi = 2*math.pi*(random.random()-0.5)
		xp_norm = pxy*math.cos(phi)
		yp_norm = pxy*math.sin(phi)
		(x,xp) = self.twissX.getU_UP(x_norm,xp_norm)
		(y,yp) = self.twissY.getU_UP(y_norm,yp_norm)
		(z,zp) = self.twissZ.getU_UP(z_norm,zp_norm)
		return (x,xp,y,yp,z,zp)

	def getTwissContainers(self):
		""" Returns the (twissX,twissY,wissZ) containers. """
		return (self.twissX,self.twissY,self.twissZ)	
	
class WaterBagDist3D:
	""" 
	Generates the Water Bag 3D distribution. The input emittance in the TwissConatainer
	is a rms emittance. The generated distribution will give the same value. Remember
	that 100% emittance is 8 times bigger for 3D WaterBag distribution. 
	"""
	def __init__(self, twissX = TwissContainer(0.,1.,1.), twissY = TwissContainer(0.,1.,1.), twissZ = TwissContainer(0.,1.,1.) ):
		""" Constructor """
		(alpha_x,beta_x,emittance_x) = twissX.getAlphaBetaEmitt()
		(alpha_y,beta_y,emittance_y) = twissY.getAlphaBetaEmitt()		
		(alpha_z,beta_z,emittance_z) = twissZ.getAlphaBetaEmitt()		
		twissX = TwissContainer(alpha_x,beta_x,8*emittance_x/6.0)
		twissY = TwissContainer(alpha_y,beta_y,8*emittance_y/6.0)
		twissZ = TwissContainer(alpha_z,beta_z,8*emittance_z/6.0)		
		self.kv_dist = KVDist3D(twissX,twissY,twissZ)
		
	def getCoordinates(self):
		""" Return (x,xp,y,yp,z,zp) distributed for the 3D WaterBag-distribution. """
		(x,xp,y,yp,z,zp) = self.kv_dist.getCoordinates()
		g = math.pow(random.random(),1./6.)
		return (g*x,g*xp,g*y,g*yp,g*z,g*zp)

	def getTwissContainers(self):
		""" Returns the (twissX,twissY,wissZ) containers. """
		return self.kv_dist.getTwissContainers()

class GaussDist1D:
	"""
	Generates the 1D Gauss distribution. exp(-x**2/(2*sigma**2)) The cut_off value is x_cutoff/sigma.
	"""
	def __init__(self, twiss = TwissContainer(0.,1.,1.), cut_off = -1.):
		""" Constructor """
		self.twiss = twiss
		self.cut_off = cut_off
		self.cut_off2 = cut_off*cut_off
		
	def getCoordinates(self):
		""" Return (u,up) distributed for the 1D Gauss distribution. """
		x_norm = random.gauss(0.,1.0)
		xp_norm = random.gauss(0.,1.0)
		if(self.cut_off > 0.):
			while((x_norm**2+xp_norm**2) > self.cut_off2):
				x_norm = random.gauss(0.,1.0)
				xp_norm = random.gauss(0.,1.0)
		return self.twiss.getU_UP(x_norm,xp_norm)
		
	def getTwissContainers(self):
		""" Returns the twiss container. """
		return self.twiss	
		
class GaussDist2D:
	"""
	Generates the 2D Gauss distribution. exp(-x**2/(2*sigma**2)) The cut_off value is x_cutoff/sigma.
	"""
	def __init__(self, twissX = TwissContainer(0.,1.,1.), twissY = TwissContainer(0.,1.,1.), cut_off = -1.):
		""" Constructor """
		self.twissX = twissX
		self.twissY = twissY
		self.gaussX =  GaussDist1D(twissX,cut_off)
		self.gaussY =  GaussDist1D(twissY,cut_off)		
		self.cut_off = cut_off
		
	def getCoordinates(self):
		""" Return (u,up) distributed for the 2D Gauss distribution. """
		(x,xp) = self.gaussX.getCoordinates()
		(y,yp) = self.gaussY.getCoordinates()
		return (x,xp,y,yp)
		
	def getTwissContainers(self):
		""" Returns the (twissX,twissY) containers. """
		return (self.twissX,self.twissY)		
		
		
class GaussDist3D:
	"""
	Generates the 3D Gauss distribution. exp(-x**2/(2*sigma**2)) The cut_off value is x_cutoff/sigma.
	"""
	def __init__(self, twissX = TwissContainer(0.,1.,1.), \
							 twissY = TwissContainer(0.,1.,1.), \
							 twissZ = TwissContainer(0.,1.,1.), \
							 cut_off = -1.):
		""" Constructor """
		self.twissX = twissX
		self.twissY = twissY
		self.twissZ = twissZ
		self.gaussX =  GaussDist1D(twissX,cut_off)
		self.gaussY =  GaussDist1D(twissY,cut_off)		
		self.gaussZ =  GaussDist1D(twissZ,cut_off)		
		self.cut_off = cut_off
		
	def getCoordinates(self):
		""" Return (u,up) distributed for the 3D Gauss distribution. """
		(x,xp) = self.gaussX.getCoordinates()
		(y,yp) = self.gaussY.getCoordinates()
		(z,zp) = self.gaussZ.getCoordinates()
		return (x,xp,y,yp,z,zp)
		
	def getTwissContainers(self):
		""" Returns the (twissX,twissY,twissZ) containers. """
		return (self.twissX,self.twissY,self.twissZ)				
	
#--------------------------------------------------
# Auxilary classes 
#--------------------------------------------------

class TwissAnalysis:
	""" 
	Calculates the rms twiss parameters for 1D,2D, and 3D distributions by 
	using the set of (x,xp), (x,xp,y,yp), and (x,xp,y,yp,z,zp) points.
	There is a c++ replacement for this class BunchTwissAnalysis in the orbit/BunchDiagnostics dir.
	"""
	
	def __init__(self, nD):
		self.nD = nD
		self.x2_avg_v = []
		self.xp2_avg_v = []
		self.x_xp_avg_v = []
		self.x_avg_v = []
		self.xp_avg_v = []
		self.xp_max_v = []
		self.x_max_v = []
		self.xp_min_v = []
		self.x_min_v = []		
		for i in range(self.nD):
			self.x2_avg_v.append(0.)
			self.xp2_avg_v.append(0.)
			self.x_xp_avg_v.append(0.)
			self.x_avg_v.append(0.)
			self.xp_avg_v.append(0.)
			self.xp_max_v.append(-1.0e+38)
			self.x_max_v.append(-1.0e+38)
			self.xp_min_v.append(1.0e+38)
			self.x_min_v.append(1.0e+38)
		self.Np = 0
		
	def init(self):
		"""
		Initilizes the analysis.
		"""
		self.x2_avg_v = []
		self.xp2_avg_v = []
		self.x_xp_avg_v = []
		self.x_avg_v = []
		self.xp_avg_v = []
		self.xp_max_v = []
		self.x_max_v = []		
		self.xp_min_v = []
		self.x_min_v = []				
		for i in range(self.nD):
			self.x2_avg_v.append(0.)
			self.xp2_avg_v.append(0.)
			self.x_xp_avg_v.append(0.)
			self.x_avg_v.append(0.)
			self.xp_avg_v.append(0.)
			self.xp_max_v.append(-1.0e+38)
			self.x_max_v.append(-1.0e+38)
			self.xp_min_v.append(1.0e+38)
			self.x_min_v.append(1.0e+38)	
		self.Np = 0

	def account(self, arr_v):
		"""
		Accounts the data. The arr_v should be a list of 2, 4, or 6 size.
		"""
		for i in range(self.nD):
			self.x_avg_v[i] = self.x_avg_v[i] + arr_v[i*2]
			self.xp_avg_v[i] = self.xp_avg_v[i] + arr_v[i*2+1]					
			self.x2_avg_v[i] = self.x2_avg_v[i] + (arr_v[i*2])**2
			self.xp2_avg_v[i] = self.xp2_avg_v[i] + (arr_v[i*2+1])**2			
			self.x_xp_avg_v[i] = self.x_xp_avg_v[i] + arr_v[i*2+1]*arr_v[i*2]
			x = arr_v[i*2]
			xp = arr_v[i*2+1]	
			if(self.x_max_v[i] < x): self.x_max_v[i] = x
			if(self.xp_max_v[i] < xp): self.xp_max_v[i] = xp
			if(self.x_min_v[i] > x): self.x_min_v[i] = x
			if(self.xp_min_v[i] > xp): self.xp_min_v[i] = xp
			
		self.Np += 1
		
	def getTwiss(self,d):
		"""
		Returns the twiss parameters in the array 
		[alpha,beta,gamma, emittance] for the dimension d.
		"""
		if(d <  0 or d >= self.nD):
			print "Dimention n="+str(d)+" does not exist!"
			sys.exit(1)
		if(self.Np == 0): return (0.,0.,0.) 
		x_avg =  self.x_avg_v[d]
		xp_avg =  self.xp_avg_v[d]		
		x2_avg =  self.x2_avg_v[d]
		xp2_avg =  self.xp2_avg_v[d]
		x_xp_avg =  self.x_xp_avg_v[d]
		x_avg = x_avg/self.Np
		xp_avg = xp_avg/self.Np		
		x2_avg = x2_avg/self.Np - x_avg*x_avg
		x_xp_avg = x_xp_avg/self.Np - x_avg*xp_avg
		xp2_avg = xp2_avg/self.Np - xp_avg*xp_avg

		emitt_rms = math.sqrt(x2_avg*xp2_avg - x_xp_avg*x_xp_avg)
		beta = x2_avg/emitt_rms
		alpha = - x_xp_avg/emitt_rms
		gamma = xp2_avg/emitt_rms		
		return (alpha,beta,gamma,emitt_rms)
		
	def getAvgU_UP(self,d):
		"""
		Returns the (u_avg,up_avg) parameters in the array 
		[u,up] for the dimension d.
		"""
		if(d <  0 or d >= self.nD):
			print "Dimention n="+str(d)+" does not exist!"
			sys.exit(1)
		if(self.Np == 0): return (0.,0.) 
		x_avg =  self.x_avg_v[d]/self.Np
		xp_avg =  self.xp_avg_v[d]/self.Np	
		return (x_avg,xp_avg)
		
	def getRmsU_UP(self,d):
		"""
		Returns the (rms u,rms up) parameters in the array 
		[u,up] for the dimension d.
		"""
		if(d <  0 or d >= self.nD):
			print "Dimention n="+str(d)+" does not exist!"
			sys.exit(1)
		if(self.Np == 0 or self.Np == 1): return (0.,0.) 
		x2_rms =  math.sqrt(math.fabs((self.x2_avg_v[d] -  (self.x_avg_v[d])**2/self.Np)/(self.Np-1)))
		xp2_rms = math.sqrt(math.fabs((self.xp2_avg_v[d] -  (self.xp_avg_v[d])**2/self.Np)/(self.Np-1))) 
		return (x2_rms,xp2_rms)
				
	def getMaxU_UP(self,d):
		"""
		Returns the (u_max,up_max) parameters in the array 
		[u,up] for the dimension d.
		"""
		if(d <  0 or d >= self.nD):
			print "Dimention n="+str(d)+" does not exist!"
			sys.exit(1)
		if(self.Np == 0): return (0.,0.) 
		return (self.x_max_v[d],self.xp_max_v[d])

	def getMinU_UP(self,d):
		"""
		Returns the (u_min,up_min) parameters in the array 
		[u,up] for the dimension d.
		"""
		if(d <  0 or d >= self.nD):
			print "Dimention n="+str(d)+" does not exist!"
			sys.exit(1)
		if(self.Np == 0): return (0.,0.) 
		return (self.x_min_v[d],self.xp_min_v[d])

