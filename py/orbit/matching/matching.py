import os
import string
import sys
from numpy import *
from scipy.optimize import fsolve
from scipy.optimize import root
from scipy.integrate import odeint
from scipy.constants import c
from matplotlib.pyplot import *

from orbit.teapot import TEAPOT_MATRIX_Lattice



class Twiss:
	# Create a simple MAD-like twiss object:  
	def __init__(self):
		self.data = { 'keyword': '',
						's': 0.0,
						'L': 0.0,
						'alfx': 0.0,
						'alfy': 0.0,
						'betx': 0.0,
						'bety': 0.0,
						'mux' : 0.0,
						'muy' : 0.0,
						'Dx': 0.0,
						'Dpx': 0.0,
						'angle': 0.0, 
						'k1': 0.0 }


class Optics:
# An container class for twiss objects:
	def __init__(self):
		self.line = []

	def __len__(self):
		return len(self.line)

	def __getitem__(self,j):
		return self.line[j]

	def __setitem__(self,j,x):
		self.line[j]=x

	def add(self, x):
		self.line.append(x)

	def print_line(self):
		for j in xrange(0,len(self.line)):
			print j, self.line[j].data['keyword'], "s:", self.line[j].data['s'], "L:", self.line[j].data['L'], 360.0*self.line[j].data['mux'],self.line[j].data['bety'],self.line[j].data['alfy'] 

	def get_element(self, s):
		Nb=len(self.line)
		if self.line[0].data['s'] >= s and s >= 0.0: 
			return 0
		for j in xrange(1,Nb):
			if self.line[j-1].data['s'] < s and self.line[j].data['s'] >=s : 
				return j
			if self.line[Nb-1].data['s'] < s : 
				return 0
		if s < 0.0 : 
			return Nb-1
		else:
			print "error: s not in range"
			print "STOP."
			sys.exit(1)

	def get_length(self):
		Nb=len(self.line)
		return self.line[Nb-1].data['s']


		
	def readtwiss_teapot(self,lattice, bunch):
		
		beamline=Optics()
		matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,bunch)

		(arrmuX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
		(arrmuY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()

		(DispersionX, DispersionXP) = matrix_lattice.getRingDispersionDataX()
		(DispersionY, DispersionYP) = matrix_lattice.getRingDispersionDataY()
		
		nodes = lattice.getNodes()
		for node in nodes:
			for j in range(len(arrPosBetaX)):
				if (round(lattice.getNodePositionsDict()[node][1],4)==round(arrPosBetaX[j][0],4)):
					muX = arrmuX[j][1]
					betaX = arrPosBetaX[j][1]
					alphaX =  arrPosAlphaX[j][1]
					dx = DispersionX[j][1]
					dmux = DispersionXP[j][1]
					muY = arrmuY[j][1]
					betaY = arrPosBetaY[j][1]
					alphaY = arrPosAlphaY[j][1]
					dmuy = DispersionYP[j][1]
			if node.getType() == "quad teapot":
				k1l = node.getParam("kq")*node.getLength()
			else:
				k1l = 0.0
			if node.getType() == "bend teapot":
				angle = node.getParam("theta")
			else:
				angle = 0.0
			beamline.add(1)
			j=len(beamline)-1
			beamline[j]=Twiss()
			beamline[j].data['keyword']=node.getName()
			beamline[j].data['marker']=node.getType()
			beamline[j].data['s']=round(lattice.getNodePositionsDict()[node][1],4)
			beamline[j].data['L']=node.getLength()
			beamline[j].data['alfx']=alphaX
			beamline[j].data['alfy']=alphaY
			beamline[j].data['betx']=betaX
			beamline[j].data['bety']=betaY
			beamline[j].data['Dx']=dx
			beamline[j].data['Dpx']=dmux
			beamline[j].data['mux']=muX
			beamline[j].data['muy']=muY
			beamline[j].data['angle']=angle
			beamline[j].data['k1']=k1l
		return beamline
#------------------------------------------------------
# Read MADX TFS file
#-------------------------------------------------------



   
#------------------------------------------------------
# Envelope solver: 
# x0, xs0, y0, ys0: initial values
# emitx/y: rms emittance
# Ksc: space charge perveance
#-------------------------------------------------------
class EnvelopeSolver:
	
	def __init__(self,beamline):
		self.beamline = beamline


	def func_odeint(self,y,s,emitx,emity,sigma_p,Ksc):
		jb=self.beamline.get_element(s)
		k1=self.beamline[jb].data['k1']
		lj=self.beamline[jb].data['L']
		anglej=self.beamline[jb].data['angle']
		f0=y[1]
		f1=-(k1/lj+(anglej/lj)**2)*y[0]+emitx**2/y[0]**3+0.5*Ksc/(y[0]+y[2])+y[4]*sigma_p**2*anglej/(y[0]*lj) 
		f2=y[3]
		f3=(k1/lj)*y[2]+emity**2/y[2]**3+0.5*Ksc/(y[0]+y[2]) # -
		f4=y[5]
		f5=-(k1/lj+(anglej/lj)**2)*y[4]+0.5*Ksc/(y[0]*(y[0]+y[2]))*y[4]+anglej/lj
		return [f0,f1,f2,f3,f4,f5]
	
	def Dfunc_odeint(self,y,s,emitx,emity,sigma_p,Ksc):
		jb=self.beamline.get_element(s)
		k1=self.beamline[jb].data['k1']
		lj=self.beamline[jb].data['L']
		anglej=self.beamline[jb].data['angle']
		a0=-(k1/lj+(anglej/lj)**2)*y[0]+emitx**2/y[0]**3+0.5*Ksc/(y[0]+y[2])+y[4]*sigma_p**2*anglej/(y[0]*lj)
		a1=-(k1/lj+(anglej/lj)**2)*y[1]-3.0*y[1]*emitx**2/y[0]**4-0.5*Ksc*(y[1]+y[3])/(y[0]+y[2])**2+y[5]*sigma_p**2*anglej/(y[0]*lj)-y[4]*y[1]*sigma_p**2*anglej/(y[0]**2*lj) 
		a2=(k1/lj)*y[2]+emity**2/y[2]**3+0.5*Ksc/(y[0]+y[2]) # -
		a3=(k1/lj)*y[3]-3.0*y[3]*emity**2/y[2]**4-0.5*Ksc*(y[1]+y[3])/(y[0]+y[2])**2 # -
		a4=-(k1/lj+(anglej/lj)**2)*y[4]+0.5*Ksc/(y[0]*(y[0]+y[2]))*y[4]+anglej/lj
		a5=-(k1/lj+(anglej/lj)**2)*y[5]+0.5*Ksc/(y[0]*(y[0]+y[2]))*y[5]-0.5*Ksc/(y[0]*(y[0]+y[2]))**2*y[4]*(y[1]*(y[0]+y[2])+y[0]*(y[1]+y[3]) )
		return [a0,a1,a2,a3,a4,a5]

	
	def envelope_odeint(self, emitx, emity, sigma_p, Ksc, x0, xs0, y0, ys0, Dx0, Dxs0):  
		Np=1000	
		Nb=len(self.beamline)
		Lb=self.beamline[Nb-1].data['s']
		s=linspace(0.0,Lb,num=Np)
		sol=odeint(self.func_odeint,[x0,xs0,y0,ys0,Dx0,Dxs0],s,args=(emitx,emity,sigma_p,Ksc),Dfun=self.Dfunc_odeint,rtol=1.0e-12,atol=1.0e-12)
		envx=sol[:,0]
		envxs=sol[:,1]
		envy=sol[:,2]
		envys=sol[:,3]
		Dx=sol[:,4]
		Dxs=sol[:,5]
		return envx,envxs,envy,envys,Dx,Dxs,s



	#------------------------------------------------------
	# Match: Periodic solution starting from MADX result
	#-------------------------------------------------------
 
	# this is the function for the root searching routine (fsolve)
 
	def func_fsolve(self,x,emitx,emity,sigma_p,Ksc):
		envx,envxs,envy,envys,Dx,Dxs,s = self.envelope_odeint(emitx,emity,sigma_p,Ksc,x[0],x[1],x[2],x[3],x[4],x[5])
		Nb=len(envx)
		return [envx[Nb-1]-x[0],envxs[Nb-1]-x[1],envy[Nb-1]-x[2],envys[Nb-1]-x[3],Dx[Nb-1]-x[4],Dxs[Nb-1]-x[5]] 

	# root searching using fsolve and initial values from MADX
	# returns matched envelopes

	def match_root(self, emitx, emity, sigma_p, Ksc):
		Nb=len(self.beamline)
		# start values
		x0=sqrt(self.beamline[Nb-1].data['betx']*emitx)
		gamx=(1.0+(self.beamline[Nb-1].data['alfx'])**2)/self.beamline[Nb-1].data['betx']
		xs0=-copysign(sqrt(gamx*emitx),self.beamline[Nb-1].data['alfx'])
		y0=sqrt(self.beamline[Nb-1].data['bety']*emity)
		gamy=(1.0+(self.beamline[Nb-1].data['alfy'])**2)/self.beamline[Nb-1].data['bety']
		ys0=-copysign(sqrt(gamy*emity),self.beamline[Nb-1].data['alfy'])
		Dx0=self.beamline[Nb-1].data['Dx']
		Dxs0=self.beamline[Nb-1].data['Dpx']
		# solver
		sol = root(self.func_fsolve, [x0,xs0,y0,ys0,Dx0,Dxs0], args=(emitx,emity,sigma_p,Ksc),method='hybr')
		x0=sol.x[0]
		xs0=sol.x[1]
		y0=sol.x[2]
		ys0=sol.x[3]
		Dx0=sol.x[4]
		Dxs0=sol.x[5]
		envx,envxs,envy,envys,Dx,Dxs,s = self.envelope_odeint(emitx,emity,sigma_p,Ksc,x0,xs0,y0,ys0,Dx0,Dxs0)
		return envx, envxs, envy, envys, Dx, Dxs, s


	# returns the matchted twiss parameter at cell entrance 

	def match_twiss(self, emitx, emity, sigma_p, Ksc):
		Nb=len(self.beamline)
		# start values
		x0=sqrt(self.beamline[Nb-1].data['betx']*emitx)
		gamx=(1.0+(self.beamline[Nb-1].data['alfx'])**2)/self.beamline[Nb-1].data['betx']
		xs0=-copysign(sqrt(gamx*emitx),self.beamline[Nb-1].data['alfx'])
		y0=sqrt(self.beamline[Nb-1].data['bety']*emity)
		gamy=(1.0+(self.beamline[Nb-1].data['alfy'])**2)/self.beamline[Nb-1].data['bety']
		ys0=-copysign(sqrt(gamy*emity),self.beamline[Nb-1].data['alfy'])
		Dx0=self.beamline[Nb-1].data['Dx']
		Dxs0=self.beamline[Nb-1].data['Dpx']
		# solver
		sol = root(self.func_fsolve, [x0,xs0,y0,ys0,Dx0,Dxs0], args=(self.beamline,emitx,emity,sigma_p,Ksc),method='hybr')
		x0=sol.x[0]
		xs0=sol.x[1]
		y0=sol.x[2]
		ys0=sol.x[3]
		Dx0=sol.x[4]
		Dxs0=sol.x[5]
		return x0**2/emitx,y0**2/emity,-copysign(sqrt(x0**2*xs0**2/emitx**2),xs0),-copysign(sqrt(y0**2*ys0**2/emity**2),ys0), Dx0, Dxs0


	#------------------------------------------------------
	# Smooth focusing 
	#-------------------------------------------------------


	def func_smooth(self,x,phase0x,phase0y,length,emitx,emity,Ksc):
		kx=(phase0x/length)**2
		ky=(phase0y/length)**2
		return[emitx**2/x[0]**3-kx*x[0]+0.5*Ksc/(x[0]+x[1]),emity**2/x[1]**3-ky*x[1]+0.5*Ksc/(x[0]+x[1])]


	def match_smooth(self,phase0x,phase0y,length,emitx,emity,Ksc):
		kx=(phase0x/length)**2
		ky=(phase0y/length)**2

		x0=(emitx**2/kx)**(1.0/4.0) 
		y0=(emity**2/ky)**(1.0/4.0) 

		sol = root(self.func_smooth,[x0,y0],args=(phase0x,phase0y,length,emitx,emity,Ksc),method='hybr')
		
		return sol.x[0]**2/emitx,sol.x[1]**2/emity     # beta functions


	#------------------------------------------------------
	# Calculate phase advance for given envelopes
	#-------------------------------------------------------

	def phase_advance(self,envx,envy,Dx,emitx,emity,sigma_p,s):
		Np=len(s)
		phasex=0.0
		phasey=0.0
		ds=s[1]-s[0]
		for j in xrange(0,Np):
			phasex+=ds*emitx/(envx[j]**2-(Dx[j]*sigma_p)**2)
			phasey+=ds*emity/envy[j]**2
		return phasex, phasey	

	# analytic phase advance depression
	# lc: length of the cell

	def phase_analytic(self,emitx,emity,Ksc,lc):
		return 0.5*Ksc*lc/(4.0*emitx), 0.5*Ksc*lc/(4.0*emity) 


	#------------------------------------------------------
	# Entropy growth rate: pre-factor
	#-------------------------------------------------------

	def entropy_rate(self,envx,envy,emitx,emity,s,beta0):
		Np=len(s)
		ratet=0.0
		ds=s[1]-s[0]
		for j in xrange(0,Np):
			Tx=envx[j]**2/emitx**2
			Ty=envy[j]**2/emity**2
			ratet+=ds/(beta0*c)*0.5*(Tx-Ty)**2/(Tx*Ty)
		return ratet


