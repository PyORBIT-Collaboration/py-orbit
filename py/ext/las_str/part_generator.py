import math,random, time, sys, os, orbit_mpi

from bunch import *
from orbit_utils import *
from orbit_mpi import mpi_comm,mpi_datatype,mpi_op


class TransverseCoordGen:
    """ Generates u and u' coordinates distributed according
    to the Gaussian with certain parameters emit_rms and
    beta: exp(-(u^2+beta^2*((alpha/beta)*u + u')^2)/(2*beta*emit_rms)) 
    alpha in [rad]
    u in [m] u' in [rad] beta in [m]
    emit_rms in [m*rad]
    """
    
    def __init__(self, alpha,beta, emit_rms, u_cutoff):
        self.alpha = alpha
        self.beta = beta
        self.emit_rms = emit_rms
        self.u_cutoff = u_cutoff
        self.hamilton_max =  u_cutoff*u_cutoff
        self.parU = math.sqrt(1./(2*beta*emit_rms))
        self.parUP = math.sqrt(beta*beta/(2*beta*emit_rms))
        
    def getCoords(self):
        """ Returns u and u' [m] and [rad] """
        res = 0
        u = 0.
        up = 0.
        while(res == 0):
            u = random.gauss(0.,math.sqrt(0.5))/self.parU
            up = random.gauss(0.,math.sqrt(0.5))/self.parUP
            if( (u*u+self.beta*self.beta*up*up) <  self.hamilton_max):
                up = up - self.alpha*u/self.beta
                res = 1
        return (u,up)
                    
class EnergyGen:
    """
    It generates momentum distributed around P0. All values in GeV.
    """
    def __init__(self,eKin,relativeSpread,mass):
        self.mass = mass
        self.eKin = eKin
        self.relativeSpreadE = relativeSpread
        self.p0 = math.sqrt(math.pow(self.mass+eKin,2) - self.mass*self.mass)
        self.spreadP = ((self.mass + eKin)*eKin/self.p0)*self.relativeSpreadE
        print self.spreadP/self.p0

        
    def getP0(self):
        return self.p0

    def getEK0(self):
        return self.eKin
        
    def getMass(self):
        return self.mass
        
    def getP(self):
        return (self.p0+random.gauss(0.,1.)*self.spreadP)
        
class ParticlesGen:
    """
    It generates (x,px,y,py,z,pz) x,y,z in [m], px,py,pz in [GeV/c].
    Dispersion D in [m] and D' in [rad]
    trGenX,trGenY transverse generators TransverseCoordGen
    pGen - EnergyGen
    """
    def __init__(self,dispD,dispDP,trGenX,trGenY,pGen):
        self.dispD = dispD
        self.dispDP = dispDP
        self.trGenX = trGenX
        self.trGenY = trGenY
        self.pGen = pGen
    
    def getCoords(self):
        p0 = self.pGen.getP0()
        pz = self.pGen.getP()
        dp = pz - p0
        dx = self.dispD*dp/p0
        dpx = self.dispDP*dp/p0
        (x,xp) = self.trGenX.getCoords()
        (y,yp) = self.trGenY.getCoords()
        x = x + dx
        px = (xp + dpx)*p0
        py = yp*p0
        return (x,px,y,py,0.,pz)
        
#-----------------------------------------------------
#Generates bunch with certain parameters
#-----------------------------------------------------


class BunchGen:
    
    def getCoords(self, alpha, beta, emit_rms, u_cutoff):
        """ Returns u and u' [m] and [rad] """

        while 1:
            u = random.gauss(0.,math.sqrt(beta*emit_rms))
            up = random.gauss(0.,math.sqrt(emit_rms/beta))
            if(u*u+beta*beta*up*up < u_cutoff*u_cutoff):
                up = up - alpha*u/beta
                return u,up
                      

    def getBunch(self,time_par):        

        rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
        random.seed((rank+1)*12571+time_par*int(time.time()))
        
        E = self.mass + self.TK
        P = math.sqrt(E*E - self.mass*self.mass)
        #vz = 299792458*P/E
        
        #beta = P/E
        #gamma = E/self.mass
        
        #bg = beta*gamma
        
        bunch = Bunch()
        bunch.charge(self.charge)
        bunch.mass(self.mass)


        for i in range(self.N_part):

            x,xp = self.getCoords(self.alphaX, self.betaX, self.emtX, math.sqrt(self.emtX*self.betaX)*3.0)
            y,yp = self.getCoords(self.alphaY, self.betaY, self.emtY, math.sqrt(self.emtY*self.betaY)*3.0)
            z,zp = self.getCoords(self.alphaZ, self.betaZ, self.emtZ, math.sqrt(self.emtZ*self.betaZ)*3.0)             
            
            x += self.dispD*zp
            xp += self.dispDP*zp           

            px = xp*P
            py = yp*P
            pz = (1 + zp)*P
            
            bunch.addParticle(x,px,y,py,z,pz)
            
        
        return bunch



    def getAutoionizationBunch(self,mult,b,time_par):

        rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
        random.seed((rank+1)*12571+time_par*int(time.time()))


        N_evol = int(b.getPartAttrDicts()['Evolution']['size'])-5

        bunch = Bunch()
        bunch.charge(1)
        bunch.mass(0.938256)
        
        bunch_unstr = Bunch()
        bunch_unstr.charge(0)
        bunch_unstr.mass(b.mass())
        

        
        for i in range(b.getSize()):
            
            dt = b.partAttrValue("Evolution",i,N_evol+1)
            x0 = b.partAttrValue("Evolution",i,N_evol+2)
            y0 = b.partAttrValue("Evolution",i,N_evol+3)
            z0 = b.partAttrValue("Evolution",i,N_evol+4)
                
            (x1, y1, z1, px, py, pz) = (b.x(i), b.y(i), b.z(i), b.px(i), b.py(i), b.pz(i))
            p1 = b.partAttrValue("Evolution",i,N_evol)

                            
            for j in range(mult):

                ran = random.random()
                if(p1 <= ran):
                    bunch_unstr.addParticle(x1,px,y1,py,z1,pz)
                else:

                    for k in range(N_evol):
                        f1 = b.partAttrValue("Evolution",i,k)
                        f2 = b.partAttrValue("Evolution",i,k+1)
                        
                        if (f1 <= ran and ran < f2):
                            coeff = (k+(ran-f1)/(f2-f1))/N_evol               
                            (x,y,z) = (x0 + (x1 - x0)*coeff, y0 + (y1 - y0)*coeff, z0 + (z1 - z0)*coeff)  
                            bunch.addParticle(x,px,y,py,z,pz)
                            break

              
        return bunch, bunch_unstr





