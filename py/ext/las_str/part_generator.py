import math,random

class TransverseCoordGen:
    """ Generates u and u' coordinates distributed according
    to the Gaussian with certain parameters emit_rms and
    beta: exp(-(u^2+beta^2*((alpha/beta)*u + u')^2)/(2*beta*emit_rms)) 
    alpha in [rad]
    u in [m] u' in [rad] beta in [m^2]
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
    def __init__(self,eKin,relativeSpread):
        self.mass = 0.938256 + 0.000511
        self.eKin = eKin
        self.relativeSpreadE = relativeSpread
        self.p0 = math.sqrt(math.pow(self.mass+eKin,2) - self.mass*self.mass)
        self.spreadP = (math.pow((self.mass + eKin),2)/self.p0)*self.relativeSpreadE
        
    def getP0(self):
        return self.p0

    def getEK0(self):
        return self.eKin
        
    def getMass(self):
        return self.mass
        
    def getP(self):
        return (self.p0+random.gauss(0.,math.sqrt(0.5))*self.spreadP)
        
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
        px = (xp - dpx)*p0
        py = yp*p0
        return (x,px,y,py,0.,pz)
        
#-----------------------------------------------------
#Generates bunch with certain parameters
#-----------------------------------------------------
