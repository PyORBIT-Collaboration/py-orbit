from starkeffect import *
from ext.las_str.SimplexModMPF import Simplex
from mpmath import *
import sys









class Stark_calcSS: 
    
    

    def __init__(self,_n1, _n2, _m, point_ground_state):
        
        
        mp.prec = 10000
        self.n1 = _n1
        self.n2 = _n2
        self.m = abs(_m)
        self.n = _n1 + _n2 + abs(_m) + 1
        self.Energy = fdiv(-1,2*self.n*self.n)
        self.Z1 = mpf(fdiv(2*(2*self.n1 + abs(self.m) + 1), self.n))
        self.Z2 = fsub(4,self.Z1)
        self.F = mpf("0.0")
        self.max_err_exp = int(-point_ground_state*point_ground_state*0.22)
        self.pointN = 0
        self.pointM = self.def_point(point_ground_state)
        self.a = FuncSS(self.n1, self.n2, self.m, self.pointM)
        self.array_E = []
        self.calc = 1
        self.step_E = mpf("0.0")
        self.countE = 0
        self.tempE = mpf("1e100")
        self.err_E = mpf("1e-20")
        self.err_E_out = mpf("1e-20")




    def def_point(self, p1_for_ground):

        
        level0 = exp(p1_for_ground*p1_for_ground*mpf("-0.5"))
        for i in range(100,0,-1):            
            level = exp(i*i*mpf("-0.5")/self.n)*power(i,self.m)*hyp1f1(-self.n1, 1 + self.m, fdiv(i*i,self.n))
            if (abs(level)>level0):
                break
        return i + 1


    def absM(self,args):
        self.Z1 = args[0]
        mm = fabs(mpf(self.a.getM(str(self.F),str(self.Energy), str(self.Z1))))
#        print  "mm = ",nstr(mm,10)
        return mm
    
    
    def defr_parameters_forF(self):
        
        self.tempE = mpf("1e100")
        self.countE = 0
        
        mp.prec = self.a.calcPrecisionForM(self.max_err_exp,str(self.F),str(self.Energy), str(self.Z1))
        self.find_Z1()
        self.Z2 = fsub(4,self.Z1)
        self.pointN = self.a.calcPrecisionForN(str(self.F),str(self.Energy), str(self.Z2) )
        
        return


    
    def def_start_EG(self):
        
        n1 = self.n1
        n2 = self.n2
        m = self.m
        n = self.n
        F = self.F
        
        self.Energy = (
        -fdiv(1,2)/(n*n)
        +F*fdiv(3,2)*n*(n1-n2)
        -F*F*fdiv(1,16)*power(n,4)*(17*n*n-3*(n1-n2)*(n1-n2)-9*m*m+19)
        +F*F*F*fdiv(3,32)*power(n,7)*(n1-n2)*(23*n*n-(n1-n2)*(n1-n2)+11*m*m+39)
        -F*F*F*F*fdiv(1,1024)*power(n,10)*(5487*power(n,4)+35182*n*n+(5754-1134*m*m+1806*n*n)*(n1-n2)*(n1-n2)-549*m*m*m*m-3402*n*n*m*m+147*power(n1-n2,4)-8622*m*m+16211)
        +F*F*F*F*F*fdiv(3,1024)*power(n,13)*(n1-n2)*(10563*n*n*n*n+90708*n*n+(n1-n2)*(n1-n2)*(220*m*m+98*n*n+780)+772*n*n*m*m-21*power(n1-n2,4)+725*m*m*m*m+830*m*m+59293)
        )
        
        
        self.step_E = fabs(-fdiv(1,2)/(n*n) - self.Energy)*mpf("1e-10") 
        self.err_E = (fdiv(1,2)/(n*n))*mpf("1e-20")
        self.err_E_out = (fdiv(1,2)/(n*n))*mpf("1e-20")
        
        return
    
    

    
    def predict_E(self,k):
        
        
        n = len(self.array_E)
        m = min(n,k)
        self.Energy = 0
        for i in range(0,m):
            self.Energy += power(-1,m+i+1)*fac(m)/(fac(m-i)*fac(i))*self.array_E[i+n-m]
            
        
        n = len(self.array_E) - 1
        m = min(n,k)
        self.step_E = 0
        for i in range(0,m):
            self.step_E += power(-1,m+i+1)*fac(m)/(fac(m-i)*fac(i))*self.array_E[i+n-m]
        self.step_E = fabs(self.step_E - self.array_E[n]) 
        
        self.err_E = self.step_E*mpf("1e-5")
        self.err_E_out = (fdiv(1,2)/(self.n*self.n))*mpf("1e-20")
        

#        print  "Epred = ",nstr(self.Energy,100),"   step_E = ",nstr(self.step_E,10),"   err_E = ",nstr(self.err_E,10)

        return
          
            
             
             
    
    def initialEG(self):
         

        self.def_start_EG()    
        
        if (len(self.array_E)>4):
            
            min = mpf("1e100000")
            for k in range(3,50):
                self.predict_E(k)
                if (self.step_E < min):
                    min = self.step_E
                    k_min = k
                    
            self.predict_E(k_min)
                      
                 
        return





    def find_Z1(self):
        
        s = Simplex(self.absM, [self.Z1], [mpf("1e-5")])
        (values, err, iter) = s.minimize(power(10,self.max_err_exp), 10000000,0)
        self.Z1 = values[0]

        
        return 

         
  
  
    
    def absB0(self,args):
        
        
        self.Energy = args[0]     
        
        self.find_Z1()
        self.Z2 = fsub(4,self.Z1)
#        print  "Z2 = ",nstr(self.Z2,20)
        B = mpf(self.a.getB(str(self.F),str(self.Energy), str(self.Z2)))
        print  "E = ",nstr(self.Energy,50),"absB = ",nstr(fabs(B),20)
        
        if(fabs(self.tempE - self.Energy) < min(self.err_E,self.err_E_out)):
            self.countE += 1
        else:
            self.countE = 0

        
        
        if (self.countE > 4):
            return mpf("1e-1000000000")
        else:
            
            self.tempE = self.Energy

            
            return fabs(B)
    
  
    def TDM(self,args):
        
        self.Energy = args
        self.find_Z1()
        self.Z2 = fsub(4,self.Z1)
        
        C = mpf(self.a.C(str(self.F),str(self.Energy), str(self.Z2)))


        M = self.M
        N = self.N
        
        if (self.m == 1):
            V = quad(lambda mu:   exp(-(mu**2)/2)*M(mu)*(mu**3)*sqrt(mu)   ,[0,inf])*quad(lambda nu:   exp(-(nu**2)/2)*N(nu)*nu*sqrt(nu)   ,[0,inf]) + quad(lambda mu:   exp(-(mu**2)/2)*M(mu)*mu*sqrt(mu)   ,[0,inf])*quad(lambda nu:   exp(-(nu**2)/2)*N(nu)*(nu**3)*sqrt(nu)   ,[0,inf])
            V *= sqrt(mpf("2.0"))
            
        if (self.m == 0):
            V = quad(lambda mu:   exp(-(mu**2)/2)*M(mu)*(mu**4)*sqrt(mu)   ,[0,inf])*quad(lambda nu:   exp(-(nu**2)/2)*N(nu)*sqrt(nu)   ,[0,inf]) - quad(lambda mu:   exp(-(mu**2)/2)*M(mu)*sqrt(mu)   ,[0,inf])*quad(lambda nu:   exp(-(nu**2)/2)*N(nu)*(nu**4)*sqrt(nu)   ,[0,inf])

###    second order DTM
#        if (self.m == 1):
#            V = quad(lambda mu:   exp(-(mu**2)/2)*(sqrt(mpf("8"))-self.F*(sqrt(mpf("2"))*mu**2+mu**4/2/sqrt(mpf("2"))))*M(mu)*(mu**3)*sqrt(mu)   ,[0,inf])*quad(lambda nu:   exp(-(nu**2)/2)*(sqrt(mpf("8"))+self.F*(sqrt(mpf("2"))*nu**2+nu**4/2/sqrt(mpf("2"))))*N(nu)*nu*sqrt(nu)   ,[0,inf]) + quad(lambda mu:   exp(-(mu**2)/2)*(sqrt(mpf("8"))-self.F*(sqrt(mpf("2"))*mu**2+mu**4/2/sqrt(mpf("2"))))*M(mu)*mu*sqrt(mu)   ,[0,inf])*quad(lambda nu:   exp(-(nu**2)/2)*(sqrt(mpf("8"))+self.F*(sqrt(mpf("2"))*nu**2+nu**4/2/sqrt(mpf("2"))))*N(nu)*(nu**3)*sqrt(nu)   ,[0,inf])
#            V *= sqrt(mpf("2.0"))
            
#        if (self.m == 0):
#            V = quad(lambda mu:   exp(-(mu**2)/2)*(sqrt(mpf("8"))-self.F*(sqrt(mpf("2"))*mu**2+mu**4/2/sqrt(mpf("2"))))*M(mu)*(mu**4)*sqrt(mu)   ,[0,inf])*quad(lambda nu:   exp(-(nu**2)/2)*(sqrt(mpf("8"))+self.F*(sqrt(mpf("2"))*nu**2+nu**4/2/sqrt(mpf("2"))))*N(nu)*sqrt(nu)   ,[0,inf]) - quad(lambda mu:   exp(-(mu**2)/2)*(sqrt(mpf("8"))-self.F*(sqrt(mpf("2"))*mu**2+mu**4/2/sqrt(mpf("2"))))*M(mu)*sqrt(mu)   ,[0,inf])*quad(lambda nu:   exp(-(nu**2)/2)*(sqrt(mpf("8"))+self.F*(sqrt(mpf("2"))*nu**2+nu**4/2/sqrt(mpf("2"))))*N(nu)*(nu**4)*sqrt(nu)   ,[0,inf])
###    

######    orthogonality
#        if (self.m == 0):
#            V = quad(lambda mu:   exp(-(mu**2)/2)*(sqrt(mpf("8"))-self.F*(sqrt(mpf("2"))*mu**2+mu**4/2/sqrt(mpf("2"))))*M(mu)*(mu**2)*sqrt(mu)   ,[0,inf])*quad(lambda nu:   exp(-(nu**2)/2)*(sqrt(mpf("8"))+self.F*(sqrt(mpf("2"))*nu**2+nu**4/2/sqrt(mpf("2"))))*N(nu)*sqrt(nu)   ,[0,inf]) + quad(lambda mu:   exp(-(mu**2)/2)*(sqrt(mpf("8"))-self.F*(sqrt(mpf("2"))*mu**2+mu**4/2/sqrt(mpf("2"))))*M(mu)*sqrt(mu)   ,[0,inf])*quad(lambda nu:   exp(-(nu**2)/2)*(sqrt(mpf("8"))+self.F*(sqrt(mpf("2"))*nu**2+nu**4/2/sqrt(mpf("2"))))*N(nu)*(nu**2)*sqrt(nu)   ,[0,inf])
            
#        if (self.m == 1):
#            V = 0
######    orthogonality        
        
#        return fabs(C*sqrt(mp.pi)*V/8)        
        return fabs(C*sqrt(mp.pi)*V)
#        return fabs(2*C*sqrt(mp.pi)*V/8)
    
    
    
    
    def find_EG(self):
        
        Ein = self.Energy

        s = Simplex(self.absB0, [self.Energy], [self.step_E])     
            
        (values, err, iter) = s.minimize(mpf("1e-1000000"), 10000000,0)

 
        if(err == mpf("0.0")):
            self.calc = 0
            self.Energy = Ein
        else:
            self.array_E.append(self.Energy)
              
        return 
    
    
    
    def M(self,mu):
        
        return mpf(self.a.M(str(mu),str(self.F),str(self.Energy), str(self.Z1)))
        
    def N(self,nu):
        
        return mpf(self.a.N(str(nu),str(self.F),str(self.Energy), str(self.Z2)))        
    
    
    
    
    def readG(self,addr):
        
        f = open(addr + '%i%i%i.txt'%(self.n1,self.n2,self.m),'r')
        lines = f.readlines()
        f.close()
                
        return mpf((lines[int(self.F/mpf((lines[1].split())[0]))].split())[2])
 




    def main(self):
        
        self.initialEG()
 
        if(self.F != mpf("0.0")):
            self.calc = 1
            self.defr_parameters_forF()
            self.find_EG()
        else:
            self.calc = 0    
        

        

            
def dE_cont(continuum_spectr):
    
    lines = open(continuum_spectr).read().splitlines()
    max = 0
    for i in xrange(len(lines)):
        E,DTM = lines[i].split()
        if(float(DTM) >= max): 
            max = float(DTM) 
            dE = float(E)

    return dE



"""

class WaveFunc:
    
    def __init__(self,n1, n2, m, maxM, maxN,str_Energy,str_Gamma,str_reZ1, str_imZ1, str_F, num):
        
        mp.dps = max(len(str_Energy),len(str_Gamma),len(str_reZ1),len(str_imZ1),len(str_F))       
        self.psi = WaveFunction(n1,n2,m,maxM,maxN,mp.dps,"("+str_Energy+" "+ str(mpf(str_Gamma)*mpf("-0.5"))+")", "("+str_reZ1+" "+str_imZ1+")",str_F)
        self.mode = self.psi.getMode()
        self.switch = 0
        n = n1 + n2 + m + 1


        if (mpf(str_F) == mpf("0.0")):
            self.sqrt_norm = 1
            self.mode = 0
            mp.dps = 15
        else:
            z0 = 1/sqrt(mpf(str_F))
            mp.dps = 15
            
            if (num < 2):
                if (self.mode == 0):
                    self.sqrt_norm = sqrt(2*pi*quad(lambda mu: quad(lambda nu:   fabs(self.psi.getMN(mu, nu))**2     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM]))
                else:
                
                    sqrt_norm_pert = sqrt(2*pi*quad(lambda mu: quad(lambda nu:   fabs(self.psi.getMN(mu, nu))**2     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM]))
                    sqrt_norm_num = sqrt(2*pi*quad(lambda mu: quad(lambda nu:   fabs(self.M(mu)*self.N(nu))**2     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM]))
                    
                    sqrt_norm_num_1 = sqrt_norm_num*sqrt(fac(n1+m)*fac(n2+m)/(fac(n1)*fac(n2)))/(sqrt(pi)*power(n,2+m)*fac(m)*fac(m))
            self.mode = 0
            self.sqrt_norm = sqrt_norm_pert
                    
                    if (fabs(sqrt_norm_num_1 - 1) < mpf("0.01")):
                        self.switch = 1

                        
            else:
                self.sqrt_norm = sqrt(2*pi*quad(lambda mu: quad(lambda nu:   fabs(self.M(mu)*self.N(nu))**2     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM]))
                self.mode = 1
"""                        


"""
class WaveFunc:
    
    def __init__(self,n1, n2, m, maxM, maxN,str_Energy,str_Gamma,str_reZ1, str_imZ1, str_F, num):
        
        mp.dps = max(len(str_Energy),len(str_Gamma),len(str_reZ1),len(str_imZ1),len(str_F))       
        self.psi = WaveFunction(n1,n2,m,maxM,maxN,mp.dps,"("+str_Energy+" "+ str(mpf(str_Gamma)*mpf("-0.5"))+")", "("+str_reZ1+" "+str_imZ1+")",str_F)
        self.mode = self.psi.getMode()
        self.switch = 0
        n = n1 + n2 + m + 1


        if (mpf(str_F) == mpf("0.0")):
            self.sqrt_norm = 1
            self.mode = 0
            mp.dps = 15
        else:
            z0 = 1/sqrt(mpf(str_F))
            mp.dps = 15
            
            if (num < 2):
                if (self.mode == 0):
                    self.sqrt_norm = sqrt(2*pi*quad(lambda mu: quad(lambda nu:   fabs(self.psi.getMN(mu, nu))**2     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM]))
                else:
                
                    sqrt_norm_pert = sqrt(2*pi*quad(lambda mu: quad(lambda nu:   fabs(self.psi.getMN(mu, nu))**2     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM]))
                    sqrt_norm_num = sqrt(2*pi*quad(lambda mu: quad(lambda nu:   fabs(self.M(mu)*self.N(nu))**2     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM]))
                    
                    sqrt_norm_num_1 = sqrt_norm_num*sqrt(fac(n1+m)*fac(n2+m)/(fac(n1)*fac(n2)))/(sqrt(pi)*power(n,2+m)*fac(m)*fac(m))
                    self.mode = 0
                    self.sqrt_norm = sqrt_norm_pert
                    
                    if (fabs(sqrt_norm_num_1 - 1) < mpf("0.01")):
                        self.switch = 1

                        
            else:
                self.sqrt_norm = sqrt(2*pi*quad(lambda mu: quad(lambda nu:   fabs(self.M(mu)*self.N(nu))**2     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM]))
                self.mode = 1
                
"""

class WaveFunc:
    
    def __init__(self,n1, n2, m, maxM, maxN,str_Energy,str_Gamma,str_reZ1, str_imZ1, str_F, num):
        
        mp.dps = max(len(str_Energy),len(str_Gamma),len(str_reZ1),len(str_imZ1),len(str_F))       
        self.psi = WaveFunction(n1,n2,m,maxM,maxN,mp.dps,"("+str_Energy+" "+ str(mpf(str_Gamma)*mpf("-0.5"))+")", "("+str_reZ1+" "+str_imZ1+")",str_F)
        self.mode = self.psi.getMode()
        self.switch = 0
        self.mode = 0
        mp.dps = 15


        F_thr = (
        ((n1,n2,m) == (0,0,0))*mpf("0.000200")+
        ((n1,n2,m) == (0,0,1))*mpf("0.000040")+
        ((n1,n2,m) == (0,1,0))*mpf("0.000040")+
        ((n1,n2,m) == (0,1,1))*mpf("0.000022")+
        ((n1,n2,m) == (0,2,0))*mpf("0.000021")+
        ((n1,n2,m) == (0,2,1))*mpf("0.000015")+
        ((n1,n2,m) == (0,3,0))*mpf("0.000014")+
        ((n1,n2,m) == (0,3,1))*mpf("0.000011")+
        ((n1,n2,m) == (0,4,0))*mpf("0.000011")+
        ((n1,n2,m) == (1,0,0))*mpf("0.000041")+
        ((n1,n2,m) == (1,0,1))*mpf("0.000026")+
        ((n1,n2,m) == (1,1,0))*mpf("0.000024")+
        ((n1,n2,m) == (1,1,1))*mpf("0.000017")+
        ((n1,n2,m) == (1,2,0))*mpf("0.000017")+
        ((n1,n2,m) == (1,2,1))*mpf("0.000012")+
        ((n1,n2,m) == (1,3,0))*mpf("0.000012")+
        ((n1,n2,m) == (2,0,0))*mpf("0.000027")+
        ((n1,n2,m) == (2,0,1))*mpf("0.000020")+
        ((n1,n2,m) == (2,1,0))*mpf("0.000019")+
        ((n1,n2,m) == (2,1,1))*mpf("0.000013")+
        ((n1,n2,m) == (2,2,0))*mpf("0.000013")+
        ((n1,n2,m) == (3,0,0))*mpf("0.000021")+
        ((n1,n2,m) == (3,0,1))*mpf("0.000016")+
        ((n1,n2,m) == (3,1,0))*mpf("0.000015")+
        ((n1,n2,m) == (4,0,0))*mpf("0.000018")
        )
           
                 
        if (mpf(str_F) == mpf("0.0")):
            self.sqrt_norm = 1
        else:
            z0 = 1/sqrt(mpf(str_F))
            
            if (mpf(str_F) > F_thr):     
                self.mode = 1
                self.switch = 2
                self.sqrt_norm = sqrt(2*pi*quad(lambda mu: quad(lambda nu:   fabs(self.M(mu)*self.N(nu))**2     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM]))
            else:
                self.sqrt_norm = sqrt(2*pi*quad(lambda mu: quad(lambda nu:   fabs(self.psi.getMN(mu,nu))**2     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM]))
    
    
    

    
    def M(self,mu):
        
        par = self.psi.M(str(mu))
        line =  par.replace("(","").replace(")","").rsplit(" ")
     
        return mpc(line[0], line[1])
    
    def N(self,nu):
        
        par = self.psi.N(str(nu))
        line =  par.replace("(","").replace(")","").rsplit(" ")
     
        return mpc(line[0], line[1])
    
    
    def MN(self,mu,nu):
        
        if (self.mode==0):
            return self.psi.getMN(mu, nu)/self.sqrt_norm
        else:
            return self.M(mu)*self.N(nu)/self.sqrt_norm
    
    
    
    
