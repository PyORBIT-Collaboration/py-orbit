#-----------------------------------------------------
#Track bunch with r and p through the external field 
# The field is 1 T and has direction (0,1,0)
#-----------------------------------------------------

import sys,math,os,orbit_mpi,random

from bunch import *
from trackerrk4 import *
from laserstripping import *
from orbit_utils import *

from orbit_mpi import mpi_comm,mpi_datatype,mpi_op



       
class SchredingerFunc:    

    
    
    
    def __init__(self,trans,n_states):
        
        
        self.TK = []
        self.n_sigma = []
        self.la = []
        self.n_step = []
        self.delta_E = []
        self.fx = []
        self.wx = []
        self.wy = []
        self.power = []
        self.data_addr_name = []
        
        self.bunch = Bunch()
        self.bunch_target = Bunch()
        self.count = 0
        
        self.Stark = HydrogenStarkParam(trans,n_states)
        self.levels = n_states*(1+n_states)*(1+2*n_states)/6
        self.By = []
        self.EMfield = []








    
    def population(self):
        


        
        bunch_target = self.bunch_target
        mpi_size = orbit_mpi.MPI_Comm_size(mpi_comm.MPI_COMM_WORLD)

        bunch = self.bunch
        N_part = bunch.getSize()
        TK = self.TK
        n_step = self.n_step
        delta_E = self.delta_E
        dip_transition = self.dip_transition
        n_sigma = self.n_sigma
        power = self.power
        levels = self.levels
        Stark = self.Stark
        fS = self.fS
        la = self.la
        fx = self.fx
        fy = self.fy
        wx = self.wx
        wy = self.wy
        ####### Here are defined parameters of the function ###############



  
        E = bunch.mass() + TK
        P = math.sqrt(E*E - bunch.mass()*bunch.mass())
        vz = 299792458*P/E
                    
        bunch_target.deleteAllParticles()
        bunch.copyBunchTo(bunch_target)
        
        la0 = 2*math.pi*5.291772108e-11/7.297352570e-3/delta_E
        kz = -1/math.sqrt(math.pow(P/(bunch.mass()*(la/la0-1)-TK),2)-1.)
    
        z0 = n_sigma*wx*math.sqrt(1+math.pow(fx*la/(wx*wx*math.pi),2))*math.sqrt(1+kz*kz)

        time_step = (2*z0/vz)/n_step
        #       alpha=360*math.acos((b.mass()*(la/la0-1)-TK)/P)/2/math.pi
        #       print alpha      
        for i in range(N_part):
            z = bunch_target.z(i)
            x = bunch_target.x(i)
            bunch_target.z(i,-z0 - kz*x)   
        #bunch_target.dumpBunch("bunch_ini"+str(count)+".dat")
        LFS = HermiteGaussianLFmode(math.sqrt(power),0,0,wx,wy,fx,fy,la)
        LFS.setLaserFieldOrientation(0.,0.,0.,   -1.,0.,kz,   1.,0.,1./kz,  0.,1.,0.)
        tracker = RungeKuttaTracker(0.000000001)
        First = SchrodingerEquation(LFS,Stark,1000.)
        tracker.track(bunch_target,0,time_step*n_step, time_step,fS,First)
#        bunch_target.dumpBunch("bunch_res"+str(1)+".dat")
        population = 0.
        population2 = 0.    
        for i in range(N_part):
            val = bunch_target.partAttrValue("Populations",i,0) - bunch_target.partAttrValue("Populations",i,1)
            population += val
            population2 += val*val
        op = mpi_op.MPI_SUM
        data_type = mpi_datatype.MPI_DOUBLE
        population = orbit_mpi.MPI_Allreduce(population,data_type,op,mpi_comm.MPI_COMM_WORLD)
        population2 = orbit_mpi.MPI_Allreduce(population2,data_type,op,mpi_comm.MPI_COMM_WORLD)
        population = population/(mpi_size*N_part)
        sigma_pop = 0.
        if(N_part*mpi_size > 1):
            sigma_pop = math.sqrt((population2 - N_part*mpi_size*population*population)/(N_part*mpi_size*(N_part*mpi_size - 1)))

        return population, sigma_pop
  
  
    def SetGroundStateBeam_ref(self,bunch_in):
        
        self.bunch = bunch_in
        bunch = self.bunch
        levels = self.levels
        
        if (bunch.hasPartAttr("Amplitudes") == 0): bunch.addPartAttr("Amplitudes",{"size":2*levels+8})
        if (bunch.hasPartAttr("Populations") == 0): bunch.addPartAttr("Populations",{"size":levels+1})
         
        for i in range(bunch.getSize()):
            bunch.partAttrValue("Amplitudes",i,1,1.0)
            
    def SetGroundStateBeam_copy(self,bunch_in):
        
        bunch_in.copyBunchTo(self.bunch)
        bunch = self.bunch
        levels = self.levels
        
        if (bunch.hasPartAttr("Amplitudes") == 0): bunch.addPartAttr("Amplitudes",{"size":2*levels+8})
        if (bunch.hasPartAttr("Populations") == 0): bunch.addPartAttr("Populations",{"size":levels+1})
         
        for i in range(bunch.getSize()):
            bunch.partAttrValue("Amplitudes",i,1,1.0)    
    



        

