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
#from ext.las_str.StarkEffectSS  import *



       
class SchredingerFunc:    

    
    
    
    def __init__(self,method):
        

        self.TK = []
        self.n_sigma = []
        self.la = []
        self.n_step = []
        self.fx = []
        self.fy = []
        self.wx = []
        self.wy = []
        self.rx = []
        self.ry = []
        self.ax = []
        self.ay = []
        self.sigma_beam = 0
        self.k_stretch = 0
        self.env_simga = 0
        self.fS = 0
        self.power = []
        self.data_addr_name = []
        self.method = method
        
        self.bunch = Bunch()
        self.bunch_target = Bunch()
        self.count = 0
        
        self.cut_par = 1000
        self.n_states = 2
        
        if (self.method == 2):
            self.St = Stark(os.environ["ORBIT_ROOT"] + "/ext/laserstripping/Hydrogen_data/",self.n_states)
        if (self.method == 3):
            self.StSf = StarkStrongField(os.environ["ORBIT_ROOT"] + "/ext/laserstripping/Hydrogen_data/",0,0,1)
        if (self.method == 4):
            self.continuum_spectr = TDMcontinuum(os.environ["ORBIT_ROOT"] + "/ext/laserstripping/Hydrogen_data/TDMcontinuum/001/")
        
        
        self.EMfield = []
        self.Nevol = 0
        self.print_file = False
        
        self.Bx = 0
        self.By = 0
        self.Bz = 0








    
    def population(self):
        


        
        bunch_target = self.bunch_target
        mpi_size = orbit_mpi.MPI_Comm_size(mpi_comm.MPI_COMM_WORLD)
        

        bunch = self.bunch
        N_part = bunch.getSize()
        TK = self.TK
        n_step = self.n_step
        n_sigma = self.n_sigma
        power = self.power
        n_states = self.n_states
        cut_par = self.cut_par
        sigma_beam = self.sigma_beam
        env_sigma = self.env_sigma
        Nevol = self.Nevol
        print_file = self.print_file 
        
        if (self.method == 2):
            St = self.St
        if (self.method == 3):
            StSf = self.StSf
        if (self.method == 4):
            continuum_spectr = self.continuum_spectr
        
        fS = self.fS
        la = self.la
        fx = self.fx
        fy = self.fy
        wx = self.wx
        wy = self.wy 
               
        rx = self.rx
        ry = self.ry
        ax = self.ax
        ay = self.ay
        
        method = self.method
        Bx = self.Bx
        By = self.By
        Bz = self.Bz
        ####### Here are defined parameters of the function ###############       

  
        E = bunch.mass() + TK
        P = math.sqrt(E*E - bunch.mass()*bunch.mass())
        vz = 299792458*P/E
        
        fS = ConstEMfield(0.,0.,0.,Bx,0.,0.)                                  
                    
        bunch_target.deleteAllParticles()
        bunch.copyBunchTo(bunch_target)
        
        
        if (method == 1):
            dip_transition = math.sqrt(256*math.pow(n_states,7)*math.pow(n_states-1,2*n_states-5)/3/math.pow(n_states+1,2*n_states+5))
            delta_E = 1./2. - 1./(2.*n_states*n_states)
            
        if (method == 2):
            delta_E = 1./2. - 1./(2.*n_states*n_states)      
            
        if (method == 3):
            delta_E = StSf.deltaE(bunch.mass(),0.,0.,0.,Bx,By,Bz,0.,0.,P)
            
        if (method == 4):
            delta_E = continuum_spectr.setField_returndE(bunch.mass(),0.,0.,0.,Bx,By,Bz,0.,0.,P)
        
        
        la0 = 2*math.pi*5.291772108e-11/7.297352570e-3/delta_E
        te = TK - bunch.mass()*(la/la0-1)
        kz = te/math.sqrt(P*P-te*te)
        print "angle = ",math.atan2(1,-kz)*360/2/math.pi

    

        zb = -5*sigma_beam*vz
        zl = zb*E/P
        time_step = (2*abs(zb)/vz)/n_step



        for i in range(N_part):
            z = bunch_target.z(i)
            bunch_target.z(i,z + zb)   
        #bunch_target.dumpBunch("bunch_ini"+str(count)+".dat")
        LFS = HermiteGaussianLFmode(math.sqrt(power),0,0,abs(wx), abs(wy),fx,fy,la,zl,env_sigma)
        LFS.setLocalParameters(abs(rx), abs(ry),ax,ay)
        
        LFS.setLaserFieldOrientation(0.,0.,0.,   -1.,0.,kz,   kz,0.,1.,  kz,0.,1.)    #perpendicular polarization
#        LFS.setLaserFieldOrientation(0.,0.,0.,   -1.,0.,kz,   kz,0.,1.,  0.,1.,0.)      #parallel polarization
        tracker = RungeKuttaTracker(0)
        
        if (method == 1):   eff = TwoLevelAtom(LFS,delta_E,dip_transition)
        if (method == 2):   eff = SchrodingerEquation(LFS,St,cut_par)
        if (method == 3):   eff = TwoLevelStrongField(LFS, StSf)
        if (method == 4):   eff = ContinuumSS(LFS,continuum_spectr)
      
        cont_eff = ExtEffectsContainer()  
        cont_eff.AddEffect(eff)
        
        if(print_file):
            pr = PrintExtEffects("Populations",n_step,os.environ["ORBIT_ROOT"]+"/ext/laserstripping/working_dir/"+"/data3.0")
            cont_eff.AddEffect(pr)
        if(Nevol != 0):
            evo = RecordEvolution("Populations",0,Nevol)
            cont_eff.AddEffect(evo)
        
        tracker.track(bunch_target,0,time_step*n_step, time_step,fS,cont_eff)
        
#        for i in range(N_part):
#            z = bunch_target.z(i)
#            bunch_target.z(i,z - vz*time_step*n_step) 
#        tracker.track(bunch_target,0,time_step*n_step, time_step,fS,cont_eff)

        population = 0.
        population2 = 0.   
        p_ioniz = 0.
        p_ioniz2 = 0.
        for i in range(N_part):
            
            if (method != 4):
                val = 1 - bunch_target.partAttrValue("Populations",i,0) - bunch_target.partAttrValue("Populations",i,1)
                p_val = bunch_target.partAttrValue("Populations",i,0)
                population += val
                population2 += val*val
                p_ioniz += p_val
                p_ioniz2 += p_val*p_val
            else:
                val = 1 - bunch_target.partAttrValue("Populations",i,0)
                p_val = 1 - bunch_target.partAttrValue("Populations",i,0)
                population += val
                population2 += val*val
                p_ioniz += p_val
                p_ioniz2 += p_val*p_val
                
            
        op = mpi_op.MPI_SUM
        data_type = mpi_datatype.MPI_DOUBLE
        population = orbit_mpi.MPI_Allreduce(population,data_type,op,mpi_comm.MPI_COMM_WORLD)
        population2 = orbit_mpi.MPI_Allreduce(population2,data_type,op,mpi_comm.MPI_COMM_WORLD)
        p_ioniz = orbit_mpi.MPI_Allreduce(p_ioniz,data_type,op,mpi_comm.MPI_COMM_WORLD)
        p_ioniz2 = orbit_mpi.MPI_Allreduce(p_ioniz2,data_type,op,mpi_comm.MPI_COMM_WORLD)
        population = population/(mpi_size*N_part)
        p_ioniz = p_ioniz/(mpi_size*N_part)
        sigma_pop = 0.
        sigma_p_ioniz = 0.
        if(N_part*mpi_size > 1):
            sigma_pop = math.sqrt((population2 - N_part*mpi_size*population*population)/(N_part*mpi_size*(N_part*mpi_size - 1)))
            sigma_p_ioniz = math.sqrt((p_ioniz2 - N_part*mpi_size*p_ioniz*p_ioniz)/(N_part*mpi_size*(N_part*mpi_size - 1)))

        return population, sigma_pop, p_ioniz, sigma_p_ioniz

