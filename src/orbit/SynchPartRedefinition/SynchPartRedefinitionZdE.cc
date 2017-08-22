#include "SynchPartRedefinitionZdE.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

#include "Bunch.hh"
#include "ParticleMacroSize.hh"
#include "SyncPart.hh"

#include "OrbitConst.hh"

/** Constructor */
SynchPartRedefinitionZdE::SynchPartRedefinitionZdE(): CppPyWrapper(NULL)
{
	
	z_dE_avg_arr = (double* ) malloc (2*sizeof(double));
	z_dE_avg_arr_MPI = (double* ) malloc (2*sizeof(double));

	for(int i = 0; i < 2; i++){
		z_dE_avg_arr[i] = 0.;
		z_dE_avg_arr_MPI[i] = 0.;		
	}
}

/** Destructor */
SynchPartRedefinitionZdE::~SynchPartRedefinitionZdE()
{
	free(z_dE_avg_arr);
	free(z_dE_avg_arr_MPI);
}

/** Performs the calculation of the z and dE averages of the bunch */		
void SynchPartRedefinitionZdE::analyzeBunch(Bunch* bunch){
	
	//initialization	
	for(int i = 0; i < 2; i++){
		z_dE_avg_arr[i] = 0.;
		z_dE_avg_arr_MPI[i] = 0.;		
	}	
	
	int count = 0;
	double total_macrosize = 0.;
	
	bunch->compress();
	double m_size = 0.;
	int nParts = bunch->getSize();
	count += nParts;
	double** part_coord_arr = bunch->coordArr();
	int has_msize = bunch->hasParticleAttributes("macrosize");
	if(has_msize > 0){
		ParticleMacroSize* macroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
		double m_size = 0.;
		for(int ip = 0; ip < nParts; ip++){
			m_size = macroSizeAttr->macrosize(ip);
			total_macrosize += m_size;
			for(int i = 0; i < 2; i++){
				z_dE_avg_arr[i] += m_size*part_coord_arr[ip][i+4];
			}		
		}	
	} else {
		m_size = 1.0;
		for(int ip = 0; ip < nParts; ip++){
			for(int i = 0; i < 2; i++){
				z_dE_avg_arr[i] += part_coord_arr[ip][i+4];
			}
		}
		total_macrosize += nParts*m_size;			
	}	
	
	int count_MPI = 0;
	ORBIT_MPI_Allreduce(&count,&count_MPI,1,MPI_INT,MPI_SUM,bunch->getMPI_Comm_Local()->comm);
	count = count_MPI;
	
	double total_macrosize_MPI = 0.;
	ORBIT_MPI_Allreduce(&total_macrosize,&total_macrosize_MPI,1,MPI_DOUBLE,MPI_SUM,bunch->getMPI_Comm_Local()->comm);
	total_macrosize = total_macrosize_MPI;
	
	ORBIT_MPI_Allreduce(z_dE_avg_arr,z_dE_avg_arr_MPI,2,MPI_DOUBLE,MPI_SUM,bunch->getMPI_Comm_Local()->comm);
	
	if(fabs(total_macrosize) > 0.){
		for(int i = 0; i < 2; i++){
			z_dE_avg_arr[i] = z_dE_avg_arr_MPI[i]/total_macrosize;
		}
	}	
}


/** Move the synch particle energy to the average energy */
void SynchPartRedefinitionZdE::centerE(Bunch* bunch){
	bunch->compress();
	int nParts = bunch->getSize();
	analyzeBunch(bunch);
	double delta_dE = getAvg_dE();
	shiftE(bunch,delta_dE);
}

/** Shift the synch particle energy */
void SynchPartRedefinitionZdE::shiftE(Bunch* bunch,double delta_dE){
	
	bunch->compress();
	
	int nParts = bunch->getSize();
	
	SyncPart* syncPart = bunch->getSyncPart();	
		
	double mass = syncPart->getMass();
		
	double momentum0 = syncPart->getMomentum();
	double beta0 = syncPart->getBeta();
	double eKin0 = syncPart->getEnergy();

	double eKin1 = eKin0 + delta_dE;
	double momentum1 = syncPart->energyToMomentum(eKin1);
	syncPart->setMomentum(momentum1);
	double beta1 = syncPart->getBeta();
	
	// the definition of xp = momentumX/momentum_synch_part
	double x_y_p_coeff = momentum0/momentum1;
	
	for(int ip = 0; ip < nParts; ip++){
		bunch->xp(ip) = x_y_p_coeff*bunch->xp(ip);
		bunch->yp(ip) = x_y_p_coeff*bunch->yp(ip);
		bunch->dE(ip) = bunch->dE(ip) - delta_dE;
	}
}

/** Move the synch particle's z position to the center of the bunch */
void SynchPartRedefinitionZdE::centerZ(Bunch* bunch){
	bunch->compress();
	int nParts = bunch->getSize();
	analyzeBunch(bunch);
	double delta_z = getAvg_Z();
	shiftZ(bunch,delta_z);
}

/** Shift the synch particle's z position */
void SynchPartRedefinitionZdE::shiftZ(Bunch* bunch,double delta_z){

	bunch->compress();
	
	int nParts = bunch->getSize();
	
	SyncPart* syncPart = bunch->getSyncPart();	
		
	double beta = syncPart->getBeta();

	// the z = c*beta*time (time is the same)
	double delta_time = delta_z/(OrbitConst::c*beta);
	
	syncPart->setTime(syncPart->getTime() - delta_time);
	
	for(int ip = 0; ip < nParts; ip++){
		bunch->z(ip) = bunch->z(ip) - delta_z;
	}
}

/** Returns the average z postion */
double SynchPartRedefinitionZdE::getAvg_Z(){
	return z_dE_avg_arr[0];
}

/** Returns the average dE value */
double SynchPartRedefinitionZdE::getAvg_dE(){
	return z_dE_avg_arr[1];
}


