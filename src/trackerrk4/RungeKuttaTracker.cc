//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   RungeKuttaTracker.cc
//
// CREATED
//    04/18/2008
//
// DESCRIPTION
//    A class that tracks relativistic particles in external magnetic
//    and electric fields by using 4-th order Runge Kutta method. 
//    The possible external interaction with something 
//    else (e.r. laser field) can be introduced.
//
//
///////////////////////////////////////////////////////////////////////////
#include "RungeKuttaTracker.hh"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace TrackerRK4;
using namespace OrbitUtils;
	
#include "OrbitConst.hh"
	
RungeKuttaTracker::RungeKuttaTracker(double lengthIn){
	
	length = lengthIn;
		
	c_light = OrbitConst::c;
		
	//spatial eps in [m]
	eps = 0.000001;
			
	//initial number of steps
	n_init = 20;
	n_steps = n_init;
	
	t_step = 0.;
	
	//set a*x + b*y + c*z + d = 0 equations for entrance and exit planes
	plEntrV = new double[4];
	plExitV = new double[4];
	
	
	plEntrV[0] = 0.;
	plEntrV[1] = 0.;
	plEntrV[2] = - 1.0/(length/2.0);
	plEntrV[3] = - 1.0;
	
	plExitV[0] = 0.;
	plExitV[1] = 0.;
	plExitV[2] = 1.0/(length/2.0);
	plExitV[3] = - 1.0;	
	
}


RungeKuttaTracker::~RungeKuttaTracker(){
	delete [] plEntrV;
	delete [] plExitV;
}

		
/** It sets the accuracy of the spatial resolution for tracking. */
void RungeKuttaTracker::setSpatialEps(double eps){
	this->eps = eps;
}

/** It returns the accuracy of the spatial resolution for tracking. */
double RungeKuttaTracker::getSpatialEps(){
	return eps;
}

/** It returns the time step size for tracking. */
double RungeKuttaTracker::getTimeStep(){
	return t_step;
}

/** It sets the number of initial steps in tracking. */
void RungeKuttaTracker::setInitialStepsNumber(int n_init){
	this->n_init = n_init;
	n_steps = n_init;
}

/** It returns the number of initial steps in tracking. */
int RungeKuttaTracker::getInitialStepsNumber(){
	return n_init;
}

/** It returns the number of steps during the synchronous particle tracking. */
int RungeKuttaTracker::getStepsNumber(){
	return n_steps;
}

/** It sets the approximate length of the region. */
void RungeKuttaTracker::setLength(double length){
	this->length = length;
}

/** It returns the approximate length of the region. */
double RungeKuttaTracker::getLength(){
	return length;
}

/** It sets (a*x+b*y+c*z+d=0) coefficients for the entrance plane. */
void RungeKuttaTracker::setEntrPlane(double a, double b, double c, double d){
	plEntrV[0] = a;
	plEntrV[1] = b;
	plEntrV[2] = c;
	plEntrV[3] = d;
}

/** It sets (a*x+b*y+c*z+d=0) coefficients for the exit plane. */
void RungeKuttaTracker::setExitPlane(double a, double b, double c, double d){
	plExitV[0] = a;
	plExitV[1] = b;
	plExitV[2] = c;
	plExitV[3] = d;
}

/** It returns (a*x+b*y+c*z+d=0) coefficients for the entrance plane. */
void RungeKuttaTracker::getEntrPlane(double& a, double& b, double& c, double& d){
	a = plEntrV[0];
	b = plEntrV[1];
	c = plEntrV[2];
	d = plEntrV[3];
}

/** It sets (a*x+b*y+c*z+d=0) coefficients for the exit plane. */
void RungeKuttaTracker::getExitPlane(double& a, double& b, double& c, double& d){
	a = plExitV[0];
	b = plExitV[1];
	c = plExitV[2];
	d = plExitV[3];
}


/** It tracks a traditional ORBIT bunch. The external effects instance could be NULL. */
void RungeKuttaTracker::trackBunch(Bunch* bunch, BaseFieldSource* fieldSource, ExternalEffects* extEff){
	n_steps = n_init;
	int nMax = 5*n_steps;
	charge = bunch->getCharge();
	mass = bunch->getMass();
	mass2 = mass*mass;
	SyncPart* syncPart = bunch->getSyncPart();
	t_step = length/(n_steps*c_light*syncPart->getBeta());
	//find the entrance point for syncPart
	double cross_t = plEntrV[0]*syncPart->getPX() + 
	                 plEntrV[1]*syncPart->getPY() + 
	                 plEntrV[2]*syncPart->getPZ();
	if(cross_t >= 0.){
		ORBIT_MPI_Finalize("RungeKuttaTracker::trackBunch - syncPart cannot cross the entrance plane.");
	}
	cross_t = -(plEntrV[0]*syncPart->getX() + 
		plEntrV[1]*syncPart->getY() + 
		plEntrV[2]*syncPart->getZ() + 
		plEntrV[3])/cross_t;
	//these are coordinates of the syncPart entrance through the entrance plane plEntrV[0...3]
	y_init_vct[0] = cross_t*syncPart->getPX() + syncPart->getX();
	y_init_vct[1] = cross_t*syncPart->getPY() + syncPart->getY();
	y_init_vct[2] = cross_t*syncPart->getPZ() + syncPart->getZ();
	//std::cout<<"debug sync part x0="<<syncPart->getX()<<" y0="<<syncPart->getY()<<" z0="<<syncPart->getZ()<<std::endl;
	//std::cout<<"debug entry     x0="<<y_init_vct[0]<<" y0="<<y_init_vct[1]<<" z0="<<y_init_vct[2]<<std::endl;
	syncPart->setXYZ(y_init_vct);
	y_init_vct[3] = syncPart->getPX();
	y_init_vct[4] = syncPart->getPY();
	y_init_vct[5] = syncPart->getPZ();
	//this is a time when synch particle crosses the entrance plane
	double t_sync_final = 0.;
	double t_previous_final = 0.;
	double y_sync_final_vct[6];
	//step size definition from synch. particle
	double t_st = t_step;
	//std::cerr<<"debug length [m]="<< length <<std::endl;
	//std::cerr<<"debug c_light="<< c_light <<std::endl;
	//std::cerr<<"debug beta="<< syncPart->getBeta() <<std::endl;
	//std::cerr<<"debug time step [sec]="<< t_st <<std::endl;
	for(int iN = 1; iN < 16; iN++){
		//std::cerr<<"debug start loop iN=" <<iN<<" time step="<<t_st<<std::endl;
		nMax = nMax*2;
		for(int i = 0; i < 6; i++){
			y_in_vct[i] = y_init_vct[i];
			y_out_vct[i] = y_init_vct[i];
		}
		int step_count = 0;
		double t = 0.;
		//std::cerr<<"debug is before exit="<< isBeforeExit(y_out_vct) <<std::endl;
		while(isBeforeExit(y_out_vct)){
			step_count = step_count + 1;
			for(int i = 0; i < 6; i++){
				y_in_vct[i] = y_out_vct[i];
			}		
			//std::cerr<<"debug step="<< step_count <<" r=" << y_in_vct[0] <<" "<< y_in_vct[1] <<" "<< y_in_vct[2] <<std::endl;
      //std::cerr<<"debug  "<< step_count <<"   " << y_in_vct[0] <<" "<< y_in_vct[1] <<" "<< y_in_vct[2] <<std::endl;
			rk4Step(t,t_st,fieldSource);
			t = t + t_st;
			if(step_count > nMax){
				ORBIT_MPI_Finalize("RungeKuttaTracker::trackBunch - syncPart cannot exit outside.");
			}
		}
		//std::cerr<<"debug Exit Pl=" << plExitV[0] <<" "<< plExitV[1] <<" "<< plExitV[2]<<" "<< plExitV[3] <<std::endl;
		//interpolate the crossing point
		cross_t = (y_out_vct[0] - y_in_vct[0])*plExitV[0] +
		          (y_out_vct[1] - y_in_vct[1])*plExitV[1] +
		          (y_out_vct[2] - y_in_vct[2])*plExitV[2];
		cross_t = -(plExitV[0]*y_in_vct[0] + 
		            plExitV[1]*y_in_vct[1] + 
			          plExitV[2]*y_in_vct[2] + plExitV[3])/cross_t;	
		for(int i = 0; i < 6; i++){
			y_final_vct[i] = cross_t* (y_out_vct[i] - y_in_vct[i])+ y_in_vct[i];
		}	
		t_sync_final = t + t_st*(cross_t - 1.);
		//std::cerr<<"debug final r=" << y_final_vct[0] <<" "<< y_final_vct[1] <<" "<< y_final_vct[2] <<std::endl;
		//std::cerr<<"debug final p=" << y_final_vct[3] <<" "<< y_final_vct[4] <<" "<< y_final_vct[5] <<std::endl;
		if(iN > 1){
			double diff = 0.;
			for(int i = 0; i < 3; i++){
				diff = diff + (y_final_vct[i]-y_sync_final_vct[i])*(y_final_vct[i]-y_sync_final_vct[i]);
			}
			//the spatial difference
			diff = sqrt(diff);
			//std::cerr<<"debug spatial diff=" <<diff<<std::endl;
			//add the time difference
			double p_out_2 = y_final_vct[3]*y_final_vct[3]+y_final_vct[4]*y_final_vct[4]+y_final_vct[5]*y_final_vct[5];
			diff += fabs(t_sync_final-t_previous_final)*c_light*sqrt(p_out_2/(mass2+p_out_2));
			//std::cerr<<"debug with time diff=" <<diff<<std::endl;
			if(diff < eps){
				//that is enough!
			  t_step = t_st;
				n_steps = step_count;
				nMax = 10*step_count;
				for(int i = 0; i < 6; i++){
					y_sync_final_vct[i] = y_final_vct[i];
				}				
				break;
			}
		}
		//remember the last one
		for(int i = 0; i < 6; i++){
			y_sync_final_vct[i] = y_final_vct[i];
		}		
		t_st = t_st/2;
		t_previous_final = t_sync_final;
	}
	//std::cerr<<"debug track of synch part FINALIZED =======================" <<std::endl;	
	syncPart->setTime(syncPart->getTime() + t_sync_final);
  //std::cerr<<"debug time step [sec]="<< t_step <<std::endl;
	//std::cerr<<"debug spatial eps    ="<< eps <<std::endl;
  //std::cerr<<"debug number of steps="<< n_steps <<std::endl;
  //std::cerr<<"debug t_sync_final="<< t_sync_final <<std::endl;	
	//------------------------------------------------
	//performs the necessary actions before tracking
	if(extEff != NULL) extEff->setupEffects(bunch);
	//------------------------------------------------
	//move from x, x', y, y', z, dE to x,px,y,py,z,pz 
	bunch->compress();
	double** partCoordArr = bunch->coordArr();
	double p_start_vct[3];
	double r_start_vct[3];
	double p_final_vct[3];
	double r_final_vct[3];
	double nx_start_vct[3];
	double ny_start_vct[3];
	double nx_final_vct[3];
	double ny_final_vct[3];
	double p_sync_start_vct[3];
	double p_norm_sync_start_vct[3];
	double p_sync_final_vct[3];
	double p_norm_sync_final_vct[3];
	double r_sync_start_vct[3];	
	double r_sync_final_vct[3];	
	double v_vct[3];
	//set up start parameters for syncPart
	double v_sync_start = c_light*syncPart->getBeta();
	double energy_sync_start = syncPart->getEnergy();
	double pSyncPart_start = syncPart->getMomentum();	
	nx_start_vct[0] = syncPart->getNormalXX();
	nx_start_vct[1] = syncPart->getNormalXY();
	nx_start_vct[2] = syncPart->getNormalXZ();
	ny_start_vct[0] = syncPart->getNormalYX();
	ny_start_vct[1] = syncPart->getNormalYY();
	ny_start_vct[2] = syncPart->getNormalYZ();
	//std::cerr<<"debug orts nx="<< nx_start_vct[0] <<" " << nx_start_vct[1] <<" " << nx_start_vct[2]<<" "  <<std::endl;
	//std::cerr<<"debug orts ny="<< ny_start_vct[0] <<" " << ny_start_vct[1] <<" " << ny_start_vct[2]<<" "  <<std::endl;
	p_sync_start_vct[0] = syncPart->getPX();
	p_sync_start_vct[1] = syncPart->getPY();
	p_sync_start_vct[2] = syncPart->getPZ();	
	p_norm_sync_start_vct[0] = p_sync_start_vct[0]/pSyncPart_start;
	p_norm_sync_start_vct[1] = p_sync_start_vct[1]/pSyncPart_start;
	p_norm_sync_start_vct[2] = p_sync_start_vct[2]/pSyncPart_start;
	r_sync_start_vct[0] = syncPart->getX();
	r_sync_start_vct[1] = syncPart->getY();
	r_sync_start_vct[2] = syncPart->getZ();	
	//set up the new SyncParticle coordinates
	syncPart->setXYZ(y_sync_final_vct[0],y_sync_final_vct[1],y_sync_final_vct[2]);
	syncPart->setPXYZ(y_sync_final_vct[3],y_sync_final_vct[4],y_sync_final_vct[5]);	
	double v_sync_final = c_light*syncPart->getBeta();
	double energy_sync_final = syncPart->getEnergy();
	double pSyncPart_final = syncPart->getMomentum();
	nx_final_vct[0] = syncPart->getNormalXX();
	nx_final_vct[1] = syncPart->getNormalXY();
	nx_final_vct[2] = syncPart->getNormalXZ();
	ny_final_vct[0] = syncPart->getNormalYX();
	ny_final_vct[1] = syncPart->getNormalYY();
	ny_final_vct[2] = syncPart->getNormalYZ();
	p_sync_final_vct[0] = syncPart->getPX();
	p_sync_final_vct[1] = syncPart->getPY();
	p_sync_final_vct[2] = syncPart->getPZ();	
	p_norm_sync_final_vct[0] = p_sync_final_vct[0]/pSyncPart_final;
	p_norm_sync_final_vct[1] = p_sync_final_vct[1]/pSyncPart_final;
	p_norm_sync_final_vct[2] = p_sync_final_vct[2]/pSyncPart_final;
	r_sync_final_vct[0] = syncPart->getX();
	r_sync_final_vct[1] = syncPart->getY();
	r_sync_final_vct[2] = syncPart->getZ();	
	v_sync_final = c_light*syncPart->getBeta();	
	double t_start = 0.;
	double t_final = 0.;
	double pz_along = 0.;
	double enrg_part = 0.;
	double p_part = 0.;
	for(int ip = 0, nParts = bunch->getSize(); ip < nParts; ip++){
		//find the transverse (relative to the synch. particle moment) moment vector
		p_start_vct[0] = pSyncPart_start*(partCoordArr[ip][1]*nx_start_vct[0] + partCoordArr[ip][3]*ny_start_vct[0]);
		p_start_vct[1] = pSyncPart_start*(partCoordArr[ip][1]*nx_start_vct[1] + partCoordArr[ip][3]*ny_start_vct[1]);
		p_start_vct[2] = pSyncPart_start*(partCoordArr[ip][1]*nx_start_vct[2] + partCoordArr[ip][3]*ny_start_vct[2]);
		pz_along = (energy_sync_start+partCoordArr[ip][5])*(energy_sync_start+partCoordArr[ip][5] + 2*mass);
		pz_along = pz_along - (partCoordArr[ip][1]*pSyncPart_start)*(partCoordArr[ip][1]*pSyncPart_start)
		                    - (partCoordArr[ip][3]*pSyncPart_start)*(partCoordArr[ip][3]*pSyncPart_start);
		pz_along = sqrt(pz_along);
		p_start_vct[0] = p_start_vct[0] + pz_along*p_norm_sync_start_vct[0];
		p_start_vct[1] = p_start_vct[1] + pz_along*p_norm_sync_start_vct[1];
		p_start_vct[2] = p_start_vct[2] + pz_along*p_norm_sync_start_vct[2];	
		//now absolute position 
		r_start_vct[0] = r_sync_start_vct[0] + partCoordArr[ip][0]*nx_start_vct[0] + partCoordArr[ip][2]*ny_start_vct[0];
		r_start_vct[1] = r_sync_start_vct[1] + partCoordArr[ip][0]*nx_start_vct[1] + partCoordArr[ip][2]*ny_start_vct[1];
		r_start_vct[2] = r_sync_start_vct[2] + partCoordArr[ip][0]*nx_start_vct[2] + partCoordArr[ip][2]*ny_start_vct[2];
		p_part = sqrt(p_start_vct[0]*p_start_vct[0] + p_start_vct[1]*p_start_vct[1] + p_start_vct[2]*p_start_vct[2]);
		enrg_part = syncPart->momentumToEnergy(p_part);
		v_vct[0] = c_light*p_start_vct[0]/enrg_part;
		v_vct[1] = c_light*p_start_vct[1]/enrg_part;
		v_vct[2] = c_light*p_start_vct[2]/enrg_part;
		//t_start has been used, we can use it again 
		//find the entrance point and time
	  cross_t = plEntrV[0]*v_vct[0] + plEntrV[1]*v_vct[1] + plEntrV[2]*v_vct[2];
	  t_start = -(plEntrV[0]*r_start_vct[0] + plEntrV[1]*r_start_vct[1] + 
		            plEntrV[2]*r_start_vct[2] + plEntrV[3])/cross_t;	
		r_start_vct[0] = t_start*v_vct[0] + r_start_vct[0];
		r_start_vct[1] = t_start*v_vct[1] + r_start_vct[1];
		r_start_vct[2] = t_start*v_vct[2] + r_start_vct[2];
		//std::cerr<<"debug t_start="<< t_start <<std::endl;	
		t_start = t_start - partCoordArr[ip][4]/v_sync_start;
		//std::cerr<<"debug r_sync_start="<< r_sync_start_vct[0] <<" "<< r_sync_start_vct[1] <<" "<< r_sync_start_vct[2] <<std::endl;
		//std::cerr<<"debug p_sync_start="<< p_sync_start_vct[0]<<" "<< p_sync_start_vct[1]<<" "<< p_sync_start_vct[2] <<std::endl;
		//std::cerr<<"debug r_start="<< r_start_vct[0]<<" "<< r_start_vct[1]<<" "<< r_start_vct[2] <<std::endl;
		//std::cerr<<"debug p_start="<< p_start_vct[0]<<" "<< p_start_vct[1]<<" "<< p_start_vct[2] <<std::endl;		
		//start tracking
	  y_in_vct[0] = r_start_vct[0];
	  y_in_vct[1] = r_start_vct[1];
	  y_in_vct[2] = r_start_vct[2];
	  y_in_vct[3] = p_start_vct[0];
	  y_in_vct[4] = p_start_vct[1];
	  y_in_vct[5] = p_start_vct[2];	
		for(int i = 0; i < 6; i++){
			y_out_vct[i] = y_in_vct[i];
		}
		int step_count = 0;
		double t = t_start;
		while(isBeforeExit(y_out_vct) && bunch->flag(ip) > 0){
			partCoordArr[ip][0] = y_in_vct[0];
			partCoordArr[ip][2] = y_in_vct[1];
			partCoordArr[ip][4] = y_in_vct[2];
			partCoordArr[ip][1] = y_in_vct[3];
			partCoordArr[ip][3] = y_in_vct[4];
			partCoordArr[ip][5] = y_in_vct[5];
			step_count = step_count + 1;
			for(int i = 0; i < 6; i++){
				y_in_vct[i] = y_out_vct[i];
			}		
			rk4Step(t,t_step,fieldSource);
			//place to account for external effeects
			if(extEff != NULL) extEff->applyEffectsForEach(bunch, ip, y_in_vct, y_out_vct, t, t_step, fieldSource, this);
			t = t + t_step;
			if(step_count > nMax){
				ORBIT_MPI_Finalize("RungeKuttaTracker::trackBunch - particle cannot exit outside.");
			}
		}
		if(bunch->flag(ip) > 0){
			//interpolate the crossing point
			cross_t = (y_out_vct[0] - y_in_vct[0])*plExitV[0] +
			          (y_out_vct[1] - y_in_vct[1])*plExitV[1] +
			          (y_out_vct[2] - y_in_vct[2])*plExitV[2];
			cross_t = -(plExitV[0]*y_in_vct[0] + 
				          plExitV[1]*y_in_vct[1] + 
				          plExitV[2]*y_in_vct[2] + plExitV[3])/cross_t;	
			for(int i = 0; i < 3; i++){
				r_final_vct[i] = cross_t* (y_out_vct[i] - y_in_vct[i])+ y_in_vct[i];
				//std::cerr<<"debug r_vct["<< i<<"]="<< r_final_vct[i]<<std::endl;
			}
			for(int i = 0; i < 3; i++){
				p_final_vct[i] = cross_t* (y_out_vct[i+3] - y_in_vct[i+3])+ y_in_vct[i+3];
				//std::cerr<<"debug p_vct["<< i<<"]="<< p_final_vct[i]<<std::endl;
			}
			//for(int i = 0; i < 6; i++){
			//	std::cerr<<"debug p_sync["<< i<<"]="<< y_sync_final_vct[i]<<std::endl;
			//}
			for(int i = 0; i < 3; i++){
				r_final_vct[i] = (r_final_vct[i] - r_sync_final_vct[i]);
			}
			t_final = t + t_step*(cross_t - 1.);
			//std::cerr<<"debug t_final="<< t_final <<std::endl;
			//tracking is finished, let's put final coordinates 
		  p_part = sqrt(p_final_vct[0]*p_final_vct[0] + p_final_vct[1]*p_final_vct[1] + p_final_vct[2]*p_final_vct[2]);
		  enrg_part = syncPart->momentumToEnergy(p_part);
			partCoordArr[ip][5] = enrg_part - energy_sync_final;
			partCoordArr[ip][1] = (p_final_vct[0]*nx_final_vct[0] + 
				                     p_final_vct[1]*nx_final_vct[1] + 
														 p_final_vct[2]*nx_final_vct[2])/pSyncPart_final;
			partCoordArr[ip][3] = (p_final_vct[0]*ny_final_vct[0] + 
				                     p_final_vct[1]*ny_final_vct[1] + 
														 p_final_vct[2]*ny_final_vct[2])/pSyncPart_final;
			partCoordArr[ip][0] = (r_final_vct[0]*nx_final_vct[0] + 
				                     r_final_vct[1]*nx_final_vct[1] + 
														 r_final_vct[2]*nx_final_vct[2]);
			partCoordArr[ip][2] = (r_final_vct[0]*ny_final_vct[0] + 
				                     r_final_vct[1]*ny_final_vct[1] + 
														 r_final_vct[2]*ny_final_vct[2]);
			partCoordArr[ip][4] = (-(t_final - t_sync_final))*v_sync_final +
				                    (r_final_vct[0]*p_norm_sync_final_vct[0] +
				                     r_final_vct[1]*p_norm_sync_final_vct[1] +
														 r_final_vct[2]*p_norm_sync_final_vct[2]);
		}
	}
	//------------------------------------------------
	//performs the necessary actions after tracking
	if(extEff != NULL) extEff->finalizeEffects(bunch);
	//------------------------------------------------	
	//remove the dead particles
	bunch->compress();
	/** Test part
  std::cerr<<"debug time step [sec]="<< t_step <<std::endl;
  std::cerr<<"debug number of steps="<< n_steps <<std::endl;
	double x = 2. ,y = 2., z = 2. , t = 3.;
	double fx,fy,fz;
	fieldSource->getElectricField(x,y,z,t,fx,fy,fz);
	std::cerr<<"debug Electr. field ="<< fx <<" "<< fy <<" "<< fz <<std::endl;
  fieldSource->getMagneticField(x,y,z,t,fx,fy,fz);
	std::cerr<<"debug Magnet. field ="<< fx <<" "<< fy <<" "<< fz <<std::endl;
	*/
}

void RungeKuttaTracker::track(Bunch* bunch,double t_begin, double t_period, double t_step_in, 
	                            BaseFieldSource* fieldSource, ExternalEffects* extEff)
{
	charge = bunch->getCharge();
	mass = bunch->getMass();
	mass2 = mass*mass;
	n_steps = int(t_period/t_step_in + 0.5);
	if(n_steps < 1){ n_steps = 1;}
	t_step = t_period/n_steps;	
	double vpt = 0.;
	double** partCoordArr = bunch->coordArr();
	int flag = 0;
	double t = t_begin;
	//------------------------------------------------
	//performs the necessary actions before tracking
	if(extEff != NULL) extEff->setupEffects(bunch);	
	//------------------------------------------------	
	for (int i = 0; i < n_steps; i++){
		
		t = t_begin + i*t_step;
		
		if(extEff != NULL) extEff->memorizeInitParams(bunch);
		
		for(int ip = 0, nParts = bunch->getSize(); ip < nParts; ip++){
			flag = bunch->flag(ip);
			if(flag > 0){
				y_in_vct[0] = partCoordArr[ip][0];
				y_in_vct[1] = partCoordArr[ip][2];
				y_in_vct[2] = partCoordArr[ip][4];
				y_in_vct[3] = partCoordArr[ip][1];
				y_in_vct[4] = partCoordArr[ip][3];
				y_in_vct[5] = partCoordArr[ip][5];
				if(!isOutside(y_in_vct)){
					rk4Step(t,t_step,fieldSource);	
				} else {
					vpt = mass2 + y_in_vct[3]*y_in_vct[3] + y_in_vct[4]*y_in_vct[4] + y_in_vct[5]*y_in_vct[5];
					vpt = c_light*t_step/sqrt(vpt);
					y_out_vct[0] = y_in_vct[0] + vpt*y_in_vct[3];
					y_out_vct[1] = y_in_vct[1] + vpt*y_in_vct[4];
					y_out_vct[2] = y_in_vct[2] + vpt*y_in_vct[5];
					y_out_vct[3] = y_in_vct[3];
					y_out_vct[4] = y_in_vct[4];
					y_out_vct[5] = y_in_vct[5];
				}
			  partCoordArr[ip][0] = y_out_vct[0];
			  partCoordArr[ip][2] = y_out_vct[1];
			  partCoordArr[ip][4] = y_out_vct[2];
			  partCoordArr[ip][1] = y_out_vct[3];
			  partCoordArr[ip][3] = y_out_vct[4];
			  partCoordArr[ip][5] = y_out_vct[5];	
			}
			if(extEff != NULL) extEff->applyEffectsForEach(bunch, ip,y_in_vct , y_out_vct, t, t_step, fieldSource, this);
		}
		//apply the external effects
		if(extEff != NULL) extEff->applyEffects(bunch, t, t_step, fieldSource, this);
		
	}
	//------------------------------------------------
	//performs the necessary actions after tracking
	if(extEff != NULL) extEff->finalizeEffects(bunch);
	//------------------------------------------------		
}

//--------------------------------------------------
// private methods of the RungeKuttaTracker class
//--------------------------------------------------
		
//return 0 if it is outside of both planes and 1 otherwise
int RungeKuttaTracker::isOutside(double* r){
	double resEntr = r[0]*plEntrV[0] + r[1]*plEntrV[1] + r[2]*plEntrV[2] + plEntrV[3];
	double resExit = r[0]*plExitV[0] + r[1]*plExitV[1] + r[2]*plExitV[2] + plExitV[3];
	if(resEntr < 0. && resExit < 0.){
		return 0;
	}
	return 1;
}

int RungeKuttaTracker::isOutside(double x, double y, double z){
	double resEntr = x*plEntrV[0] + y*plEntrV[1] + z*plEntrV[2] + plEntrV[3];
	double resExit = x*plExitV[0] + y*plExitV[1] + z*plExitV[2] + plExitV[3];
	if(resEntr < 0. && resExit < 0.){
		return 0;
	}
	return 1;
}

int RungeKuttaTracker::isAfterEntrance(double* r){
	double resEntr = r[0]*plEntrV[0] + r[1]*plEntrV[1] + r[2]*plEntrV[2] + plEntrV[3];
	if(resEntr < 0.){
		return 1;
	}
	return 0;	
}

int RungeKuttaTracker::isAfterEntrance(double x, double y, double z){
	double resEntr = x*plEntrV[0] + y*plEntrV[1] + z*plEntrV[2] + plEntrV[3];
	if(resEntr < 0.){
		return 1;
	}
	return 0;		
}

int RungeKuttaTracker::isBeforeExit(double* r){
	double resExit = r[0]*plExitV[0] + r[1]*plExitV[1] + r[2]*plExitV[2] + plExitV[3];
	if(resExit < 0.){
		return 1;
	}
	return 0;
}

int RungeKuttaTracker::isBeforeExit(double x, double y, double z){
	double resExit = x*plExitV[0] + y*plExitV[1] + z*plExitV[2] + plExitV[3];
	if(resExit < 0.){
		return 1;
	}
	return 0;	
}

void RungeKuttaTracker::rk4Step(double t, double t_st, BaseFieldSource* fieldSource){
	for(int i = 0; i < 6; i++){
		y_vct[i] = y_in_vct[i];
	}
	calculateRightSideODE(t,fieldSource);
	for(int i = 0; i < 6; i++){
		k1_vct[i] = t_st*ff_vct[i];
		y_vct[i] =  y_in_vct[i] + k1_vct[i]*0.5;
	}
	calculateRightSideODE(t+t_st*0.5,fieldSource);
	for(int i = 0; i < 6; i++){
		k2_vct[i] = t_st*ff_vct[i];
		y_vct[i] =  y_in_vct[i] + k2_vct[i]*0.5;
	}
	calculateRightSideODE(t+t_st*0.5,fieldSource);
	for(int i = 0; i < 6; i++){
		k3_vct[i] = t_st*ff_vct[i];
		y_vct[i] =  y_in_vct[i] + k3_vct[i];
	}
	calculateRightSideODE(t+t_st,fieldSource);
	for(int i = 0; i < 6; i++){
		k4_vct[i] = t_st*ff_vct[i];
		y_out_vct[i] =  y_in_vct[i] + (k1_vct[i] + 2*(k2_vct[i] + k3_vct[i]) + k4_vct[i])/6.0;
	}
}

void RungeKuttaTracker::calculateRightSideODE(double t, BaseFieldSource* fieldSource){
	fieldSource->getElectricMagneticField(y_vct[0],y_vct[1],y_vct[2],t,
		                                    e_vct[0],e_vct[1],e_vct[2],
																				b_vct[0],b_vct[1],b_vct[2]);
	//std::cerr<<"debug E=" << e_vct[0] <<" "<< e_vct[1] <<" "<< e_vct[2] <<std::endl;
	//std::cerr<<"debug B=" << b_vct[0] <<" "<< b_vct[1] <<" "<< b_vct[2] <<std::endl;
	double coef = mass2 + y_vct[3]*y_vct[3] + y_vct[4]*y_vct[4] + y_vct[5]*y_vct[5];
	coef = c_light/sqrt(coef);
	//std::cerr<<"debug coef=" << coef <<std::endl;
	ff_vct[0] = y_vct[3]*coef;
	ff_vct[1] = y_vct[4]*coef;
	ff_vct[2] = y_vct[5]*coef;
	coef *= charge;
	ff_vct[3] = c_light*(charge*e_vct[0] + coef*(y_vct[4]*b_vct[2] - y_vct[5]*b_vct[1]));
	ff_vct[4] = c_light*(charge*e_vct[1] + coef*(y_vct[5]*b_vct[0] - y_vct[3]*b_vct[2]));
	ff_vct[5] = c_light*(charge*e_vct[2] + coef*(y_vct[3]*b_vct[1] - y_vct[4]*b_vct[0]));
	ff_vct[3] = ff_vct[3]/1.0e+9;
	ff_vct[4] = ff_vct[4]/1.0e+9;
	ff_vct[5] = ff_vct[5]/1.0e+9;
}

