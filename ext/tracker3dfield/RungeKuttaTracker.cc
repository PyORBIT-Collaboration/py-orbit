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

using namespace Tracker3DField;
using namespace OrbitUtils;
	
#include "OrbitConst.hh"
	
RungeKuttaTracker::RungeKuttaTracker(double lengthIn){
	
	length = lengthIn;
	
	//spatial eps in [m]
	eps = 0.00001;
			
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


/** It tracks the bunch. The external effects instance could be NULL. */
void RungeKuttaTracker::trackBunch(Bunch* bunch, BaseFieldSource* fieldSource, ExternalEffects* extEff){
	int nMax = 5*n_steps;
	c_light = OrbitConst::c;
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
		ORBIT_MPI_Finalize("RungeKuttaTracker::trackBunch - syncPart cannot cross the entrance plane (0).");
	}
	cross_t = -(plEntrV[0]*syncPart->getX() + 
		          plEntrV[1]*syncPart->getY() + 
							plEntrV[2]*syncPart->getZ() + plEntrV[3])/cross_t;
	double pSyncPart = syncPart->getMomentum();
	y_init_vct[0] = cross_t*syncPart->getPX() + syncPart->getX() + eps*syncPart->getPX()/pSyncPart;
	y_init_vct[1] = cross_t*syncPart->getPY() + syncPart->getY() + eps*syncPart->getPY()/pSyncPart;
	y_init_vct[2] = cross_t*syncPart->getPZ() + syncPart->getZ() + eps*syncPart->getPZ()/pSyncPart;
	y_init_vct[3] = syncPart->getPX();
	y_init_vct[4] = syncPart->getPY();
	y_init_vct[5] = syncPart->getPZ();
	//this is a time when synch particle crosses the entrance plane
	double t_sync_start = (cross_t*pSyncPart + eps)/(c_light*syncPart->getBeta());
	//std::cerr<<"debug start time="<< t_sync_start <<std::endl;
	if(isOutside(y_in_vct)){
		ORBIT_MPI_Finalize("RungeKuttaTracker::trackBunch - syncPart cannot cross the entrance plane (1).");
	}
	double y_sync_final_vct[6];
	//step size definition from synch. particle
	double t_st = t_step;
	//std::cerr<<"debug length [m]="<< length <<std::endl;
	//std::cerr<<"debug c_light="<< c_light <<std::endl;
	//std::cerr<<"debug beta="<< syncPart->getBeta() <<std::endl;
	//std::cerr<<"debug time step [sec]="<< t_st <<std::endl;
	for(int iN = 1; iN < 16; iN++){
		nMax = nMax*2;
		for(int i = 0; i < 6; i++){
			y_in_vct[i] = y_init_vct[i];
			y_out_vct[i] = y_init_vct[i];
		}
		int step_count = 0;
		double t = 0.;
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
		//std::cerr<<"debug out r=" << y_out_vct[0] <<" "<< y_out_vct[1] <<" "<< y_out_vct[2] <<std::endl;
		//std::cerr<<"debug Exit Pl=" << plExitV[0] <<" "<< plExitV[1] <<" "<< plExitV[2]<<" "<< plExitV[3] <<std::endl;
		cross_t = (y_out_vct[0] - y_in_vct[0])*plExitV[0] +
		(y_out_vct[1] - y_in_vct[1])*plExitV[1] +
		(y_out_vct[2] - y_in_vct[2])*plExitV[2];
		cross_t = -(plExitV[0]*y_in_vct[0] + 
		            plExitV[1]*y_in_vct[1] + 
			          plExitV[2]*y_in_vct[2] + plExitV[3])/cross_t;	
		for(int i = 0; i < 6; i++){
			y_final_vct[i] = cross_t* (y_out_vct[i] - y_in_vct[i])+ y_in_vct[i];
		}
		//std::cerr<<"debug final r=" << y_final_vct[0] <<" "<< y_final_vct[1] <<" "<< y_final_vct[2] <<std::endl;
		if(iN > 1){
			double diff = 0.;
			for(int i = 0; i < 3; i++){
				diff = diff + (y_final_vct[i]-y_sync_final_vct[i])*(y_final_vct[i]-y_sync_final_vct[i]);
			}
			diff = sqrt(diff);
			//remember the last one
			for(int i = 0; i < 6; i++){
				y_sync_final_vct[i] = y_final_vct[i];
			}
			t_st = t_st/2;
			if(diff < eps){
				//that is enough!
			  t_step = t_st;
				n_steps = step_count;
				nMax = 10*n_steps;
				break;
			}
		}
	}
	//integration parameters are ready
  std::cerr<<"debug time step [sec]="<< t_step <<std::endl;
  std::cerr<<"debug number of steps="<< n_steps <<std::endl;
	double x = 2. ,y = 2., z = 2. , t = 3.;
	double fx,fy,fz;
	fieldSource->getElectricField(x,y,z,t,fx,fy,fz);
	std::cerr<<"debug Electr. field ="<< fx <<" "<< fy <<" "<< fz <<std::endl;
  fieldSource->getMagneticField(x,y,z,t,fx,fy,fz);
	std::cerr<<"debug Magnet. field ="<< fx <<" "<< fy <<" "<< fz <<std::endl;
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
	fieldSource->getElectricField(y_vct[0],y_vct[1],y_vct[2],t,e_vct[0],e_vct[1],e_vct[2]);
	fieldSource->getMagneticField(y_vct[0],y_vct[1],y_vct[2],t,b_vct[0],b_vct[1],b_vct[2]);
	//std::cerr<<"debug E=" << e_vct[0] <<" "<< e_vct[1] <<" "<< e_vct[2] <<std::endl;
	//std::cerr<<"debug B=" << b_vct[0] <<" "<< b_vct[1] <<" "<< b_vct[2] <<std::endl;
	double coef = mass2 + y_vct[3]*y_vct[3] + y_vct[4]*y_vct[4] + y_vct[5]*y_vct[5];
	coef = c_light/sqrt(coef);
	//std::cerr<<"debug coef=" << coef <<std::endl;
	ff_vct[0] = y_vct[3]*coef;
	ff_vct[1] = y_vct[4]*coef;
	ff_vct[2] = y_vct[5]*coef;
	coef *= charge;
	ff_vct[3] = c_light*(e_vct[0] + coef*(y_vct[4]*b_vct[2] - y_vct[5]*b_vct[1]));
	ff_vct[4] = c_light*(e_vct[1] + coef*(y_vct[5]*b_vct[0] - y_vct[3]*b_vct[2]));
	ff_vct[5] = c_light*(e_vct[2] + coef*(y_vct[3]*b_vct[1] - y_vct[4]*b_vct[0]));
	ff_vct[3] = ff_vct[3]/1.0e+9;
	ff_vct[4] = ff_vct[4]/1.0e+9;
	ff_vct[5] = ff_vct[5]/1.0e+9;
}

