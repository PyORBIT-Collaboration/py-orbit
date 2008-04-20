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

		
RungeKuttaTracker::RungeKuttaTracker(double length){
	
	//spatial eps in [m]
	eps = 0.00003;
			
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
}

//--------------------------------------------------
// private methods of the RungeKuttaTracker class
//--------------------------------------------------
		
//return 0 if it is outside of both planes and 1 otherwise
int RungeKuttaTracker::isOutside(double* r){
	double resEntr = r[0]*plEntrV[0] + r[1]*plEntrV[1] + r[2]*plEntrV[2] + plEntrV[3];
	double resExit = r[0]*plExitV[0] + r[1]*plExitV[1] + r[2]*plExitV[2] + plExitV[3];
	if(resEntr < 0. && resExit < 0.){
		return 1;
	}
	return 0;
}

int RungeKuttaTracker::isOutside(double x, double y, double z){
	double resEntr = x*plEntrV[0] + y*plEntrV[1] + z*plEntrV[2] + plEntrV[3];
	double resExit = x*plExitV[0] + y*plExitV[1] + z*plExitV[2] + plExitV[3];
	if(resEntr < 0. && resExit < 0.){
		return 1;
	}
	return 0;
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


