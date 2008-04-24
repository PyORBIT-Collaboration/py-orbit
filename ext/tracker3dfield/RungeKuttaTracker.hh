//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   RungeKuttaTracker.hh
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
#ifndef RUNGE_KUTTA_3D_TRACKER_H
#define RUNGE_KUTTA_3D_TRACKER_H

#include "Bunch.hh"
#include "BaseFieldSource.hh"
#include "ExternalEffects.hh"

namespace Tracker3DField{
	
	class RungeKuttaTracker
	{
		//--------------------------------------------------
		// public methods of the RungeKuttaTracker class
		//--------------------------------------------------
	public:
		
		RungeKuttaTracker(double length);
		~RungeKuttaTracker();
		
		//-----------------------------------
		//  public data members
		//-----------------------------------
		
		/** It sets the accuracy of the spatial resolution for tracking. */
		void setSpatialEps(double eps);

		/** It returns the accuracy of the spatial resolution for tracking. */
		double getSpatialEps();
		
		/** It sets the number of initial steps in tracking. */
		void setInitialStepsNumber(int n_init);
		
		/** It returns the number of initial steps in tracking. */
		int getInitialStepsNumber();
		
		/** It returns the number of steps during the synchronous particle tracking. */
		int getStepsNumber();
		
		/** It sets the approximate length of the region. */
		void setLength(double length);
		
		/** It returns the approximate length of the region. */
		double getLength();
		
		/** It sets (a*x+b*y+c*z+d=0) coefficients for the entrance plane. */
		void setEntrPlane(double a, double b, double c, double d);
		
		/** It sets (a*x+b*y+c*z+d=0) coefficients for the exit plane. */
		void setExitPlane(double a, double b, double c, double d);

		/** It returns (a*x+b*y+c*z+d=0) coefficients for the entrance plane. */
		void getEntrPlane(double& a, double& b, double& c, double& d);
		
		/** It returns (a*x+b*y+c*z+d=0) coefficients for the exit plane. */
		void getExitPlane(double& a, double& b, double& c, double& d);
		
		/** It tracks the bunch. The external effects instance could be NULL. */
		void trackBunch(Bunch* bunch, OrbitUtils::BaseFieldSource* fieldSource, ExternalEffects* extEff);
		
		
	private:
		//--------------------------------------------------
		// private methods of the RungeKuttaTracker class
		//--------------------------------------------------
		
		//return 0 if it is outside of both planes and 1 otherwise
		int isOutside(double* r);
		int isOutside(double x, double y, double z);
		
		int isAfterEntrance(double* r);
		int isAfterEntrance(double x, double y, double z);
		
		int isBeforeExit(double* r);
		int isBeforeExit(double x, double y, double z);
		
		//calculates ff_vct[6] - right side of ODE system
		//d(r)/d(t) = c*p/sqrt(p^2+m^2)
		//d(p)/d(t) = (c*E + c*q*[pxB]/sqrt(p^2+m^2))/(10^9)
		//p and m in GeV/c and GeV, c in [m/c], r in [m]
		//E in [V/m] and B in [T]
		void calculateRightSideODE(double t, OrbitUtils::BaseFieldSource* fieldSource);
		void rk4Step(double t, double t_st, OrbitUtils::BaseFieldSource* fieldSource);
		
		//-----------------------------------
		//  private data members
		//-----------------------------------
		
		//spatial eps in [m]
		double eps;
		
		//set initial number of steps
		int n_init;
		
		//steps during synch. part. tracking
		int n_steps;
		
		//time step
		double t_step;
		
		//approximate length of the tracker space [m]
		double length;
		
		
		//entry and exit planes coefficients a*x + b*y + c*z + d = 0
		//which is really n*(r-r0) = 0
		//the direction of a normal vector (a,b,c) vector is outside
		//the region 
		//(if (a,b,c)*r + d < 0 for both planes the point is inside region)
		double* plEntrV;
		double* plExitV;
		
		//vectors for ODE system of equations
		//y[6] r(position) and p(momentum) vectors
    double ff_vct[6];
		double y_init_vct[6];
		double y_in_vct[6];
		double y_vct[6];
		double y_out_vct[6];
		double y_final_vct[6];
		double e_vct[3];
		double b_vct[3];
		double c_light;
		double charge;
		double mass;
		double mass2;
		
		double k1_vct[6];
		double k2_vct[6];
		double k3_vct[6];
		double k4_vct[6];
	};
	
}; // end of Tracker3DField name-space


#endif
