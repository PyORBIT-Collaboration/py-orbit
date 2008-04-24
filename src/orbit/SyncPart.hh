//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    SyncPart.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    12/14/2005
//
// DESCRIPTION
//    Specification and inline functions for a  synchronous macro-particle.
//    The synchronous particle bears the coordinate system with z-axis 
//    directing along pz vector, and x-axis perpendicular to z. The y-axis
//    ort direction is defined as ny = [nz*nx].
//
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <string>
#include <set>
#include <map>
#include <vector>

using namespace std;

#ifndef SYNC_PARTICLE_H
#define SYNC_PARTICLE_H
///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    SyncPart
//
///////////////////////////////////////////////////////////////////////////
class Bunch;

class  SyncPart
{
	public:
  //--------------------------------------
  //the public methods of the SyncPart class
  //--------------------------------------
  virtual ~SyncPart();

	void setPyWrapper(PyObject* py_wrapper_In);
	PyObject* getPyWrapper();

	/**
	  Kinetic energy in GeV
	*/
	double getEnergy();

	/**
	  Rest mass in GeV
	*/
	double getMass();

	/**
	  time in seconds
	*/
	void setTime(double time);
	double getTime();

	/**
	   positions in meters
	*/
	void setXYZ(const double* xyz);
	void setXYZ(double x, double y, double z);
	void setX(double x);
	void setY(double y);
	void setZ(double z);
	double getX();
	double getY();
	double getZ();

	/**
	   momentum in GeV/c
	*/
	void setPXYZ(const double* pxyz);
	void setPXYZ(double px, double py, double pz);
	void setPX(double px);
	void setPY(double py);
	void setPZ(double pz);
	void setMomentum(double p);
	double getPX();
	double getPY();
	double getPZ();
	double getMomentum();

	/**
	  normal unit vector of x-axis
	*/
	void setNormalX(double x, double y, double z);
	double getNormalXX();
	double getNormalXY();
	double getNormalXZ();
	
	double getBeta();
	double getGamma();

	double momentumToEnergy(double p);
	double energyToMomentum(double ek);

protected:

  //Initializes the different data for SyncPartes.
  void init();

private:

  void updateKinematics();

  //---------------------------------------
  //the private methods of the SyncPart class
  //---------------------------------------

  friend class Bunch;

  SyncPart(Bunch* bunch);

	//initilaze the sync. particle from file
	void readSyncPart(const char* fileName);

	//print data abount synch. particle into the stream
	void print(std::ostream& Out);

  Bunch* bunch;
	PyObject* py_wrapper;

	//--------------------------------------------
	//parameters
	//--------------------------------------------

	//energy in GeV and momentum in GeV/c
	double energy;
	double p_abs;
	double beta;
	double gamma;

	//time in sec
	double time;

	//position [m] and momentum [GeV/c]
	double xyz[3];
	double pxyz[3];
	
	//vector of local x-axis
	double x_normal[3];
	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
