
//The base class for Apertures. It defines the interface for Aperture
#ifndef FIELDTRACKER_H
#define FIELDTRACKER_H

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "Bunch.hh"
#include <Grid3D.hh>

using namespace std;

/**
 The FieldTracker class is used to define how a particle/bunch  propogates through an 
 3-dimensional magnetic field
 */

class FieldTracker: public OrbitUtils::CppPyWrapper
{
public:
	
	/** FieldTracker */
    FieldTracker(double a);
    
	/** Routine for transfering particles through a aperture */
	void trackBunch(Bunch* b);

	double * BGrid3D(double xField3D,double yField3D,double zField3D,
			double XGrid[],double YGrid[],double ZGrid[],
			int nXgrid,int nYGrid,int nZGrid,
			double BxField3D,double ByField3D,double BzField3D,
			Grid3D BXGrid, Grid3D BYGrid,Grid3D BZGrid,
			int zsymmetry);

    
protected:

};

//end of FieldTracker_H ifdef
#endif
