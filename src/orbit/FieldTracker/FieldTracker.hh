
//The base class for Apertures. It defines the interface for Aperture
#ifndef FIELDTRACKER_H
#define FIELDTRACKER_H

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "Bunch.hh"
#include "Grid3D.hh"

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
			int nXgrid,int nYGrid,int nZGrid,
			double BxField3D,double ByField3D,double BzField3D,
			int zsymmetry);

	void ParseGrid3D(const string &fileName,
	                 const double &xmin, const double &xmax,
	                 const double &ymin, const double &ymax,
	                 const double &zmin, const double &zmax,
	                 const int &skipX,
	                 const int &skipY,
	                 const int &skipZ);

    
protected:
	double * XGrid;
	double * YGrid;
	double * ZGrid;
	Grid3D* BXGrid;
	Grid3D* BYGrid;
	Grid3D* BZGrid;
	Grid3D* BMagGrid;
	int nXGrid, nYGrid, nZGrid;


};

//end of FieldTracker_H ifdef
#endif
