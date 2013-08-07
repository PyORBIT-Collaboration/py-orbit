
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
    FieldTracker( const double &bx, const double &by,
 	       const double &ax, const double &ay,
 	       const double &ex, const double &epx,
 	       const double &l,
 	       const double &zi, const double &zf,
 	       const double &ds, const int &niters,
 	       const double &resid,
 	       const double &xrefi, const double &yrefi,
 	       const double &eulerai, const double &eulerbi,
 	       const double &eulergi, Bunch* b, string &filename);
    
	/** Routine for transfering particles through a aperture */
	void trackBunch(Bunch* b);

	void BGrid3D();

	void ParseGrid3D(const string &fileName,
	                 const double &xmin, const double &xmax,
	                 const double &ymin, const double &ymax,
	                 const double &zmin, const double &zmax,
	                 const int &skipX,
	                 const int &skipY,
	                 const int &skipZ);

	void initVars();
	void nodeCalculator(Bunch* b);
	void setPathVariable(int i);
    
protected:
	double * XGrid;
	double * YGrid;
	double * ZGrid;
	Grid3D* BXGrid;
	Grid3D* BYGrid;
	Grid3D* BZGrid;
	Grid3D* BMagGrid;
	int  nXGrid;
	int  nYGrid;
	int  nZGrid;
	double  BxField3D;
	double  ByField3D;
	double  BzField3D;
	double xField3D;
	double yField3D;
	double zField3D;

	double _length;
	double _xrefi;
	double _xreff;
	double _yrefi;
	double _yreff;
	double _eulerai;
	double _eulerbi;
	double _eulergi;
	double _euleraf;
	double _eulerbf;
	double _eulergf;
	double _zi;
	double _zf;
	double _ds;
	double _sref;
	double rmax;

	// NODE CALCULATOR STUFF
	double _betaX;
	double _betaY;
	double _alphaX;
	double _alphaY;
	double _etaX;
	double _etaPX;

	  double xFoilMin;
	  double xFoilMax;
	  double  yFoilMin;
	  double yFoilMax;
	  double zFoilMin;
	  double zFoilMax;
	  int doRefPath;
	//

	double _niters;

	double rMax;
	double foilAngle;

	double _resid;
	double sMax;
	int getPath;
	double BScale;
	int zsymmetry;





};

//end of FieldTracker_H ifdef
#endif
