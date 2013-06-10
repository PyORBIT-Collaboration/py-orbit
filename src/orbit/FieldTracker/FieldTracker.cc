#include "FieldTracker.hh"



#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>



//////////////////////////////////////////////////////////////////////////
//
// Comments and Stuff
//
//
//
//
//////////////////////////////////////////////////////////////////////////

//Constructor

FieldTracker::FieldTracker(double a)
{
    
    cerr<<"Instantiating the 3D field track class \n";
    
}

void FieldTracker::trackBunch(Bunch* b)
{
    
	cerr<<"Ready to track in the 3D field!\n";
    
}

///////////////////////////////////////////////////////////////////////////
//
//
//
//
//
//
//
///////////////////////////////////////////////////////////////////////////

double * FieldTracker::BGrid3D(double xField3D,double yField3D,double zField3D,
		double XGrid[],double YGrid[],double ZGrid[],
		int nXGrid,int nYGrid,int nZGrid,
		double BxField3D,double ByField3D,double BzField3D,
		Grid3D BXGrid, Grid3D BYGrid,Grid3D BZGrid,
		int zsymmetry)
{
	  double xF, yF, zF;

	  xF = xField3D;
	  yF = yField3D;
	  zF = zField3D;

	  if((zsymmetry != 0) && (zF < ZGrid[0]))
	  {
	    zF = ZGrid[0] - zF + ZGrid[0];
	  }

	  int ig, jg, kg, igp, jgp, kgp;
	  double dx, dy, dz, xfrac, yfrac, zfrac, xfracc, yfracc, zfracc;

	  dx = (XGrid[nXGrid-1] - XGrid[0]) / double(nXGrid - 1);
	  dy = (YGrid[nYGrid-1] - YGrid[0]) / double(nYGrid - 1);
	  dz = (ZGrid[nZGrid-1] - ZGrid[0]) / double(nZGrid - 1);

	  ig = int((xF - XGrid[0]) / dx) + 1;
	  jg = int((yF - YGrid[0]) / dy) + 1;
	  kg = int((zF - ZGrid[0]) / dz) + 1;

	  if(ig < 0) ig = 0;
	  if(jg < 0) jg = 0;
	  if(kg < 0) kg = 0;
	  //could need to be changed to -2 do to array differences in C and SC
	  if(ig > nXGrid - 1) ig = nXGrid - 1;
	  if(jg > nYGrid - 1) jg = nYGrid - 1;
	  if(kg > nZGrid - 1) kg = nZGrid - 1;

	  igp = ig + 1;
	  jgp = jg + 1;
	  kgp = kg + 1;

	  xfrac = (xF - XGrid[ig]) / dx;
	  yfrac = (yF - YGrid[jg]) / dy;
	  zfrac = (zF - ZGrid[kg]) / dz;

	  xfracc = 1.0 - xfrac;
	  yfracc = 1.0 - yfrac;
	  zfracc = 1.0 - zfrac;

	  BxField3D =   BXGrid.getValueOnGrid(ig , jg , kg ) * xfracc * yfracc * zfracc +
	              + BXGrid.getValueOnGrid(ig , jg , kgp) * xfracc * yfracc * zfrac  +
	              + BXGrid.getValueOnGrid(ig , jgp, kg ) * xfracc * yfrac  * zfracc +
	              + BXGrid.getValueOnGrid(ig , jgp, kgp) * xfracc * yfrac  * zfrac  +
	              + BXGrid.getValueOnGrid(igp, jg , kg ) * xfrac  * yfracc * zfracc +
	              + BXGrid.getValueOnGrid(igp, jg , kgp) * xfrac  * yfracc * zfrac  +
	              + BXGrid.getValueOnGrid(igp, jgp, kg ) * xfrac  * yfrac  * zfracc +
	              + BXGrid.getValueOnGrid(igp, jgp, kgp) * xfrac  * yfrac  * zfrac;

	  ByField3D =   BYGrid.getValueOnGrid(ig , jg , kg ) * xfracc * yfracc * zfracc +
	              + BYGrid.getValueOnGrid(ig , jg , kgp) * xfracc * yfracc * zfrac  +
	              + BYGrid.getValueOnGrid(ig , jgp, kg ) * xfracc * yfrac  * zfracc +
	              + BYGrid.getValueOnGrid(ig , jgp, kgp) * xfracc * yfrac  * zfrac  +
	              + BYGrid.getValueOnGrid(igp, jg , kg ) * xfrac  * yfracc * zfracc +
	              + BYGrid.getValueOnGrid(igp, jg , kgp) * xfrac  * yfracc * zfrac  +
	              + BYGrid.getValueOnGrid(igp, jgp, kg ) * xfrac  * yfrac  * zfracc +
	              + BYGrid.getValueOnGrid(igp, jgp, kgp) * xfrac  * yfrac  * zfrac;

	  BzField3D =   BZGrid.getValueOnGrid(ig , jg , kg ) * xfracc * yfracc * zfracc +
	              + BZGrid.getValueOnGrid(ig , jg , kgp) * xfracc * yfracc * zfrac  +
	              + BZGrid.getValueOnGrid(ig , jgp, kg ) * xfracc * yfrac  * zfracc +
	              + BZGrid.getValueOnGrid(ig , jgp, kgp) * xfracc * yfrac  * zfrac  +
	              + BZGrid.getValueOnGrid(igp, jg , kg ) * xfrac  * yfracc * zfracc +
	              + BZGrid.getValueOnGrid(igp, jg , kgp) * xfrac  * yfracc * zfrac  +
	              + BZGrid.getValueOnGrid(igp, jgp, kg ) * xfrac  * yfrac  * zfracc +
	              + BZGrid.getValueOnGrid(igp, jgp, kgp) * xfrac  * yfrac  * zfrac;
	  if((zsymmetry != 0) && (zF < ZGrid[1]))
	  {
	    BzField3D = -BzField3D;
	  }
	  double FieldArr[3];
	  FieldArr[0] = BxField3D;
	  FieldArr[1] = ByField3D;
	  FieldArr[2] = BzField3D;
	  return FieldArr;


}

