#ifndef CONVEX_APERTURE_SHAPE_H
#define CONVEX_APERTURE_SHAPE_H

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Bunch.hh"
#include "BaseApertureShape.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   ConvexApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
///////////////////////////////////////////////////////////////////////////

/** 
  ConvexApertureShape is an implementation of BaseApertureShape class
  and a collection of (x,y) points describing convex shape. User should 
  define these points.	
*/
    
class ConvexApertureShape: public BaseApertureShape
{
public:
	
	/** ConvexApertureShape constructor */
	ConvexApertureShape();
	
	/** ConvexApertureShape decstructor */
	virtual ~ConvexApertureShape();
	
	/** Return 1 if the particular macro-particle is inside this shape */
	int inside(Bunch* bunch, int count);
	
	/** Adds the new aperture shape to the collection */ 
	void addPoint(double x, double y);
	
	/** Returns vector of x-coordinates of points */ 
	std::vector<double> getPointsX();
	
	/** Returns vector of y-coordinates of points */ 
	std::vector<double> getPointsY();
	
	/** Removes all (x,y) points */ 
	void removeAllPoints();
	
	/** Checks that we have a convex shape*/ 
	int checkAllPoints();
	
	/** Checks that a point with index i that it is at clockwise half-plain relative to 0->1 line */ 
	int checkOnePoint(int i0, int i1, int i);
	
	/** Checks that a point (x,y) that it is at clockwise half-plain relative to 0->1 line */ 
	int checkOnePoint(int i0, int i1, double x, double y);

protected:
	
	std::vector<double>  convexX;
	std::vector<double>  convexY;

};

//end of COMPOSITE_APERTURE_SHAPE_H ifdef
#endif