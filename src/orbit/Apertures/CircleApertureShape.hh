#ifndef CIRCULAR_APERTURE_SHAPE_H
#define CIRCULAR_APERTURE_SHAPE_H

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Bunch.hh"
#include "BaseApertureShape.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   CircleApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
///////////////////////////////////////////////////////////////////////////

/** 
	CircleApertureShape is an implementation of BaseApertureShape class.
*/
    
class CircleApertureShape: public BaseApertureShape
{
public:
	
	/** CircleApertureShape constructor */
	CircleApertureShape();
	
	/** CircleApertureShape decstructor */
	virtual ~CircleApertureShape();
	
	/** Return 1 if the particular macro-particle is inside this shape */
	int inside(Bunch* bunch, int count);
	
	/** Sets the radius of the circle in [m] */ 
	void setRadius(double radius);
	
	/** Returns the radius of the circle in [m] */ 
	double getRadius();
	
protected:
	
	double radius;
	double radius2;

};

//end of CIRCULAR_APERTURE_SHAPE_H ifdef
#endif