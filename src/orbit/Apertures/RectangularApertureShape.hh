#ifndef RECTANGULAR_APERTURE_SHAPE_H
#define RECTANGULAR_APERTURE_SHAPE_H

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Bunch.hh"
#include "BaseApertureShape.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   RectangularApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
///////////////////////////////////////////////////////////////////////////

/** 
	RectangularApertureShape is an implementation of BaseApertureShape class.
*/
    
class RectangularApertureShape: public BaseApertureShape
{
public:
	
	/** RectangularApertureShape constructor */
	RectangularApertureShape();
	
	/** RectangularApertureShape decstructor */
	virtual ~RectangularApertureShape();
	
	/** Return 1 if the particular macro-particle is inside this shape */
	int inside(Bunch* bunch, int count);
	
	/** Sets the half size in X-direction in [m] */ 
	void setHalfX(double x_half_size);
	
	/** Returns the half size in X-direction in [m]  */ 
	double getHalfX();
	
	/** Sets the half size in Y-direction in [m] */ 
	void setHalfY(double y_half_size);
	
	/** Returns the half size in Y-direction in [m]  */ 
	double getHalfY();
	
protected:
	
	double x_half_size;
	double y_half_size;

};

//end of RECTANGULAR_APERTURE_SHAPE_H ifdef
#endif