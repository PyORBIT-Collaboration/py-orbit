#ifndef ELLIPSE_APERTURE_SHAPE_H
#define ELLIPSE_APERTURE_SHAPE_H

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Bunch.hh"
#include "BaseApertureShape.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   EllipseApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
///////////////////////////////////////////////////////////////////////////

/** 
	EllipseApertureShape is an implementation of BaseApertureShape class.
*/
    
class EllipseApertureShape: public BaseApertureShape
{
public:
	
	/** EllipseApertureShape constructor */
	EllipseApertureShape();
	
	/** EllipseApertureShape decstructor */
	virtual ~EllipseApertureShape();
	
	/** Return 1 if the particular macro-particle is inside this shape */
	int inside(Bunch* bunch, int count);
	
	/** Sets the half axis of ellipse in X-direction in [m] */ 
	void setHalfAxisX(double x_half_axis);
	
	/** Returns the half axis of ellipse in X-direction in [m]  */ 
	double getHalfAxisX();
	
	/** Sets the half axis of ellipse in Y-direction in [m] */ 
	void setHalfAxisY(double y_half_axis);
	
	/** Returns the half axis of ellipse in Y-direction in [m]  */ 
	double getHalfAxisY();
	
protected:
	
	double x_half_axis;
	double y_half_axis;

};

//end of ELLIPSE_APERTURE_SHAPE_H ifdef
#endif