#include "EllipseApertureShape.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   EllipseApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
//   EllipseApertureShape is an implementation of BaseApertureShape class.
//
///////////////////////////////////////////////////////////////////////////

/** EllipseApertureShape constructor */
EllipseApertureShape::EllipseApertureShape(): BaseApertureShape()
{
		shapeName = "ellipse";
		typeName = "ellipse";
		x_half_axis = 0.;
		y_half_axis = 0.;
}

/** EllipseApertureShape decstructor */
EllipseApertureShape::~EllipseApertureShape()
{
}

/** Return 1 if the particular macro-particle is inside this shape */
int EllipseApertureShape::inside(Bunch* bunch, int count){
	
	if(x_half_axis == 0. || y_half_axis == 0.){
		return 0;
	}
	
	double** coord = bunch->coordArr();
	
	double x = (coord[count][0] - x_center)/x_half_axis;
	double y = (coord[count][2] - y_center)/y_half_axis;
	double r2 = x*x + y*y;
	if(r2 <= 1.0){
		return 1;
	}
	return 0;
}

/** Sets the half axis of ellipse in X-direction in [m] */ 
void EllipseApertureShape::setHalfAxisX(double x_half_axis)
{
	this->x_half_axis = x_half_axis;
}

/** Returns the half axis of ellipse in X-direction in [m]  */ 
double EllipseApertureShape::getHalfAxisX()
{
	return x_half_axis;
}

/** Sets the half axis of ellipse in Y-direction in [m] */ 
void EllipseApertureShape::setHalfAxisY(double y_half_axis)
{
	this->y_half_axis = y_half_axis;
}

/** Returns the half axis of ellipse in Y-direction in [m]  */ 
double EllipseApertureShape::getHalfAxisY()
{
	return y_half_axis;
}