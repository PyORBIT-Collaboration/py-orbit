#include "RectangularApertureShape.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   RectangularApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
//   RectangularApertureShape is an implementation of BaseApertureShape class.
//
///////////////////////////////////////////////////////////////////////////

/** RectangularApertureShape constructor */
RectangularApertureShape::RectangularApertureShape(): BaseApertureShape()
{
		shapeName = "rectangular";
		typeName = "rectangular";
		x_half_size = 0.;
		y_half_size = 0.;
}

/** RectangularApertureShape decstructor */
RectangularApertureShape::~RectangularApertureShape()
{
}

/** Return 1 if the particular macro-particle is inside this shape */
int RectangularApertureShape::inside(Bunch* bunch, int count){
	
	if(x_half_size == 0. || y_half_size == 0.){
		return 0;
	}
	
	double** coord = bunch->coordArr();
	
	double x = fabs((coord[count][0] - x_center));
	double y = fabs((coord[count][2] - y_center));
	if(x <= x_half_size && y < y_half_size){
		return 1;
	}
	return 0;
}

/** Sets the half size in X-direction in [m] */ 
void RectangularApertureShape::setHalfX(double x_half_size)
{
	this->x_half_size = x_half_size;
}

/** Returns the half size in X-direction in [m]  */ 
double RectangularApertureShape::getHalfX()
{
	return x_half_size;
}

/** Sets the half size in Y-direction in [m] */ 
void RectangularApertureShape::setHalfY(double y_half_size)
{
	this->y_half_size = y_half_size;
}

/** Returns the half size in Y-direction in [m]  */ 
double RectangularApertureShape::getHalfY()
{
	return y_half_size;
}