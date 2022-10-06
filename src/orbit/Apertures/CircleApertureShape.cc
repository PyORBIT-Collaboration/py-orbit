#include "CircleApertureShape.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   CircleApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
//   CircleApertureShape is an implementation of BaseApertureShape class.
//
///////////////////////////////////////////////////////////////////////////

/** CircleApertureShape constructor */
CircleApertureShape::CircleApertureShape(): BaseApertureShape()
{
		shapeName = "circle";
		typeName = "circle";
		radius = 0.;
		radius2 = 0.;		
}

/** CircleApertureShape decstructor */
CircleApertureShape::~CircleApertureShape()
{
}

/** Return 1 if the particular macro-particle is inside this shape */
int CircleApertureShape::inside(Bunch* bunch, int count){
	
	double** coord = bunch->coordArr();
	
	double x = coord[count][0] - x_center;
	double y = coord[count][2] - y_center;
	double r2 = x*x + y*y;
	if(r2 <= radius2){
		return 1;
	}
	return 0;
}

/** Sets the radius of the circle in [m] */ 
void CircleApertureShape::setRadius(double radius)
{
	this->radius = radius;
	radius2 = radius*radius;
}

/** Returns the radius of the circle in [m] */ 
double CircleApertureShape::getRadius()
{
	return radius;
}
