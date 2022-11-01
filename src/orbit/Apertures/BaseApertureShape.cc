#include "BaseApertureShape.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   BaseApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
//   BaseApertureShape defines the interface for BaseApertureShape subclasses.
//   It decides if particular particle in the bunch is inside the aperture shape.
//
///////////////////////////////////////////////////////////////////////////

/** BaseApertureShape constructor */
BaseApertureShape::BaseApertureShape(): CppPyWrapper(NULL)
{
	shapeName = "no_shape";
	typeName = "no_type";
	x_center = 0.;
	y_center = 0.;
}

/** BaseApertureShape decstructor */
BaseApertureShape::~BaseApertureShape()
{
}

/** Return 1 if the particular macro-particle is inside this shape */
int BaseApertureShape::inside(Bunch* bunch, int count)
{
	return 0;
}

/** Sets the center of shape in X direction */
void BaseApertureShape::setCenterX(double x_center)
{
	this->x_center = x_center;
}

/** Returns the center of shape in X direction */
double BaseApertureShape::getCenterX()
{
	return x_center;
}

/** Sets the center of shape in Y direction */
void BaseApertureShape::setCenterY(double y_center)
{
	this->y_center = y_center;
}

/** Returns the center of shape in Y direction */
double BaseApertureShape::getCenterY()
{
	return y_center;
}


/** Returns the shape name */
string BaseApertureShape::getName()
{
	return shapeName;
}
	
/** Sets the shape name */
void BaseApertureShape::setName(string shapeNameIn)
{
	shapeName = shapeNameIn;
}

/** Returns the shape type name */
string BaseApertureShape::getTypeName()
{
	return typeName;
}
	

