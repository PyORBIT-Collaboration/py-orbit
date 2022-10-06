//The base class for BaseApertureShapes. It defines the interface for BaseApertureShape
#ifndef BASE_APERTURE_SHAPE_H
#define BASE_APERTURE_SHAPE_H

#include "Bunch.hh"
#include "BaseApertureShape.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   BaseApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
///////////////////////////////////////////////////////////////////////////

/** 
   The base class for BaseApertureShapes. 
   It defines the interface for BaseApertureShape subclasses.
*/
    
class BaseApertureShape: public OrbitUtils::CppPyWrapper
{
public:
	
	/** BaseApertureShape constructor */
	BaseApertureShape();
	
	/** BaseApertureShape decstructor */
	virtual ~BaseApertureShape();
	
	/** Return 1 if the particular macro-particle is inside this shape */
	virtual int inside(Bunch* bunch, int count);
	
	/** Sets the center of shape in X direction */
	void setCenterX(double x_center);
	
	/** Returns the center of shape in X direction */
	double getCenterX();
	
	/** Sets the center of shape in Y direction */
	void setCenterY(double y_center);
	
	/** Returns the center of shape in Y direction */
	double getCenterY();	

	/** Returns the shape name */
	string getName();
	
	/** Sets the shape name */
	void setName(string shapeName);
	
	/** Returns the shape type name */
	string getTypeName();
	
protected:
	
	string shapeName;
	
	string typeName;
	
	double x_center, y_center;

};

//end of BASE_APERTURE_SHAPE_H ifdef
#endif
