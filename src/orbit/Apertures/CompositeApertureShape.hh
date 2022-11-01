#ifndef COMPOSITE_APERTURE_SHAPE_H
#define COMPOSITE_APERTURE_SHAPE_H

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Bunch.hh"
#include "BaseApertureShape.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   CompositeApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
//   CompositeApertureShape is an implementation of BaseApertureShape class
//   to represent a logical union of several shapes dtored in collection. 
//   To get "isinside" method result 1 (Yes) this class will go through 
//   all BaseApertureShape class instances in the collection and 
//   should get 1 from at least one shape.  
//
///////////////////////////////////////////////////////////////////////////
    
class CompositeApertureShape: public BaseApertureShape
{
public:
	
	/** CompositeApertureShape constructor */
	CompositeApertureShape();
	
	/** CompositeApertureShape decstructor */
	virtual ~CompositeApertureShape();
	
	/** Return 1 if the particular macro-particle is inside this shape */
	int inside(Bunch* bunch, int count);
	
	/** Adds the new aperture shape to the collection */ 
	void addApertureShape(BaseApertureShape* apertureShape);
	
	
	/** Returns array of aperture shapes that are in this collection */ 
	std::vector<BaseApertureShape*> getApertureShape();

protected:
	
	 std::vector<BaseApertureShape*> apertureShapes;

};

//end of COMPOSITE_APERTURE_SHAPE_H ifdef
#endif
