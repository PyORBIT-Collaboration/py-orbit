#include "CompositeApertureShape.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

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

/** CompositeApertureShape constructor */
CompositeApertureShape::CompositeApertureShape(): BaseApertureShape()
{
		shapeName = "composite";
		typeName = "composite";
}

/** CompositeApertureShape decstructor */
CompositeApertureShape::~CompositeApertureShape()
{	
	int n_shapes = apertureShapes.size();
	for(int ind = 0; ind < n_shapes; ind++){
		if(apertureShapes[ind]->getPyWrapper() != NULL){
			Py_XDECREF((PyObject*) apertureShapes[ind]->getPyWrapper());
		}
	}
}

/** Return 1 if the particular macro-particle is inside this shape */
int CompositeApertureShape::inside(Bunch* bunch, int count){
	int n_shapes = apertureShapes.size();
	if(n_shapes == 0){
		return 0;
	}	
	for(int ind = 0; ind < n_shapes; ind++){
		int isIn = apertureShapes[ind]->inside(bunch,count);
		if(isIn > 0){;
			return 1;
		}
	}
	return 0;
}

/** Adds the new aperture shape to the collection */ 
void CompositeApertureShape::addApertureShape(BaseApertureShape* apertureShape)
{
	apertureShapes.push_back(apertureShape);
	Py_INCREF((PyObject*) apertureShape->getPyWrapper());
}

/** Returns vector of pointers to aperture shapes that are in this collection */ 
std::vector<BaseApertureShape*> CompositeApertureShape::getApertureShape()
{
	return apertureShapes;
}
