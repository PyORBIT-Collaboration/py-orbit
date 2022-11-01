#include "PyBaseApertureShape.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   PyBaseApertureShape
//
// AUTHOR: 
//   Andrei Shishlo October 2022
//
//   PyBaseApertureShape is an implementation of BaseApertureShape class which 
//   will be exposed to Python level. It uses the Python methods to calculate 
//   the result of int PyBaseApertureShape::inside(...) method, so it is very 
//   slow and should be used only for testing, developing prototypes etc.
//
///////////////////////////////////////////////////////////////////////////

/** PyBaseApertureShape constructor */
PyBaseApertureShape::PyBaseApertureShape(): BaseApertureShape()
{
		shapeName = "python_class_shape";
		typeName = "python_class_shape";
}

/** PyBaseApertureShape decstructor */
PyBaseApertureShape::~PyBaseApertureShape()
{
}

/** Return 1 if the particular macro-particle is inside this shape */
int PyBaseApertureShape::inside(Bunch* bunch, int count){
	
	double** coord = bunch->coordArr();
	
	PyObject* py_wrp = getPyWrapper();
	PyObject* py_bunch = bunch->getPyWrapper();
	
	int res_isinside = 0;
	
	PyObject* py_res = PyObject_CallMethod(py_wrp,const_cast<char*>("inside"),const_cast<char*>("Oi"),py_bunch,count);
	
	res_isinside = (int) PyInt_AS_LONG(py_res);

	Py_DECREF(py_res);
	return res_isinside;
}

	

