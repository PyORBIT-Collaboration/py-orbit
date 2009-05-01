//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    PyExternalEffects.cc
//
// CREATED
//    05/14/2008
//
// DESCRIPTION
//    The base class for Python implementation of a external effects 
//    during the transport of particles through the external field. 
//    It should be sub-classed on Python level and implement
//    setupEffects(Bunch* bunch)
//    finalizeEffects(Bunch* bunch)
//    applyEffects(Bunch* bunch, int index, 
//	                            double* y_in_vct, double* y_out_vct, 
//														  double t, double t_step, 
//														  OrbitUtils::BaseFieldSource* fieldSource)
//    methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//    ATTENTION: Using this class in real calculations is not wise! 
//               It is slow, because it delegates the field calculations
//               to the Python level. It is for prototyping and 
//               debugging only.
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include <iostream>

#include "PyExternalEffects.hh"
#include "RungeKuttaTracker.hh"

using namespace TrackerRK4;
using namespace OrbitUtils;

PyExternalEffects::PyExternalEffects(PyObject* py_wrapperIn)
{
	setName("PyExternalEffects");
	setPyWrapper(py_wrapperIn);
}

PyExternalEffects::~PyExternalEffects()
{
}

void PyExternalEffects::setupEffects(Bunch* bunch){
	PyObject* py_wrp = getPyWrapper();
	PyObject* py_bunch = bunch->getPyWrapper();
	PyObject* res_tuple = PyObject_CallMethod(py_wrp,"setupEffects","O",py_bunch);
	Py_DECREF(res_tuple);
}

void PyExternalEffects::memorizeInitParams(Bunch* bunch){
	PyObject* py_wrp = getPyWrapper();
	PyObject* py_bunch = bunch->getPyWrapper();
	PyObject* res_tuple = PyObject_CallMethod(py_wrp,"memorizeInitParams","O",py_bunch);
	Py_DECREF(res_tuple);
}
		
	
void PyExternalEffects::finalizeEffects(Bunch* bunch){
	PyObject* py_wrp = getPyWrapper();
	PyObject* py_bunch = bunch->getPyWrapper();
	PyObject* res_tuple = PyObject_CallMethod(py_wrp,"finalizeEffects","O",py_bunch);
	Py_DECREF(res_tuple);	
}

void PyExternalEffects::applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)
{
	PyObject* py_wrp = getPyWrapper();
	PyObject* py_bunch = bunch->getPyWrapper();
	PyObject* py_field = fieldSource->getPyWrapper();
	PyObject* py_tracker = tracker->getPyWrapper();
	PyObject* pyInVct = Py_BuildValue("(dddddd)",y_in_vct[0],y_in_vct[1],y_in_vct[2],
		                                           y_in_vct[3],y_in_vct[4],y_in_vct[5]);
	PyObject* pyOutVct = Py_BuildValue("(dddddd)",y_out_vct[0],y_out_vct[1],y_out_vct[2],
		                                           y_out_vct[3],y_out_vct[4],y_out_vct[5]);
	PyObject* res_tuple = PyObject_CallMethod(py_wrp,"applyEffects","OiOOddOO",
		py_bunch,
		index,
		pyInVct,pyOutVct,
		t,t_step,
		py_field,
		py_tracker);
	Py_DECREF(res_tuple);
	Py_DECREF(pyInVct);
	Py_DECREF(pyOutVct);
}

