//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    PyBaseFieldSource.cc
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for Python implementation of a field source. 
//    It should be sub-classed on Python level and implements 
//    getElectricField(x,y,z,t) and getMagneticField (x,y,z,t) methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//    ATTENTION: Using this class in real calculations is not wise! 
//               It is slow, because it delegates the field calculations
//               to the Python level. It is only for prototyping and 
//               debugging only.
//
///////////////////////////////////////////////////////////////////////////
#include "PyBaseFieldSource.hh"

#include "orbit_mpi.hh"
#include <iostream>

using namespace OrbitUtils;

PyBaseFieldSource::PyBaseFieldSource(PyObject* py_wrapperIn): CppPyWrapper(py_wrapperIn)
{
}

PyBaseFieldSource::~PyBaseFieldSource()
{
}

void PyBaseFieldSource::getElectricField(double x, double y, double z, double t, double& f_x, double& f_y, double& f_z)
{	  
	  PyObject* py_wrp = getPyWrapper();
  	PyObject* ef_tuple = PyObject_CallMethod(py_wrp,"getElectricField","dddd",x,y,z,t);
    if(!PyArg_ParseTuple(	ef_tuple,"ddd:electric_field",&f_x,&f_y,&f_z)){
      ORBIT_MPI_Finalize("PyBaseFieldSource - getElectricField(x,y,z,t0 method does not work!");
    }	
		//the references should be decreased because they were created as "new reference"
		Py_DECREF(ef_tuple);
}

void PyBaseFieldSource::getMagneticField(double x, double y, double z, double t, double& f_x, double& f_y, double& f_z)
{
	  PyObject* py_wrp = getPyWrapper();
  	PyObject* ef_tuple = PyObject_CallMethod(py_wrp,"getMagneticField","dddd",x,y,z,t);
    if(!PyArg_ParseTuple(	ef_tuple,"ddd:electric_field",&f_x,&f_y,&f_z)){
      ORBIT_MPI_Finalize("PyBaseFieldSource - getElectricField(x,y,z,t0 method does not work!");
    }	
		//the references should be decreased because they were created as "new reference"
		Py_DECREF(ef_tuple);
}

