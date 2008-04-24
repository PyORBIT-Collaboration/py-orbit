//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppPyWrapper.cc
//
// CREATED
//    04/22/2003
//
// DESCRIPTION
//    The base class that provides capability to keep the reference to 
//    the Python class instance that is wrapping the c++ subclassing
//    this class.
//
///////////////////////////////////////////////////////////////////////////
#include "CppPyWrapper.hh"

using namespace OrbitUtils;

CppPyWrapper::CppPyWrapper(PyObject* py_wrapperIn)
{
	cpp_py_wrapper = py_wrapperIn;
}

CppPyWrapper::CppPyWrapper()
{
	cpp_py_wrapper = NULL;
}

CppPyWrapper::~CppPyWrapper()
{
}

void CppPyWrapper::setPyWrapper(PyObject* py_wrapperIn){
	cpp_py_wrapper = py_wrapperIn;
}

PyObject* CppPyWrapper::getPyWrapper(){
	return cpp_py_wrapper;
}
