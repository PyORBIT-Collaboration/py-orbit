//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppPyWrapper.hh
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
#ifndef CPP_PY_WRAPPER_H
#define CPP_PY_WRAPPER_H

#include "Python.h"

namespace OrbitUtils{
	
  /**    
	  The base class that provides capability to keep the reference to 
    the Python class instance that is wrapping the c++ subclassing
    this class.	
	*/
	
	class  CppPyWrapper
	{
		public:
		
			/** Constructor of the CppPyWrapper class with reference to Python class instance. */
			CppPyWrapper(PyObject* py_wrapperIn);
			
			/** Constructor of the CppPyWrapper class with NULL reference to Python class. */
			CppPyWrapper();
			
			/** Destrictor. It is empty. */
			~CppPyWrapper();
		
			/** Sets the reference to Python class instance. */
			void setPyWrapper(PyObject* py_wrapperIn);
			
			/** Returns the reference to Python class instance. */
			PyObject* getPyWrapper();
			
			private:
				
				PyObject* cpp_py_wrapper;

	};
};

#endif
