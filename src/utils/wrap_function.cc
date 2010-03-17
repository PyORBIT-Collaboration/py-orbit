#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_function.hh"

#include <iostream>
#include <string>

#include "OU_Function.hh"

using namespace OrbitUtils;

namespace wrap_function{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	    Constructor for python class wrapping c++ Function instance.
      It never will be called directly.
	*/
	static PyObject* Function_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int Function_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  self->cpp_obj =  new Function();	  
	  ((Function*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
	/** It will add (x,y) pair to the Function instance */
  static PyObject* Function_add(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
	  double x,y;
		if(!PyArg_ParseTuple(	args,"dd:",&x,&y))
			error("pyFunction add(x,y) - parameters are needed");
		else {
			cpp_Function->add(x,y);
		}
		Py_INCREF(Py_None);
		return Py_None;
  }
	
 	/** It will return the number of (x,y) pairs in the Function instance */
  static PyObject* Function_getSize(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
	  int size = cpp_Function->getSize();
		return Py_BuildValue("i",size);
  }
	
 	/** It will return x for a particular index ind */
  static PyObject* Function_x(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(	args,"i:",&ind)){
			error("pyFunction x(index) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_Function->x(ind));
  }
	
 	/** It will return y for a particular index ind */
  static PyObject* Function_y(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(	args,"i:",&ind)){
			error("pyFunction y(index) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_Function->y(ind));
  }
	
 	/** It will return (x,y) for a particular index ind */
  static PyObject* Function_xy(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(	args,"i:",&ind)){
			error("pyFunction xy(index) - parameter is needed");
		}	
		return Py_BuildValue("(dd)",cpp_Function->x(ind),cpp_Function->y(ind));
  }	
	
  //-----------------------------------------------------
  //destructor for python Function class (__del__ method).
  //-----------------------------------------------------
  static void Function_del(pyORBIT_Object* self){
		//std::cerr<<"The Function __del__ has been called!"<<std::endl;
		delete ((Function*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python Function wrapper class
	// they will be vailable from python level
  static PyMethodDef FunctionClassMethods[] = {
		{ "add",				 Function_add,    	    METH_VARARGS,"Adds (x,y) to the Function container."},
		{ "getSize",		 Function_getSize,    	METH_VARARGS,"Returns the number of (x,y) in Function"},
		{ "x",				   Function_x,          	METH_VARARGS,"Returns x value for a point with a particular index"},
 		{ "y",				   Function_y,    	      METH_VARARGS,"Returns y value for a point with a particular index"},
 		{ "xy",				   Function_xy,    	      METH_VARARGS,"Returns (x,y) value for a point with a particular index"},
    {NULL}
  };
	
	// defenition of the memebers of the python Function wrapper class
	// they will be vailable from python level
	static PyMemberDef FunctionClassMembers [] = {
		{NULL}
	};
	
	//new python Function wrapper type definition
	static PyTypeObject pyORBIT_Function_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Function", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Function_del , /*tp_dealloc*/
		0, /*tp_print*/
		0, /*tp_getattr*/
		0, /*tp_setattr*/
		0, /*tp_compare*/
		0, /*tp_repr*/
		0, /*tp_as_number*/
		0, /*tp_as_sequence*/
		0, /*tp_as_mapping*/
		0, /*tp_hash */
		0, /*tp_call*/
		0, /*tp_str*/
		0, /*tp_getattro*/
		0, /*tp_setattro*/
		0, /*tp_as_buffer*/
		Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
		"The Function python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		FunctionClassMethods, /* tp_methods */
		FunctionClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Function_init, /* tp_init */
		0, /* tp_alloc */
		Function_new, /* tp_new */
	};	
	
	
	
	//--------------------------------------------------
	//Initialization function of the pyFunction class
	//--------------------------------------------------
  void initFunction(PyObject* module){
		if (PyType_Ready(&pyORBIT_Function_Type) < 0) return;
		Py_INCREF(&pyORBIT_Function_Type);
		PyModule_AddObject(module, "Function", (PyObject *)&pyORBIT_Function_Type);
		
	}

#ifdef __cplusplus
}
#endif


}
