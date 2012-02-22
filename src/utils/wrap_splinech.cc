#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_function.hh"

#include <iostream>
#include <string>

#include "OU_SplineCH.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_splinech{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	    Constructor for python class wrapping c++ SplineCH instance.
      It never will be called directly.
	*/
	static PyObject* SplineCH_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int SplineCH_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  self->cpp_obj =  new SplineCH();	  
	  ((SplineCH*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
	/** It will caluclate the SplineCH instance for the function */
  static PyObject* SplineCH_compile(PyObject *self, PyObject *args){
	  SplineCH* cpp_SplineCH = (SplineCH*)((pyORBIT_Object*) self)->cpp_obj;
	  PyObject* pyF;
		Function* f = NULL;
		int inf = -1;		
		if(!PyArg_ParseTuple(	args,"O:",&pyF))
			error("pySplineCH compile(Function F) - parameter is needed");
		else {
			PyObject* pyORBIT_Function_Type = getOrbitUtilsType("Function");
			if(!PyObject_IsInstance(pyF,pyORBIT_Function_Type)){
				error("pySplineCH - compile(Function F) - Function parameter is needed.");
			}			
			f= (Function*) ((pyORBIT_Object*) pyF)->cpp_obj;
			inf = cpp_SplineCH->compile(f);
		}		
		return Py_BuildValue("i",inf);
  }
	
 	/** It will return the number of (x,y) pairs in the SplineCH instance */
  static PyObject* SplineCH_getSize(PyObject *self, PyObject *args){
	  SplineCH* cpp_SplineCH = (SplineCH*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("i",cpp_SplineCH->getSize());
  }
	
 	/** It will return x for a particular index ind */
  static PyObject* SplineCH_x(PyObject *self, PyObject *args){
	  SplineCH* cpp_SplineCH = (SplineCH*)((pyORBIT_Object*) self)->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(	args,"i:",&ind)){
			error("pySplineCH x(index) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_SplineCH->x(ind));
  }
	
 	/** It will return y for a particular index ind */
  static PyObject* SplineCH_y(PyObject *self, PyObject *args){
	  SplineCH* cpp_SplineCH = (SplineCH*)((pyORBIT_Object*) self)->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(	args,"i:",&ind)){
			error("pySplineCH y(index) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_SplineCH->y(ind));
  }

 	/** It will return y for a specified x value */
  static PyObject* SplineCH_getY(PyObject *self, PyObject *args){
	  SplineCH* cpp_SplineCH = (SplineCH*)((pyORBIT_Object*) self)->cpp_obj;
		double val = 0.;
		if(!PyArg_ParseTuple(	args,"d:",&val)){
			error("pySplineCH getY(x) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_SplineCH->getY(val));
  }
	
 	/** It will return y for a specified x value */
  static PyObject* SplineCH_getYP(PyObject *self, PyObject *args){
	  SplineCH* cpp_SplineCH = (SplineCH*)((pyORBIT_Object*) self)->cpp_obj;
		double val = 0.;
		if(!PyArg_ParseTuple(	args,"d:",&val)){
			error("pySplineCH getYP(x) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_SplineCH->getYP(val));
  }
	
  //Prints SplineCH into the std::cout stream or file
  static PyObject* SplineCH_dump(PyObject *self, PyObject *args){
		SplineCH* cpp_SplineCH = (SplineCH*)((pyORBIT_Object*) self)->cpp_obj;
    //if nVars == 0 print into std::cout
    //if nVars == 1 print into the file
    int nVars = PyTuple_Size(args);
    const char* file_name = NULL;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        cpp_SplineCH->print(std::cout);
      }
      else{
        if(!PyArg_ParseTuple(	args,"s:dump",&file_name)){
          error("pySplineCH - dump(fileName) - a file name is needed");
        }
        cpp_SplineCH->print(file_name);
      }
    }
    else{
      error("pySplineCH. You should call dump() or dump(file_name)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }	
	

  //-----------------------------------------------------
  //destructor for python SplineCH class (__del__ method).
  //-----------------------------------------------------
  static void SplineCH_del(pyORBIT_Object* self){
		//std::cerr<<"The SplineCH __del__ has been called!"<<std::endl;
		delete ((SplineCH*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python SplineCH wrapper class
	// they will be vailable from python level
  static PyMethodDef SplineCHClassMethods[] = {
		{ "compile",			 SplineCH_compile,    	METH_VARARGS,"Creates the spline from Function."},
		{ "getSize",		 	 SplineCH_getSize,    	METH_VARARGS,"Returns the number of (x,y) in SplineCH"},
		{ "x",				     SplineCH_x,          	METH_VARARGS,"Returns x value for a point with a particular index"},
 		{ "y",				     SplineCH_y,    	      METH_VARARGS,"Returns y value for a point with a particular index"},
 		{ "getY",				 	 SplineCH_getY,    	    METH_VARARGS,"Returns y for a specified x value "},
 		{ "getYP",				 SplineCH_getYP,    	  METH_VARARGS,"Returns y' for a specified x value "},
 		{ "dump",				   SplineCH_dump,    	    METH_VARARGS,"Prints SplineCH into the std::cout stream or file"},
    {NULL}
  };
	
	// defenition of the memebers of the python SplineCH wrapper class
	// they will be vailable from python level
	static PyMemberDef SplineCHClassMembers [] = {
		{NULL}
	};
	
	//new python SplineCH wrapper type definition
	static PyTypeObject pyORBIT_SplineCH_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SplineCH", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SplineCH_del , /*tp_dealloc*/
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
		"The SplineCH python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SplineCHClassMethods, /* tp_methods */
		SplineCHClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SplineCH_init, /* tp_init */
		0, /* tp_alloc */
		SplineCH_new, /* tp_new */
	};	
	
	
	
	//--------------------------------------------------
	//Initialization function of the pySplineCH class
	//--------------------------------------------------
  void initSplineCH(PyObject* module){
		if (PyType_Ready(&pyORBIT_SplineCH_Type) < 0) return;
		Py_INCREF(&pyORBIT_SplineCH_Type);
		PyModule_AddObject(module, "SplineCH", (PyObject *)&pyORBIT_SplineCH_Type);
	}

#ifdef __cplusplus
}
#endif


}
