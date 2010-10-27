#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_function.hh"

#include <iostream>
#include <string>

#include "OU_Function.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

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
	
 	/** It will return minimal x value in the Function */
  static PyObject* Function_getMinX(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_Function->getMinX());
  }	
	
 	/** It will return maximal x value in the Function */
  static PyObject* Function_getMaxX(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_Function->getMaxX());
  }	
	
 	/** It will return minimal y value in the Function */
  static PyObject* Function_getMinY(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_Function->getMinY());
  }	
	
 	/** It will return maximal y value in the Function */
  static PyObject* Function_getMaxY(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_Function->getMaxY());
  }	
		
 	/** It will remove all points in the Function */
  static PyObject* Function_clean(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		cpp_Function->clean();
	 	Py_INCREF(Py_None);
		return Py_None; 
  }	
	
 	/** It will free the memeory and will remove all points in the Function */
  static PyObject* Function_cleanMemory(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		cpp_Function->cleanMemory();
	 	Py_INCREF(Py_None);
		return Py_None; 
  }		
	
 	/** It will return y for a specified x value */
  static PyObject* Function_getY(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		double val = 0.;
		if(!PyArg_ParseTuple(	args,"d:",&val)){
			error("pyFunction getY(x) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_Function->getY(val));
  }
	
 	/** It will return x for a specified y value */
  static PyObject* Function_getX(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		double val = 0.;
		if(!PyArg_ParseTuple(	args,"d:",&val)){
			error("pyFunction getX(y) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_Function->getX(val));
  }
	
 	/** It will set the constant step flag to 1 if it is possible */
  static PyObject* Function_setConstStep(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		int inf = -1;
		if(!PyArg_ParseTuple(	args,"i:",&inf)){
			error("pyFunction setConstStep(inf) - parameter is needed");
		}	
		return Py_BuildValue("i",cpp_Function->setConstStep(inf));
  }	
	
 	/** It will return 1 if the step is const and 0 otherwise */
  static PyObject* Function_isStepConst(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("i",cpp_Function->isStepConst());
  }	
	
 	/** It will build the reverse Function if it is possible and return 1 or 0 */
  static PyObject* Function_setInverse(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		Function* rf = NULL;
	  PyObject* pyF;
		int inf = -1;
		if(!PyArg_ParseTuple(	args,"O:",&pyF))
			error("pyFunction setInverse(pyFunction F) - parameter is needed");
		else {
			PyObject* pyORBIT_Function_Type = getOrbitUtilsType("Function");
			if(!PyObject_IsInstance(pyF,pyORBIT_Function_Type)){
				error("pyFunction - setInverse(pyFunction F) - pyFunction parameter is needed.");
			}			
			rf= (Function*) ((pyORBIT_Object*) pyF)->cpp_obj;
			inf = cpp_Function->setInverse(rf);
		}
		return Py_BuildValue("i",inf);
  }	

  //Prints Function into the std::cout stream or file
  static PyObject* Function_dump(PyObject *self, PyObject *args){
		Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
    //if nVars == 0 print into std::cout
    //if nVars == 1 print into the file
    int nVars = PyTuple_Size(args);
    const char* file_name = NULL;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        cpp_Function->print(std::cout);
      }
      else{
        if(!PyArg_ParseTuple(	args,"s:dump",&file_name)){
          error("pyFunction - dump(fileName) - a file name is needed");
        }
        cpp_Function->print(file_name);
      }
    }
    else{
      error("pyFunction. You should call dump() or dump(file_name)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }	
	
 	/** It will return 1 if it is success and 0 otherwise */
  static PyObject* Function_normalize(PyObject *self, PyObject *args){
	  Function* cpp_Function = (Function*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("i",cpp_Function->normalize());
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
		{ "add",				 	 Function_add,    	    METH_VARARGS,"Adds (x,y) to the Function container."},
		{ "getSize",		 	 Function_getSize,    	METH_VARARGS,"Returns the number of (x,y) in Function"},
		{ "x",				     Function_x,          	METH_VARARGS,"Returns x value for a point with a particular index"},
 		{ "y",				     Function_y,    	      METH_VARARGS,"Returns y value for a point with a particular index"},
 		{ "xy",				     Function_xy,    	      METH_VARARGS,"Returns (x,y) value for a point with a particular index"},
 		{ "getMinX",		 	 Function_getMinX,    	METH_VARARGS,"Returns the minimal x value in the Function"},
 		{ "getMaxX",		 	 Function_getMaxX,    	METH_VARARGS,"Returns the maximal x value in the Function"},
 		{ "getMinY",		 	 Function_getMinY,    	METH_VARARGS,"Returns the minimal y value in the Function"},
 		{ "getMaxY",		 	 Function_getMaxY,    	METH_VARARGS,"Returns the maximal y value in the Function"},
 		{ "clean",			 	 Function_clean,    	  METH_VARARGS,"It will remove all points in the Function"},
 		{ "cleanMemory",	 Function_cleanMemory, METH_VARARGS,"It will free the memory and remove all points in the Function"},
 		{ "getY",				 	 Function_getY,    	    METH_VARARGS,"Returns y for a specified x value "},
 		{ "getX",				 	 Function_getX,    	    METH_VARARGS,"Returns x for a specified y value "},
 		{ "setConstStep", Function_setConstStep, METH_VARARGS,"It will set the constant step flag to 1 if it is possible"},
 		{ "isStepConst", 	 Function_isStepConst,  METH_VARARGS,"It will return 1 if the step is const and 0 otherwise"},
 		{ "setInverse",		 Function_setInverse,   METH_VARARGS,"It will build the reverse Function if it is possible and return 1 or 0"},
 		{ "dump",				   Function_dump,    	    METH_VARARGS,"Prints Function into the std::cout stream or file"},
 		{ "normalize",	   Function_normalize,   METH_VARARGS,"It will return 1 if it is success and 0 otherwise"},
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
