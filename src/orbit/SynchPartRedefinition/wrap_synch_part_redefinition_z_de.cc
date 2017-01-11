#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_synch_part_redefinition_z_de.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "SynchPartRedefinitionZdE.hh"

namespace wrap_synch_part_redefinition{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	    Constructor for python class wrapping c++ SynchPartRedefinitionZdE instance.
      It never will be called directly.
	*/
	static PyObject* SynchPartRedefinitionZdE_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int SynchPartRedefinitionZdE_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj =  new SynchPartRedefinitionZdE();
	  ((SynchPartRedefinitionZdE*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
 /** Performs the calculation of the z and dE averages of the bunch */
  static PyObject* SynchPartRedefinitionZdE_analyzeBunch(PyObject *self, PyObject *args){
	  SynchPartRedefinitionZdE* cpp_SynchPartRedefinitionZdE = (SynchPartRedefinitionZdE*)((pyORBIT_Object*) self)->cpp_obj;
		PyObject* pyBunch;
		if(!PyArg_ParseTuple(args,"O:analyzeBunch",&pyBunch)){
			ORBIT_MPI_Finalize("SynchPartRedefinitionZdE - analyzeBunch(Bunch* bunch) - parameter are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("SynchPartRedefinitionZdE - analyzeBunch(Bunch* bunch) - method needs a Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		cpp_SynchPartRedefinitionZdE->analyzeBunch(cpp_bunch);
		Py_INCREF(Py_None);
		return Py_None;
  }
	
  /** Transforms the synch particle parameters and the coordinates of the particles */
  static PyObject* SynchPartRedefinitionZdE_transformBunch(PyObject *self, PyObject *args){
	  SynchPartRedefinitionZdE* cpp_SynchPartRedefinitionZdE = (SynchPartRedefinitionZdE*)((pyORBIT_Object*) self)->cpp_obj;
		PyObject* pyBunch;
		if(!PyArg_ParseTuple(args,"O:transformBunch",&pyBunch)){
			ORBIT_MPI_Finalize("SynchPartRedefinitionZdE - transformBunch(Bunch* bunch) - parameter are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("SynchPartRedefinitionZdE - transformBunch(Bunch* bunch) - method needs a Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		cpp_SynchPartRedefinitionZdE->transformBunch(cpp_bunch);
		Py_INCREF(Py_None);
		return Py_None;
  } 
  
	/** Returns the average z postion */
	static PyObject* SynchPartRedefinitionZdE_getAvg_Z(PyObject *self, PyObject *args){
		SynchPartRedefinitionZdE* cpp_SynchPartRedefinitionZdE = (SynchPartRedefinitionZdE*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_SynchPartRedefinitionZdE->getAvg_Z());
	}
	
	/** Returns the average dE value */
	static PyObject* SynchPartRedefinitionZdE_getAvg_dE(PyObject *self, PyObject *args){
		SynchPartRedefinitionZdE* cpp_SynchPartRedefinitionZdE = (SynchPartRedefinitionZdE*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_SynchPartRedefinitionZdE->getAvg_dE());
	}	
		
  //-----------------------------------------------------
  //destructor for python SynchPartRedefinitionZdE class (__del__ method).
  //-----------------------------------------------------
  static void SynchPartRedefinitionZdE_del(pyORBIT_Object* self){
		//std::cerr<<"The SynchPartRedefinitionZdE __del__ has been called!"<<std::endl;
		delete ((SynchPartRedefinitionZdE*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python SynchPartRedefinitionZdE wrapper class
	// they will be vailable from python level
  static PyMethodDef SynchPartRedefinitionZdEClassMethods[] = {
		{ "analyzeBunch",	SynchPartRedefinitionZdE_analyzeBunch,		METH_VARARGS,"Calculates of the z and dE averages of the bunch."},
		{ "transformBunch",SynchPartRedefinitionZdE_transformBunch,	METH_VARARGS,"Transforms the synch part. params and particles' coordinates."},
 		{ "getAvg_Z",			SynchPartRedefinitionZdE_getAvg_Z,    		METH_VARARGS,"Returns the average z postion."},
 		{ "getAvg_dE",		SynchPartRedefinitionZdE_getAvg_dE,    		METH_VARARGS,"Returns the average dE value."},			
		{NULL}
  };
	
	// defenition of the memebers of the python SynchPartRedefinitionZdE wrapper class
	// they will be vailable from python level
	static PyMemberDef SynchPartRedefinitionZdEClassMembers [] = {
		{NULL}
	};
	
	//new python SynchPartRedefinitionZdE wrapper type definition
	static PyTypeObject pyORBIT_SynchPartRedefinitionZdE_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SynchPartRedefinitionZdE", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SynchPartRedefinitionZdE_del , /*tp_dealloc*/
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
		"The SynchPartRedefinitionZdE python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SynchPartRedefinitionZdEClassMethods, /* tp_methods */
		SynchPartRedefinitionZdEClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SynchPartRedefinitionZdE_init, /* tp_init */
		0, /* tp_alloc */
		SynchPartRedefinitionZdE_new, /* tp_new */
	};	
	
	//--------------------------------------------------
	//Initialization SynchPartRedefinitionZdE of the pySynchPartRedefinitionZdE class
	//--------------------------------------------------
  void initsynchpartredefinition(PyObject* module){
		if (PyType_Ready(&pyORBIT_SynchPartRedefinitionZdE_Type) < 0) return;
		Py_INCREF(&pyORBIT_SynchPartRedefinitionZdE_Type);
		PyModule_AddObject(module, "SynchPartRedefinitionZdE", (PyObject *)&pyORBIT_SynchPartRedefinitionZdE_Type);
	}

#ifdef __cplusplus
}
#endif


}
