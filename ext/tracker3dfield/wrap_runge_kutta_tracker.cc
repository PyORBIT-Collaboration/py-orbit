#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_runge_kutta_tracker.hh"

#include <iostream>

#include "RungeKuttaTracker.hh"

using namespace Tracker3DField;
using namespace OrbitUtils;

namespace wrap_tracker3dfield{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python RungeKuttaTracker class definition
	//---------------------------------------------------------

	//constructor for python class wrapping RungeKuttaTracker instance
	//It never will be called directly
	static PyObject* RungeKuttaTracker_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The RungeKuttaTracker new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  RungeKuttaTracker class
  //this is implementation of the __init__ method
  static int RungeKuttaTracker_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		double length = 0.;
		if(!PyArg_ParseTuple(args,"d:__init__",&length)){
			error("PyRungeKuttaTracker - RungeKuttaTracker(length[m]) - constructor needs a parameter.");
		}
		self->cpp_obj = new RungeKuttaTracker(length);
		//std::cerr<<"The RungeKuttaTracker __init__ has been called!"<<std::endl;
		return 0;
  }

	//  track(bunch,field_source,exteranl_effects) - track bunch
  static PyObject* RungeKuttaTracker_track(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		PyObject *pyBunch = NULL;
		PyObject *pyFieldSource = NULL;
		PyObject *pyExtEffects = NULL;
		if(!PyArg_ParseTuple(args,"OO|O:track",&pyBunch,&pyFieldSource,&pyExtEffects)){
			error("PyRungeKuttaTracker - track(bunch,field_source[,exteranl_effects]) - parameters are needed.");
		}
		Bunch* bunch = (Bunch*) ((pyORBIT_Object*) pyBunch)->cpp_obj;
		BaseFieldSource* fs = (BaseFieldSource*) ((pyORBIT_Object*) pyFieldSource)->cpp_obj;
		ExternalEffects* extEf = NULL;
		if(pyExtEffects != NULL){
			extEf = (ExternalEffects*) ((pyORBIT_Object*) pyExtEffects)->cpp_obj;
		}
		cpp_RungeKuttaTracker->trackBunch(bunch,fs,extEf);
		Py_INCREF(Py_None);
    return Py_None;	
  }	
	
  //-----------------------------------------------------
  //destructor for python RungeKuttaTracker class (__del__ method).
  //-----------------------------------------------------
  static void RungeKuttaTracker_del(pyORBIT_Object* self){
		//std::cerr<<"The RungeKuttaTracker __del__ has been called!"<<std::endl;
		delete ((RungeKuttaTracker*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python RungeKuttaTracker wrapper class
	// they will be vailable from python level
  static PyMethodDef RungeKuttaTrackerClassMethods[] = {
		{ "track",     RungeKuttaTracker_track   ,METH_VARARGS,"Tracks a bunch through the fields and effects."},
    {NULL}
  };

	// defenition of the memebers of the python RungeKuttaTracker wrapper class
	// they will be vailable from python level
	static PyMemberDef RungeKuttaTrackerClassMembers [] = {
		{NULL}
	};

	//new python RungeKuttaTracker wrapper type definition
	static PyTypeObject pyORBIT_RungeKuttaTracker_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"RungeKuttaTracker", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) RungeKuttaTracker_del , /*tp_dealloc*/
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
		"The RungeKuttaTracker python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		RungeKuttaTrackerClassMethods, /* tp_methods */
		RungeKuttaTrackerClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) RungeKuttaTracker_init, /* tp_init */
		0, /* tp_alloc */
		RungeKuttaTracker_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyRungeKuttaTracker class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initRungeKuttaTracker(PyObject* module){
		if (PyType_Ready(&pyORBIT_RungeKuttaTracker_Type) < 0) return;
		Py_INCREF(&pyORBIT_RungeKuttaTracker_Type);
		PyModule_AddObject(module, "RungeKuttaTracker", (PyObject *)&pyORBIT_RungeKuttaTracker_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_tracker3dfield
}
