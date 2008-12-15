#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_runge_kutta_tracker.hh"

#include <iostream>

#include "RungeKuttaTracker.hh"

using namespace TrackerRK4;
using namespace OrbitUtils;

namespace wrap_trackerrk4{

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
		((RungeKuttaTracker*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The RungeKuttaTracker __init__ has been called!"<<std::endl;
		return 0;
  }

	// spatialEps([eps]) - set or get the spatial accuracy for integration
  static PyObject* RungeKuttaTracker_spatialEps(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		double eps = - 1.0;
		if(!PyArg_ParseTuple(args,"|d:spatialEps",&eps)){
			error("PyRungeKuttaTracker - the call should be - spatialEps([eps in meters]).");
		}
		if(eps > 0.){
			cpp_RungeKuttaTracker->setSpatialEps(eps);
		}
		return Py_BuildValue("d",cpp_RungeKuttaTracker->getSpatialEps());
  }	
	
	// length([l in [m]]) - set or get the approximate length of element
  static PyObject* RungeKuttaTracker_length(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		double l = - 1.0;
		if(!PyArg_ParseTuple(args,"|d:length",&l)){
			error("PyRungeKuttaTracker - the call should be - length([l in meters]).");
		}
		if(l > 0.){
			cpp_RungeKuttaTracker->setLength(l);
		}
		return Py_BuildValue("d",cpp_RungeKuttaTracker->getLength());
  }	
	
	// timeStep() - returns the time step [seconds] in Runge-Kutta integration 
  static PyObject* RungeKuttaTracker_timeStep(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		return Py_BuildValue("d",cpp_RungeKuttaTracker->getTimeStep());
  }		
		
	// stepNumber([n_step]) - set or get the number of step during integration
  static PyObject* RungeKuttaTracker_stepsNumber(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		int n = - 1;
		if(!PyArg_ParseTuple(args,"|i:stepsNumber",&n)){
			error("PyRungeKuttaTracker - the call should be - stepsNumber([n]).");
		}
		if(n > 0){
			cpp_RungeKuttaTracker->setInitialStepsNumber(n);
			return Py_BuildValue("i",cpp_RungeKuttaTracker->getInitialStepsNumber());
		}
		return Py_BuildValue("i",cpp_RungeKuttaTracker->getStepsNumber());
  }	
	
	// entrancePlane([(a,b,c,d)]) - set or get the coeff. in a*x+b*y+c*z+d = 0 entrance plane
  static PyObject* RungeKuttaTracker_entrancePlane(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		double a = 0., b = 0., c = 0., d = 0.;
		if(!PyArg_ParseTuple(args,"|dddd:entrancePlane",&a,&b,&c,&d)){
			error("PyRungeKuttaTracker - the call should be - entrancePlane([(a,b,c,d)]).");
		}
		if(a == 0. && b == 0. && c == 0. && d == 0.){
			cpp_RungeKuttaTracker->getEntrPlane(a,b,c,d);
			return Py_BuildValue("(dddd)",a,b,c,d);
		}
		cpp_RungeKuttaTracker->setEntrPlane(a,b,c,d);
		return Py_BuildValue("(dddd)",a,b,c,d);
  }		
	
	// exitPlane([(a,b,c,d)]) - set or get the coeff. in a*x+b*y+c*z+d = 0 exit plane
  static PyObject* RungeKuttaTracker_exitPlane(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		double a = 0., b = 0., c = 0., d = 0.;
		if(!PyArg_ParseTuple(args,"|dddd:exitPlane",&a,&b,&c,&d)){
			error("PyRungeKuttaTracker - the call should be - exitPlane([(a,b,c,d)]).");
		}
		if(a == 0. && b == 0. && c == 0. && d == 0.){
			cpp_RungeKuttaTracker->getExitPlane(a,b,c,d);
			return Py_BuildValue("(dddd)",a,b,c,d);
		}
		cpp_RungeKuttaTracker->setExitPlane(a,b,c,d);
		return Py_BuildValue("(dddd)",a,b,c,d);
  }		
	
	// isOutside(x,y,z) - check if coordinates is inside the tracker and returns 0 or 1 
  static PyObject* RungeKuttaTracker_isOutside(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		double x, y, z;
		if(!PyArg_ParseTuple(args,"ddd:isOutside",&x,&y,&z)){
			error("PyRungeKuttaTracker - the call should be - isOutside(x,y,z).");
		}
		return Py_BuildValue("i",cpp_RungeKuttaTracker->isOutside(x,y,z));
  }		
	
	// isAfterEntrance(x,y,z) - check if coordinates is after entrance plane of the tracker and returns 0 or 1 
  static PyObject* RungeKuttaTracker_isAfterEntrance(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		double x, y, z;
		if(!PyArg_ParseTuple(args,"ddd:isAfterEntrance",&x,&y,&z)){
			error("PyRungeKuttaTracker - the call should be - isAfterEntrance(x,y,z).");
		}
		return Py_BuildValue("i",cpp_RungeKuttaTracker->isAfterEntrance(x,y,z));
  }	

	// isBeforeExit(x,y,z) - check if coordinates is before exit plane of the tracker and returns 0 or 1 
  static PyObject* RungeKuttaTracker_isBeforeExit(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		double x, y, z;
		if(!PyArg_ParseTuple(args,"ddd:isBeforeExit",&x,&y,&z)){
			error("PyRungeKuttaTracker - the call should be - isBeforeExit(x,y,z).");
		}
		return Py_BuildValue("i",cpp_RungeKuttaTracker->isBeforeExit(x,y,z));
  }	
	
	// trackBunch(bunch,field_source,exteranl_effects) - track ORBIT traditional bunch
  static PyObject* RungeKuttaTracker_trackBunch(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		PyObject *pyBunch = NULL;
		PyObject *pyFieldSource = NULL;
		PyObject *pyExtEffects = NULL;
		if(!PyArg_ParseTuple(args,"OO|O:trackBunch",&pyBunch,&pyFieldSource,&pyExtEffects)){
			error("PyRungeKuttaTracker - trackBunch(bunch,field_source[,exteranl_effects]) - parameters are needed.");
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
	
	// track(bunch,time, time_step, field_source,exteranl_effects) - track bunch with r and p vectors
  static PyObject* RungeKuttaTracker_track(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRungeKuttaTracker = (pyORBIT_Object*) self;
		RungeKuttaTracker* cpp_RungeKuttaTracker = (RungeKuttaTracker*) pyRungeKuttaTracker->cpp_obj;
		PyObject *pyBunch = NULL;
		PyObject *pyFieldSource = NULL;
		PyObject *pyExtEffects = NULL;
		double tm = 0.;
		double t_step = 0.;
		double t_begin = 0.;
		if(!PyArg_ParseTuple(args,"OdddO|O:track",&pyBunch,&t_begin,&tm,&t_step, &pyFieldSource,&pyExtEffects)){
			error("PyRungeKuttaTracker - track(bunch,time, time_step, field_source[,exteranl_effects]) - parameters are needed.");
		}
		Bunch* bunch = (Bunch*) ((pyORBIT_Object*) pyBunch)->cpp_obj;
		BaseFieldSource* fs = (BaseFieldSource*) ((pyORBIT_Object*) pyFieldSource)->cpp_obj;
		ExternalEffects* extEf = NULL;
		if(pyExtEffects != NULL){
			extEf = (ExternalEffects*) ((pyORBIT_Object*) pyExtEffects)->cpp_obj;
		}    
		cpp_RungeKuttaTracker->track(bunch,t_begin,tm,t_step,fs,extEf);
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
		{ "spatialEps",      RungeKuttaTracker_spatialEps,      METH_VARARGS,"spatialEps([eps in meter]) - spatial accuracy."},
		{ "length",          RungeKuttaTracker_length,          METH_VARARGS,"length([L in meter]) - approximate length."},
		{ "timeStep",        RungeKuttaTracker_timeStep,        METH_VARARGS,"returns time step in seconds."},
		{ "stepsNumber",     RungeKuttaTracker_stepsNumber,     METH_VARARGS,"sets or returns the number of steps in integration."},
		{ "entrancePlane",   RungeKuttaTracker_entrancePlane,   METH_VARARGS,"sets or returns (a,b,c,d) in a*x+b*y+c*z+d=0 for entrance."}, 
		{ "exitPlane",       RungeKuttaTracker_exitPlane,       METH_VARARGS,"sets or returns (a,b,c,d) in a*x+b*y+c*z+d=0 for exit."}, 
		{ "isOutside",       RungeKuttaTracker_isOutside,       METH_VARARGS,"returns 0 or 1"}, 
		{ "isAfterEntrance", RungeKuttaTracker_isAfterEntrance, METH_VARARGS,"returns 0 or 1"}, 
		{ "isBeforeExit",    RungeKuttaTracker_isBeforeExit,    METH_VARARGS,"returns 0 or 1"}, 
		{ "trackBunch",      RungeKuttaTracker_trackBunch,      METH_VARARGS,"Tracks a ORBIT bunch through the fields and effects."},
		{ "track",           RungeKuttaTracker_track,           METH_VARARGS,"Tracks a bunch r and p vectors through the fields and effects."},
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

//end of namespace wrap_trackerrk4
}
