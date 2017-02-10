#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_RfGapTTF_slow.hh"
#include "wrap_linacmodule.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "wrap_utils.hh"
#include "RfGapTTF_slow.hh"
#include "OU_Polynomial.hh"

using namespace OrbitUtils;

namespace wrap_linac{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python RfGapTTF_slow class definition
	//---------------------------------------------------------

	//constructor for python class wrapping RfGapTTF_slow instance
	//It never will be called directly
	static PyObject* RfGapTTF_slow_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The RfGapTTF_slow new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  RfGapTTF_slow class
  //this is implementation of the __init__ method
  static int RfGapTTF_slow_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new RfGapTTF_slow();	
		((RfGapTTF_slow*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		return 0;
  }
			
	//trackBunch(Bunch* bunch)
  static PyObject* RfGapTTF_slow_trackBunch(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRfGapTTF_slow = (pyORBIT_Object*) self;
		RfGapTTF_slow* cpp_RfGapTTF_slow = (RfGapTTF_slow*) pyRfGapTTF_slow->cpp_obj;
		PyObject* pyBunch;
		PyObject* pyPolyT;
		PyObject* pyPolyTP;
		PyObject* pyPolyS;
		PyObject* pyPolySP;
	  double frequency, E0L, phase;		
		if(!PyArg_ParseTuple(args,"OdddOOOO:trackBunch",&pyBunch,&frequency,&E0L,&phase,&pyPolyT,&pyPolyS,&pyPolyTP,&pyPolySP)){
			ORBIT_MPI_Finalize("PyRfGapTTF_slow - trackBunch(Bunch* bunch, frequency, E0L, phase, TTFs TSTpSp) - parameters are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PyRfGapTTF_slow - trackBunch(Bunch* bunch, frequency, E0L, phase, TTFs TSTpSp) - first param. should be a Bunch.");
		}
		PyObject* pyORBIT_Polynomial_Type = wrap_orbit_utils::getOrbitUtilsType("Polynomial");
		if(!PyObject_IsInstance(pyPolyT,pyORBIT_Polynomial_Type)){
			ORBIT_MPI_Finalize("PyRfGapTTF_slow - trackBunch(...) - the last parameters should be a Polynomial.");
		}			
		if(!PyObject_IsInstance(pyPolyS,pyORBIT_Polynomial_Type)){
			ORBIT_MPI_Finalize("PyRfGapTTF_slow - trackBunch(...) - the last parameters should be a Polynomial.");
		}			
		if(!PyObject_IsInstance(pyPolyTP,pyORBIT_Polynomial_Type)){
			ORBIT_MPI_Finalize("PyRfGapTTF_slow - trackBunch(...) - the last parameters should be a Polynomial.");
		}			
		if(!PyObject_IsInstance(pyPolySP,pyORBIT_Polynomial_Type)){
			ORBIT_MPI_Finalize("PyRfGapTTF_slow - trackBunch(...) - the last parameters should be a Polynomial.");
		}			
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		Polynomial* cpp_polyT = (Polynomial*) ((pyORBIT_Object*)pyPolyT)->cpp_obj;
		Polynomial* cpp_polyS = (Polynomial*) ((pyORBIT_Object*)pyPolyS)->cpp_obj;
		Polynomial* cpp_polyTp = (Polynomial*) ((pyORBIT_Object*)pyPolyTP)->cpp_obj;
		Polynomial* cpp_polySp = (Polynomial*) ((pyORBIT_Object*)pyPolySP)->cpp_obj;		
		cpp_RfGapTTF_slow->trackBunch(cpp_bunch,frequency,E0L,phase,cpp_polyT,cpp_polyS,cpp_polyTp,cpp_polySp);
		Py_INCREF(Py_None);
    return Py_None;	
	}		
		
  //-----------------------------------------------------
  //destructor for python RfGapTTF_slow class (__del__ method).
  //-----------------------------------------------------
  static void RfGapTTF_slow_del(pyORBIT_Object* self){
		//std::cerr<<"The RfGapTTF_slow __del__ has been called!"<<std::endl;
		RfGapTTF_slow* cpp_RfGapTTF_slow = (RfGapTTF_slow*) self->cpp_obj;	
		delete cpp_RfGapTTF_slow;
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python RfGapTTF_slow wrapper class
	// they will be vailable from python level
  static PyMethodDef RfGapTTF_slowClassMethods[] = {
		{ "trackBunch",     RfGapTTF_slow_trackBunch,    METH_VARARGS,"tracks the Bunch through the RF gap trackBunch(bunch,E0,phase)."},
    {NULL}
  };

	// defenition of the memebers of the python RfGapTTF_slow wrapper class
	// they will be vailable from python level
	static PyMemberDef RfGapTTF_slowClassMembers [] = {
		{NULL}
	};

	//new python RfGapTTF_slow wrapper type definition
	static PyTypeObject pyORBIT_RfGapTTF_slow_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"RfGapTTF_slow", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) RfGapTTF_slow_del , /*tp_dealloc*/
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
		"The RfGapTTF_slow python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		RfGapTTF_slowClassMethods, /* tp_methods */
		RfGapTTF_slowClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) RfGapTTF_slow_init, /* tp_init */
		0, /* tp_alloc */
		RfGapTTF_slow_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyRfGapTTF_slow class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initRfGapTTF_slow(PyObject* module){
		if (PyType_Ready(&pyORBIT_RfGapTTF_slow_Type) < 0) return;
		Py_INCREF(&pyORBIT_RfGapTTF_slow_Type);
		PyModule_AddObject(module, "RfGapTTF_slow", (PyObject *)&pyORBIT_RfGapTTF_slow_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_linac
}
