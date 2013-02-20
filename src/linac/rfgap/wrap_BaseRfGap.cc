#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_BaseRfGap.hh"
#include "wrap_linacmodule.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "BaseRfGap.hh"

using namespace OrbitUtils;

namespace wrap_linac{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python BaseRfGap class definition
	//---------------------------------------------------------

	//constructor for python class wrapping BaseRfGap instance
	//It never will be called directly
	static PyObject* BaseRfGap_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The BaseRfGap new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  BaseRfGap class
  //this is implementation of the __init__ method
  static int BaseRfGap_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new BaseRfGap();	
		((BaseRfGap*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		return 0;
  }
			
	//trackBunch(Bunch* bunch)
  static PyObject* BaseRfGap_trackBunch(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBaseRfGap = (pyORBIT_Object*) self;
		BaseRfGap* cpp_BaseRfGap = (BaseRfGap*) pyBaseRfGap->cpp_obj;
		PyObject* pyBunch;
	  double frequency, e0tl, phase, ampl;		
		if(!PyArg_ParseTuple(args,"Odddd:trackBunch",&pyBunch,&frequency,&ampl,&e0tl,&phase)){
			ORBIT_MPI_Finalize("PyBaseRfGap - trackBunch(Bunch* bunch, freq, ampl, E0TL, phase) - parameters are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PyBaseRfGap - trackBunch(Bunch* bunch, freq, ampl, E0TL, phase) - first param. should be a Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		cpp_BaseRfGap->trackBunch(cpp_bunch,frequency,ampl,e0tl,phase);
		Py_INCREF(Py_None);
    return Py_None;	
	}		
	
  //-----------------------------------------------------
  //destructor for python BaseRfGap class (__del__ method).
  //-----------------------------------------------------
  static void BaseRfGap_del(pyORBIT_Object* self){
		//std::cerr<<"The BaseRfGap __del__ has been called!"<<std::endl;
		BaseRfGap* cpp_BaseRfGap = (BaseRfGap*) self->cpp_obj;
		delete cpp_BaseRfGap;
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python BaseRfGap wrapper class
	// they will be vailable from python level
  static PyMethodDef BaseRfGapClassMethods[] = {
		{ "trackBunch",     BaseRfGap_trackBunch,    METH_VARARGS,"tracks the Bunch through the RF gap."},
    {NULL}
  };

	// defenition of the memebers of the python BaseRfGap wrapper class
	// they will be vailable from python level
	static PyMemberDef BaseRfGapClassMembers [] = {
		{NULL}
	};

	//new python BaseRfGap wrapper type definition
	static PyTypeObject pyORBIT_BaseRfGap_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"BaseRfGap", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) BaseRfGap_del , /*tp_dealloc*/
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
		"The BaseRfGap python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		BaseRfGapClassMethods, /* tp_methods */
		BaseRfGapClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) BaseRfGap_init, /* tp_init */
		0, /* tp_alloc */
		BaseRfGap_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyBaseRfGap class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initBaseRfGap(PyObject* module){
		if (PyType_Ready(&pyORBIT_BaseRfGap_Type) < 0) return;
		Py_INCREF(&pyORBIT_BaseRfGap_Type);
		PyModule_AddObject(module, "BaseRfGap", (PyObject *)&pyORBIT_BaseRfGap_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
