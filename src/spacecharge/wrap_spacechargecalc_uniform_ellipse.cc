#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_spacechargecalc_uniform_ellipse.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "SpaceChargeCalcUnifEllipse.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python SpaceChargeCalcUnifEllipse class definition
	//---------------------------------------------------------

	//constructor for python class wrapping SpaceChargeCalcUnifEllipse instance
	//It never will be called directly

	static PyObject* SpaceChargeCalcUnifEllipse_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  //initializator for python SpaceChargeCalcUnifEllipse class
  //this is implementation of the __init__ method SpaceChargeCalcUnifEllipse(nEllipses = 1])
  static int SpaceChargeCalcUnifEllipse_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
  	int nEllipses = 1;
		if(!PyArg_ParseTuple(args,"|i:__init__",&nEllipses)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcUnifEllipse - SpaceChargeCalcUnifEllipse([nEllipses = 1]) - constructor needs parameters.");
		}
		self->cpp_obj = new SpaceChargeCalcUnifEllipse(nEllipses);
		return 0;
	}
	
  //trackBunch(Bunch* bunch, double length)
  static PyObject* SpaceChargeCalcUnifEllipse_trackBunch(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalcUnifEllipse = (pyORBIT_Object*) self;
		SpaceChargeCalcUnifEllipse* cpp_SpaceChargeCalcUnifEllipse = (SpaceChargeCalcUnifEllipse*) pySpaceChargeCalcUnifEllipse->cpp_obj;
		PyObject* pyBunch;
		double length;
		if(!PyArg_ParseTuple(args,"Od:trackBunch",&pyBunch,&length)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcUnifEllipse.trackBunch(pyBunch,length) - method needs parameters.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcUnifEllipse.trackBunch(pyBunch,length) - pyBunch is not Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;		
		cpp_SpaceChargeCalcUnifEllipse->trackBunch(cpp_bunch,length);
		Py_INCREF(Py_None);
		return Py_None;
  }

  //-----------------------------------------------------
  //destructor for python SpaceChargeCalcUnifEllipse class (__del__ method).
  //-----------------------------------------------------
  static void SpaceChargeCalcUnifEllipse_del(pyORBIT_Object* self){
		SpaceChargeCalcUnifEllipse* cpp_SpaceChargeCalcUnifEllipse = (SpaceChargeCalcUnifEllipse*) self->cpp_obj;
		if(cpp_SpaceChargeCalcUnifEllipse != NULL){
			delete cpp_SpaceChargeCalcUnifEllipse;
		}
		self->ob_type->tp_free((PyObject*)self);
  }	
  
  // defenition of the methods of the python SpaceChargeCalcUnifEllipse wrapper class
  // they will be vailable from python level
  static PyMethodDef SpaceChargeCalcUnifEllipseClassMethods[] = {
		{ "trackBunch",  SpaceChargeCalcUnifEllipse_trackBunch, METH_VARARGS,"track the bunch - trackBunch(pyBunch,length)"},
		{NULL}
  };
  
  // defenition of the memebers of the python SpaceChargeCalcUnifEllipse wrapper class
  // they will be vailable from python level
  static PyMemberDef SpaceChargeCalcUnifEllipseClassMembers [] = {
		{NULL}
  };

	//new python SpaceChargeCalcUnifEllipse wrapper type definition
	static PyTypeObject pyORBIT_SpaceChargeCalcUnifEllipse_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SpaceChargeCalcUnifEllipse", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SpaceChargeCalcUnifEllipse_del , /*tp_dealloc*/
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
		"The SpaceChargeCalcUnifEllipse python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SpaceChargeCalcUnifEllipseClassMethods, /* tp_methods */
		SpaceChargeCalcUnifEllipseClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SpaceChargeCalcUnifEllipse_init, /* tp_init */
		0, /* tp_alloc */
		SpaceChargeCalcUnifEllipse_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pySpaceChargeCalcUnifEllipse class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initSpaceChargeCalcUniformEllipse(PyObject* module){
		if (PyType_Ready(&pyORBIT_SpaceChargeCalcUnifEllipse_Type) < 0) return;
		Py_INCREF(&pyORBIT_SpaceChargeCalcUnifEllipse_Type);
		PyModule_AddObject(module, "SpaceChargeCalcUnifEllipse", (PyObject *)&pyORBIT_SpaceChargeCalcUnifEllipse_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
