#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_MatrixRfGap.hh"
#include "wrap_linacmodule.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "MatrixRfGap.hh"

using namespace OrbitUtils;

namespace wrap_linac{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python MatrixRfGap class definition
	//---------------------------------------------------------

	//constructor for python class wrapping MatrixRfGap instance
	//It never will be called directly
	static PyObject* MatrixRfGap_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  MatrixRfGap class
  //this is implementation of the __init__ method
  static int MatrixRfGap_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new MatrixRfGap();	
		((MatrixRfGap*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		return 0;
  }
			
	//trackBunch(Bunch* bunch)
  static PyObject* MatrixRfGap_trackBunch(PyObject *self, PyObject *args){
    pyORBIT_Object* pyMatrixRfGap = (pyORBIT_Object*) self;
		MatrixRfGap* cpp_MatrixRfGap = (MatrixRfGap*) pyMatrixRfGap->cpp_obj;
		PyObject* pyBunch;
	  double frequency, e0tl, phase;		
		if(!PyArg_ParseTuple(args,"Oddd:trackBunch",&pyBunch,&frequency,&e0tl,&phase)){
			ORBIT_MPI_Finalize("PyMatrixRfGap - trackBunch(Bunch* bunch, freq, E0TL, phase) - parameters are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PyMatrixRfGap - trackBunch(Bunch* bunch, freq, E0TL, phase) - first param. should be a Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		cpp_MatrixRfGap->trackBunch(cpp_bunch,frequency,e0tl,phase);
		Py_INCREF(Py_None);
    return Py_None;	
	}		
	
  //-----------------------------------------------------
  //destructor for python MatrixRfGap class (__del__ method).
  //-----------------------------------------------------
  static void MatrixRfGap_del(pyORBIT_Object* self){
		//std::cerr<<"The MatrixRfGap __del__ has been called!"<<std::endl;
		MatrixRfGap* cpp_MatrixRfGap = (MatrixRfGap*) self->cpp_obj;
		delete cpp_MatrixRfGap;
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python MatrixRfGap wrapper class
	// they will be vailable from python level
  static PyMethodDef MatrixRfGapClassMethods[] = {
		{ "trackBunch",     MatrixRfGap_trackBunch,    METH_VARARGS,"tracks the Bunch through the RF gap."},
    {NULL}
  };

	// defenition of the memebers of the python MatrixRfGap wrapper class
	// they will be vailable from python level
	static PyMemberDef MatrixRfGapClassMembers [] = {
		{NULL}
	};

	//new python MatrixRfGap wrapper type definition
	static PyTypeObject pyORBIT_MatrixRfGap_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"MatrixRfGap", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) MatrixRfGap_del , /*tp_dealloc*/
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
		"The MatrixRfGap python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		MatrixRfGapClassMethods, /* tp_methods */
		MatrixRfGapClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) MatrixRfGap_init, /* tp_init */
		0, /* tp_alloc */
		MatrixRfGap_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyMatrixRfGap class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initMatrixRfGap(PyObject* module){
		if (PyType_Ready(&pyORBIT_MatrixRfGap_Type) < 0) return;
		Py_INCREF(&pyORBIT_MatrixRfGap_Type);
		PyModule_AddObject(module, "MatrixRfGap", (PyObject *)&pyORBIT_MatrixRfGap_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
