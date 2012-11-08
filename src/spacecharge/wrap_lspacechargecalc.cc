#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#
#include "wrap_lspacechargecalc.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "LSpaceChargeCalc.hh"

using namespace OrbitUtils;

namespace wrap_LSpaceChargeCalc{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python LSpaceChargeCalc class definition
	//---------------------------------------------------------

	//constructor for python class wrapping LSpaceChargeCalc instance
	//It never will be called directly

	static PyObject* LSpaceChargeCalc_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The LSpaceChargeCalc new has been called!"<<std::endl;
		return (PyObject *) self;
	}
	
  //initializator for python LSpaceChargeCalc class
  //this is implementation of the __init__ method LSpaceChargeCalc(double b_a, double length, int nMacrosMin, int useSpaceCharge, int nBins)
  static int LSpaceChargeCalc_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  
		pyORBIT_Object* pyLSpaceChargeCalc = (pyORBIT_Object*) self;
		LSpaceChargeCalc* cpp_LSpaceChargeCalc = (LSpaceChargeCalc*) pyLSpaceChargeCalc->cpp_obj;
	  
		double b_a = 1.0;
		double length = 1.0;
		int nMacrosMin;
		int useSpaceCharge;
		int nBins;
	  
		if(!PyArg_ParseTuple(args,"ddiii:arguments",&b_a,&length,&nMacrosMin,&useSpaceCharge,&nBins)){
			ORBIT_MPI_Finalize("PyLSpaceChargeCalc - LSpaceChargeCalc(b_a, length, nMacrosMin, useSpaceCharge, nBins) - constructor needs parameters.");
		}
	  
		self->cpp_obj = new LSpaceChargeCalc(b_a, length, nMacrosMin, useSpaceCharge, nBins);
		
		((LSpaceChargeCalc*) self->cpp_obj)->setPyWrapper((PyObject*) self);

		return 0;
	}

	
	//assignImpedance  A routine to import a python complex tuple and convert to c++ impedance array
	static PyObject* assignImpedance(PyObject *self, PyObject *args){
		
		pyORBIT_Object* pyLSpaceChargeCalc = (pyORBIT_Object*) self;
		LSpaceChargeCalc* cpp_LSpaceChargeCalc = (LSpaceChargeCalc*) pyLSpaceChargeCalc->cpp_obj;
		PyObject* py_cmplx_arr;
		
		if(!PyArg_ParseTuple(args,"O:get_complex_arr",&py_cmplx_arr)){
			ORBIT_MPI_Finalize("ERROR! You have to specify a parameter - array of complex numbers!");
		}
		if(PySequence_Check(py_cmplx_arr) != 1){
			ORBIT_MPI_Finalize("ERROR! You have to specify a parameter - array of complex numbers!");
		}		
		int size = PySequence_Size(py_cmplx_arr);
		Py_complex cmplx;
		PyObject* py_cmplx;	
		double real,imag;
		for(int i = 0; i < size; i++){
			py_cmplx = PySequence_Fast_GET_ITEM(py_cmplx_arr, i);
			if(!PyComplex_Check(py_cmplx)){
				ORBIT_MPI_Finalize("ERROR! No complex numbers!");
			}
			cmplx = PyComplex_AsCComplex(py_cmplx);
			real = cmplx.real;
			imag = cmplx.imag;
			cpp_LSpaceChargeCalc->assignImpedanceValue(i,real,imag);
		}
		
		Py_INCREF(Py_None);
		return Py_None;  
		
	}


	
	//assignImpedanceValue(int, real, real).  Wraps the LongSpaceChargeCalc routine assigning an impedance mode
	static PyObject* LSpaceChargeCalc_assignImpedanceValue(PyObject *self, PyObject *args){
		
		pyORBIT_Object* pyLSpaceChargeCalc = (pyORBIT_Object*) self;
		LSpaceChargeCalc* cpp_LSpaceChargeCalc = (LSpaceChargeCalc*) pyLSpaceChargeCalc->cpp_obj;
		
		int i = 0;
		double real= 0.0;
		double imag = 0.0;
		
		if(!PyArg_ParseTuple(args,"idd:arguments",&i,&real,&imag)){
			ORBIT_MPI_Finalize("PyLSpaceChargeCalc - assignImpedanceValue(i, real,imag) - constructor needs parameters.");
		}
		cpp_LSpaceChargeCalc->assignImpedanceValue(i,real,imag);

		Py_INCREF(Py_None);
		return Py_None;  
	}

	
//trackBunchBunch(Bunch* bunch)
  static PyObject* LSpaceChargeCalc_trackBunch(PyObject *self, PyObject *args){
		int nVars = PyTuple_Size(args);
		pyORBIT_Object* pyLSpaceChargeCalc = (pyORBIT_Object*) self;
		LSpaceChargeCalc* cpp_LSpaceChargeCalc = (LSpaceChargeCalc*) pyLSpaceChargeCalc->cpp_obj;
		PyObject* pyBunch;
		
		if(!PyArg_ParseTuple(args,"O:trackBunch",&pyBunch)){
			ORBIT_MPI_Finalize("PyLSpaceChargeCalc.trackBunch(pyBunch) - method needs parameters.");
		}
		
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PyLSpaceChargeCalc.trackBunch(pyBunch) - pyBunch is not Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;		
		cpp_LSpaceChargeCalc->trackBunch(cpp_bunch);
		
		Py_INCREF(Py_None);
		return Py_None;  
  }

  //-----------------------------------------------------
  //destructor for python LSpaceChargeCalc class (__del__ method).
  //-----------------------------------------------------
  static void LSpaceChargeCalc_del(pyORBIT_Object* self){
		LSpaceChargeCalc* cpp_LSpaceChargeCalc = (LSpaceChargeCalc*) self->cpp_obj;
		if(cpp_LSpaceChargeCalc != NULL){
			delete cpp_LSpaceChargeCalc;
		}
		self->ob_type->tp_free((PyObject*)self);
  }	
  
  // defenition of the methods of the python LSpaceChargeCalc wrapper class
  // they will be vailable from python level
  static PyMethodDef LSpaceChargeCalcClassMethods[] = {
		{ "trackBunch",  LSpaceChargeCalc_trackBunch, METH_VARARGS,"trackBunch the bunch - trackBunch(pyBunch)"},
		{ "assignImpedanceValue",  LSpaceChargeCalc_assignImpedanceValue, METH_VARARGS,"assigne the impedance for the ith mode - assignImpedanceValue(i,real,imag)"},
		{ "assignImpedance",  assignImpedance, METH_VARARGS,"assigne the impedance for the ith mode - assignImpedance(Z))"},
		{NULL}
  };
  
  // defenition of the members of the python LSpaceChargeCalc wrapper class
  // they will be vailable from python level
  static PyMemberDef LSpaceChargeCalcClassMembers [] = {
		{NULL}
  };

	//new python LSpaceChargeCalc wrapper type definition
	static PyTypeObject pyORBIT_LSpaceChargeCalc_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"LSpaceChargeCalc", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) LSpaceChargeCalc_del , /*tp_dealloc*/
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
		"The LSpaceChargeCalc python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		LSpaceChargeCalcClassMethods, /* tp_methods */
		LSpaceChargeCalcClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) LSpaceChargeCalc_init, /* tp_init */
		0, /* tp_alloc */
		LSpaceChargeCalc_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyLSpaceChargeCalc class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initLSpaceChargeCalc(PyObject* module){
		if (PyType_Ready(&pyORBIT_LSpaceChargeCalc_Type) < 0) return;
		Py_INCREF(&pyORBIT_LSpaceChargeCalc_Type);
		PyModule_AddObject(module, "LSpaceChargeCalc", (PyObject *)&pyORBIT_LSpaceChargeCalc_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
