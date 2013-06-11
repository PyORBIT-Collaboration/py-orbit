#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_RfGapThreePointTTF.hh"
#include "wrap_linacmodule.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "wrap_utils.hh"
#include "RfGapThreePointTTF.hh"
#include "OU_Polynomial.hh"

using namespace OrbitUtils;

namespace wrap_linac{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python RfGapThreePointTTF class definition
	//---------------------------------------------------------

	//constructor for python class wrapping RfGapThreePointTTF instance
	//It never will be called directly
	static PyObject* RfGapThreePointTTF_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The RfGapThreePointTTF new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  RfGapThreePointTTF class
  //this is implementation of the __init__ method
  static int RfGapThreePointTTF_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new RfGapThreePointTTF();	
  	((RfGapThreePointTTF*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		return 0;
  }
			
	//trackBunch(Bunch* bunch)
  static PyObject* RfGapThreePointTTF_trackBunch(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRfGapThreePointTTF = (pyORBIT_Object*) self;
		RfGapThreePointTTF* cpp_RfGapThreePointTTF = (RfGapThreePointTTF*) pyRfGapThreePointTTF->cpp_obj;
		PyObject* pyBunch;
	  double dz, Em, E0, Ep, rf_frequency, phase;		
		if(!PyArg_ParseTuple(args,"Odddddd:trackBunch",&pyBunch,&dz,&Em,&E0,&Ep,&rf_frequency,&phase)){
			ORBIT_MPI_Finalize("PyRfGapThreePointTTF - trackBunch(Bunch* bunch, dz, Em, E0, Ep, frequency, phase) - parameters are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PyRfGapThreePointTTF - trackBunch(Bunch* bunch, dz, Em, E0, Ep, frequency, phase ) - first param. should be a Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		cpp_RfGapThreePointTTF->trackBunch(cpp_bunch,dz,Em,E0,Ep,rf_frequency,phase);
		Py_INCREF(Py_None);
    return Py_None;	
	}		
	
	//getT_TTF() returns the symmetrical part of the three-point transit time factor
  static PyObject* RfGapThreePointTTF_getT_TTF(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRfGapThreePointTTF = (pyORBIT_Object*) self;
    RfGapThreePointTTF* cpp_RfGapThreePointTTF = (RfGapThreePointTTF*) pyRfGapThreePointTTF->cpp_obj;    
    double dz, a, b, kappa;
		if(!PyArg_ParseTuple(args,"dddd:getT_TTF",&dz,&a,&b,&kappa)){
			ORBIT_MPI_Finalize("PyRfGapThreePointTTF - getT_TTF(dz,a,b,kappa) - parameters are needed.");
		}		
		double ttf = cpp_RfGapThreePointTTF->Tttf(dz,a,b,kappa);
   return  Py_BuildValue("d",ttf);
	}	
	
	//getS_TTF() returns the asymmetrical part of the three-point transit time factor
  static PyObject* RfGapThreePointTTF_getS_TTF(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRfGapThreePointTTF = (pyORBIT_Object*) self;
    RfGapThreePointTTF* cpp_RfGapThreePointTTF = (RfGapThreePointTTF*) pyRfGapThreePointTTF->cpp_obj;    
    double dz, a, b, kappa;
		if(!PyArg_ParseTuple(args,"dddd:getS_TTF",&dz,&a,&b,&kappa)){
			ORBIT_MPI_Finalize("PyRfGapThreePointTTF - getS_TTF(dz,a,b,kappa) - parameters are needed.");
		}		
		double ttf = cpp_RfGapThreePointTTF->Sttf(dz,a,b,kappa);
   return  Py_BuildValue("d",ttf);
	}	
	
	//getTp_TTF() returns the symmetrical part of the three-point transit time factor
  static PyObject* RfGapThreePointTTF_getTp_TTF(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRfGapThreePointTTF = (pyORBIT_Object*) self;
    RfGapThreePointTTF* cpp_RfGapThreePointTTF = (RfGapThreePointTTF*) pyRfGapThreePointTTF->cpp_obj;    
    double dz, a, b, kappa;
		if(!PyArg_ParseTuple(args,"dddd:getTp_TTF",&dz,&a,&b,&kappa)){
			ORBIT_MPI_Finalize("PyRfGapThreePointTTF - getTp_TTF(dz,a,b,kappa) - parameters are needed.");
		}		
		double ttf = cpp_RfGapThreePointTTF->Tpttf(dz,a,b,kappa);
   return  Py_BuildValue("d",ttf);
	}	
	
	//getSp_TTF() returns the asymmetrical part of the three-point transit time factor
  static PyObject* RfGapThreePointTTF_getSp_TTF(PyObject *self, PyObject *args){
    pyORBIT_Object* pyRfGapThreePointTTF = (pyORBIT_Object*) self;
    RfGapThreePointTTF* cpp_RfGapThreePointTTF = (RfGapThreePointTTF*) pyRfGapThreePointTTF->cpp_obj;    
    double dz, a, b, kappa;
		if(!PyArg_ParseTuple(args,"dddd:getSp_TTF",&dz,&a,&b,&kappa)){
			ORBIT_MPI_Finalize("PyRfGapThreePointTTF - getSp_TTF(dz,a,b,kappa) - parameters are needed.");
		}		
		double ttf = cpp_RfGapThreePointTTF->Spttf(dz,a,b,kappa);
   return  Py_BuildValue("d",ttf);
	}	
	
  //-----------------------------------------------------
  //destructor for python RfGapThreePointTTF class (__del__ method).
  //-----------------------------------------------------
  static void RfGapThreePointTTF_del(pyORBIT_Object* self){
		//std::cerr<<"The RfGapThreePointTTF __del__ has been called!"<<std::endl;
		RfGapThreePointTTF* cpp_RfGapThreePointTTF = (RfGapThreePointTTF*) self->cpp_obj;	
		delete cpp_RfGapThreePointTTF;
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python RfGapThreePointTTF wrapper class
	// they will be vailable from python level
  static PyMethodDef RfGapThreePointTTFClassMethods[] = {
		{ "trackBunch",     RfGapThreePointTTF_trackBunch,    METH_VARARGS,"tracks the Bunch through the RF gap trackBunch(bunch,dz,Em,E0,Ep,freq,phase)."},
		{ "getT_TTF",       RfGapThreePointTTF_getT_TTF,      METH_VARARGS,"get the T TTF value of the gap model."},
		{ "getS_TTF",       RfGapThreePointTTF_getS_TTF,      METH_VARARGS,"get the S TTF value of the gap model."},
		{ "getTp_TTF",      RfGapThreePointTTF_getTp_TTF,     METH_VARARGS,"get the Tp TTF value of the gap model."},
		{ "getSp_TTF",      RfGapThreePointTTF_getSp_TTF,     METH_VARARGS,"get the Sp TTF value of the gap model."},
    {NULL}
  };

	// defenition of the memebers of the python RfGapThreePointTTF wrapper class
	// they will be vailable from python level
	static PyMemberDef RfGapThreePointTTFClassMembers [] = {
		{NULL}
	};

	//new python RfGapThreePointTTF wrapper type definition
	static PyTypeObject pyORBIT_RfGapThreePointTTF_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"RfGapThreePointTTF", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) RfGapThreePointTTF_del , /*tp_dealloc*/
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
		"The RfGapThreePointTTF python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		RfGapThreePointTTFClassMethods, /* tp_methods */
		RfGapThreePointTTFClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) RfGapThreePointTTF_init, /* tp_init */
		0, /* tp_alloc */
		RfGapThreePointTTF_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyRfGapThreePointTTF class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initRfGapThreePointTTF(PyObject* module){
		if (PyType_Ready(&pyORBIT_RfGapThreePointTTF_Type) < 0) return;
		Py_INCREF(&pyORBIT_RfGapThreePointTTF_Type);
		PyModule_AddObject(module, "RfGapThreePointTTF", (PyObject *)&pyORBIT_RfGapThreePointTTF_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_linac
}
