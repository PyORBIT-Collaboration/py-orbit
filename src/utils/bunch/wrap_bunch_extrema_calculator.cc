#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_matrix.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "BunchExtremaCalculator.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_utils_bunch{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python BunchExtremaCalculator class definition
	//---------------------------------------------------------

	//constructor for python class wrapping BunchExtremaCalculator instance
	//It never will be called directly
	static PyObject* BunchExtremaCalculator_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  BunchExtremaCalculator class
  //this is implementation of the __init__ method
  static int BunchExtremaCalculator_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new BunchExtremaCalculator();
    return 0;
  }

	//Calculates xMin,xMax,  yMin,yMax, zMin,zMax,
  static PyObject* BunchExtremaCalculator_extremaXYZ(PyObject *self, PyObject *args){
		PyObject *pyIn;
		if(!PyArg_ParseTuple(args,"O:extremaXYZ",&pyIn)){
			error("PyBunchExtremaCalculator - extremaXYZ(Bunch) - Bunch is needed.");
		}			
		PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
		if((!PyObject_IsInstance(pyIn,pyBunchType))){
			error("PyBunchExtremaCalculator - extremaXYZ(Bunch) - input parameter is not Bunch");
		}		
		Bunch* bunch = (Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj;
		BunchExtremaCalculator* cpp_BunchExtremaCalculator = (BunchExtremaCalculator*) (((pyORBIT_Object*) self)->cpp_obj);
		double xMin,xMax,yMin,yMax,zMin,zMax;
		cpp_BunchExtremaCalculator->getExtremaXYZ(bunch, xMin, xMax, yMin, yMax, zMin, zMax);	
		return Py_BuildValue("(d,d,d,d,d,d)", xMin, xMax, yMin, yMax, zMin, zMax);
  }

  //-----------------------------------------------------
  //destructor for python BunchExtremaCalculator class (__del__ method).
  //-----------------------------------------------------
  static void BunchExtremaCalculator_del(pyORBIT_Object* self){
		delete ((BunchExtremaCalculator*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python BunchExtremaCalculator wrapper class
	// they will be vailable from python level
  static PyMethodDef BunchExtremaCalculatorClassMethods[] = {
    { "extremaXYZ",   BunchExtremaCalculator_extremaXYZ  ,METH_VARARGS,"Returns tuple with (xMin, xMax, yMin, yMax, zMin, zMax)"},
    {NULL}
  };

	// defenition of the memebers of the python BunchExtremaCalculator wrapper class
	// they will be vailable from python level
	static PyMemberDef BunchExtremaCalculatorClassMembers [] = {
		{NULL}
	};

	//new python BunchExtremaCalculator wrapper type definition
	static PyTypeObject pyORBIT_BunchExtremaCalculator_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"BunchExtremaCalculator", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) BunchExtremaCalculator_del , /*tp_dealloc*/
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
		"The BunchExtremaCalculator python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		BunchExtremaCalculatorClassMethods, /* tp_methods */
		BunchExtremaCalculatorClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) BunchExtremaCalculator_init, /* tp_init */
		0, /* tp_alloc */
		BunchExtremaCalculator_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the pyBunchExtremaCalculator class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initBunchExtremaCalculator(PyObject* module){
		if (PyType_Ready(&pyORBIT_BunchExtremaCalculator_Type) < 0) return;
		Py_INCREF(&pyORBIT_BunchExtremaCalculator_Type);
		PyModule_AddObject(module, "BunchExtremaCalculator", (PyObject *)&pyORBIT_BunchExtremaCalculator_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_martix
}
