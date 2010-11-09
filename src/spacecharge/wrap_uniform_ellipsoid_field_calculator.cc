#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "UniformEllipsoidFieldCalculator.hh"

#include "wrap_uniform_ellipsoid_field_calculator.hh"
#include "wrap_spacecharge.hh"

#include <iostream>

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python UniformEllipsoidFieldCalculator class definition
	//---------------------------------------------------------

	//constructor for python class wrapping UniformEllipsoidFieldCalculator instance
	//It never will be called directly
	static PyObject* UniformEllipsoidFieldCalculator_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  UniformEllipsoidFieldCalculator class
  //this is implementation of the __init__ method
  static int UniformEllipsoidFieldCalculator_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new UniformEllipsoidFieldCalculator();
		((UniformEllipsoidFieldCalculator*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		return 0;
  }
		
	/** Sets the half-axis of the ellipsoid and maximal values of radius*/
	static PyObject* UniformEllipsoidFieldCalculator_setEllipsoid(PyObject *self, PyObject *args){
		pyORBIT_Object* pyUniformEllipsoidFieldCalculator = (pyORBIT_Object*) self;
		UniformEllipsoidFieldCalculator* cpp_UniformEllipsoidFieldCalculator = (UniformEllipsoidFieldCalculator*) pyUniformEllipsoidFieldCalculator->cpp_obj;
		double a,b,c,r_max;
		if(!PyArg_ParseTuple(args,"dddd:setEllipsoid",&a,&b,&c,&r_max)){		
			ORBIT_MPI_Finalize("PyUniformEllipsoidFieldCalculator.setEllipsoid(a,b,c,r_max) - method needs parameters.");
		}			
		cpp_UniformEllipsoidFieldCalculator->setEllipsoid(a,b,c,r_max);
		Py_INCREF(Py_None);
    return Py_None;
	}
	
	/** Calculates the field components */
	static PyObject* UniformEllipsoidFieldCalculator_calcField(PyObject *self, PyObject *args){
		pyORBIT_Object* pyUniformEllipsoidFieldCalculator = (pyORBIT_Object*) self;
		UniformEllipsoidFieldCalculator* cpp_UniformEllipsoidFieldCalculator = (UniformEllipsoidFieldCalculator*) pyUniformEllipsoidFieldCalculator->cpp_obj;
		double x,y,z;
		if(!PyArg_ParseTuple(args,"ddd:calcField",&x,&y,&z)){		
			ORBIT_MPI_Finalize("PyUniformEllipsoidFieldCalculator.calcField(x,y,z) - method needs parameters.");
		}	
		double x2,y2,z2;
		double ex,ey,ez;
		x2 = x*x; y2 = y*y; z2 = z*z;
		cpp_UniformEllipsoidFieldCalculator->calcField(x,y,z,x2,y2,z2,ex,ey,ez);
		return Py_BuildValue("(ddd)",ex,ey,ez);
	}
	
  //-----------------------------------------------------
  //destructor for python UniformEllipsoidFieldCalculator class (__del__ method).
  //-----------------------------------------------------
  static void UniformEllipsoidFieldCalculator_del(pyORBIT_Object* self){
		UniformEllipsoidFieldCalculator* cpp_UniformEllipsoidFieldCalculator = (UniformEllipsoidFieldCalculator*) self->cpp_obj;
		if(cpp_UniformEllipsoidFieldCalculator != NULL){
			delete cpp_UniformEllipsoidFieldCalculator;
		}
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python UniformEllipsoidFieldCalculator wrapper class
	// they will be vailable from python level
  static PyMethodDef UniformEllipsoidFieldCalculatorClassMethods[] = {
		{ "setEllipsoid",        UniformEllipsoidFieldCalculator_setEllipsoid,        METH_VARARGS,"sets the half-axis of the ellipsoid and maximal values of radius"},
		{ "calcField",           UniformEllipsoidFieldCalculator_calcField,           METH_VARARGS,"returns (ex,ey,ez) for (x,y,z) input"},
    {NULL}
  };

	// defenition of the memebers of the python UniformEllipsoidFieldCalculator wrapper class
	// they will be vailable from python level
	static PyMemberDef UniformEllipsoidFieldCalculatorClassMembers [] = {
		{NULL}
	};

	//new python UniformEllipsoidFieldCalculator wrapper type definition
	static PyTypeObject pyORBIT_UniformEllipsoidFieldCalculator_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"UniformEllipsoidFieldCalculator", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) UniformEllipsoidFieldCalculator_del , /*tp_dealloc*/
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
		"The UniformEllipsoidFieldCalculator python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		UniformEllipsoidFieldCalculatorClassMethods, /* tp_methods */
		UniformEllipsoidFieldCalculatorClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) UniformEllipsoidFieldCalculator_init, /* tp_init */
		0, /* tp_alloc */
		UniformEllipsoidFieldCalculator_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyUniformEllipsoidFieldCalculator class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initUniformEllipsoidFieldCalculator(PyObject* module){
		if (PyType_Ready(&pyORBIT_UniformEllipsoidFieldCalculator_Type) < 0) return;
		Py_INCREF(&pyORBIT_UniformEllipsoidFieldCalculator_Type);
		PyModule_AddObject(module, "UniformEllipsoidFieldCalculator", (PyObject *)&pyORBIT_UniformEllipsoidFieldCalculator_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
