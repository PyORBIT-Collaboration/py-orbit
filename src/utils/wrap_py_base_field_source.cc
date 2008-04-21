#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_py_base_field_source.hh"

#include <iostream>

#include "PyBaseFieldSource.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_utils_py_base_field_source{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python PyBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping PyBaseFieldSource instance
	//It never will be called directly
	static PyObject* PyBaseFieldSource_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  PyBaseFieldSource class
  //this is implementation of the __init__ method
  static int PyBaseFieldSource_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new PyBaseFieldSource((PyObject*) self);
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void PyBaseFieldSource_del(pyORBIT_Object* self){
		//std::cerr<<"The PyBaseFieldSource __del__ has been called!"<<std::endl;
		delete ((PyBaseFieldSource*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef PyBaseFieldSourceClassMethods[] = {
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef PyBaseFieldSourceClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_PyBaseFieldSource_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"PyBaseFieldSource", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) PyBaseFieldSource_del , /*tp_dealloc*/
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
		"The PyBaseFieldSource python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		PyBaseFieldSourceClassMethods, /* tp_methods */
		PyBaseFieldSourceClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) PyBaseFieldSource_init, /* tp_init */
		0, /* tp_alloc */
		PyBaseFieldSource_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initPyBaseFieldSource(PyObject* module){
		if (PyType_Ready(&pyORBIT_PyBaseFieldSource_Type) < 0) return;
		Py_INCREF(&pyORBIT_PyBaseFieldSource_Type);
		PyModule_AddObject(module, "PyBaseFieldSource", (PyObject *)&pyORBIT_PyBaseFieldSource_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_py_base_field_source
}
