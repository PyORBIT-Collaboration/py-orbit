#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_cpp_base_field_source.hh"

#include <iostream>
#include <string>

#include "BaseFieldSource.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_tracker3dfield_cpp_base_field_source{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python BaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping BaseFieldSource instance
	//It never will be called directly
	static PyObject* BaseFieldSource_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  BaseFieldSource class
  //this is implementation of the __init__ method
  static int BaseFieldSource_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new BaseFieldSource();
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void BaseFieldSource_del(pyORBIT_Object* self){
		//std::cerr<<"The BaseFieldSource __del__ has been called!"<<std::endl;
		delete ((BaseFieldSource*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef BaseFieldSourceClassMethods[] = {
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef BaseFieldSourceClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_BaseFieldSource_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"BaseFieldSource", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) BaseFieldSource_del , /*tp_dealloc*/
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
		"The BaseFieldSource python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		BaseFieldSourceClassMethods, /* tp_methods */
		BaseFieldSourceClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) BaseFieldSource_init, /* tp_init */
		0, /* tp_alloc */
		BaseFieldSource_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initBaseFieldSource(PyObject* module){
		if (PyType_Ready(&pyORBIT_BaseFieldSource_Type) < 0) return;
		Py_INCREF(&pyORBIT_BaseFieldSource_Type);
		PyModule_AddObject(module, "BaseFieldSource", (PyObject *)&pyORBIT_BaseFieldSource_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
