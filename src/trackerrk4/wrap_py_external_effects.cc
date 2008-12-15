#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_py_external_effects.hh"

#include <iostream>
#include <string>

#include "PyExternalEffects.hh"

using namespace OrbitUtils;
using namespace Tracker3DField;

namespace wrap_tracker3dfield_py_external_effects{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python PyExternalEffects class definition
	//---------------------------------------------------------

	//constructor for python class wrapping PyExternalEffects instance
	//It never will be called directly
	static PyObject* PyExternalEffects_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppExternalEffects class
  //this is implementation of the __init__ method
  static int PyExternalEffects_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new PyExternalEffects((PyObject*) self);
    return 0;
  }
  
  
  
	
	// name([name]) - sets or returns the name of the External Effeects class 
  static PyObject* PyExternalEffects_name(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPyExternalEffects= (pyORBIT_Object*) self;
		PyExternalEffects* cpp_PyExternalEffects = (PyExternalEffects*) pyPyExternalEffects->cpp_obj;
    const char* name = NULL;
    if(!PyArg_ParseTuple(	args,"|s:name",&name)){
      error("PyExternalEffects - call should be - name([name]).");
    }
		if(name != NULL){
      std::string name_str(name);
      cpp_PyExternalEffects->setName(name_str);
		}
		return Py_BuildValue("s",cpp_PyExternalEffects->getName().c_str());
  }	
	
  //-----------------------------------------------------
  //destructor for python PyExternalEffects class (__del__ method).
  //-----------------------------------------------------
  static void PyExternalEffects_del(pyORBIT_Object* self){
		//std::cerr<<"The PyExternalEffects __del__ has been called!"<<std::endl;
		delete ((PyExternalEffects*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyExternalEffects wrapper class
	// they will be vailable from python level
  static PyMethodDef PyExternalEffectsClassMethods[] = {
		{ "name",         PyExternalEffects_name,         METH_VARARGS,"Sets or returns the name of effects."},
    {NULL}
  };

	// defenition of the memebers of the python PyExternalEffects wrapper class
	// they will be vailable from python level
	static PyMemberDef PyExternalEffectsClassMembers [] = {
		{NULL}
	};

	//new python PyExternalEffects wrapper type definition
	static PyTypeObject pyORBIT_PyExternalEffects_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"PyExternalEffects", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) PyExternalEffects_del , /*tp_dealloc*/
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
		"The PyExternalEffects python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		PyExternalEffectsClassMethods, /* tp_methods */
		PyExternalEffectsClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) PyExternalEffects_init, /* tp_init */
		0, /* tp_alloc */
		PyExternalEffects_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyExternalEffects class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initPyExternalEffects(PyObject* module){
		if (PyType_Ready(&pyORBIT_PyExternalEffects_Type) < 0) return;
		Py_INCREF(&pyORBIT_PyExternalEffects_Type);
		PyModule_AddObject(module, "PyExternalEffects", (PyObject *)&pyORBIT_PyExternalEffects_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_tracker3dfield
}
