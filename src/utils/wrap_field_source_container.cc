#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_field_source_container.hh"

#include <iostream>
#include <string>

#include "FieldSourceContainer.hh"
#include "BaseFieldSource.hh"


using namespace OrbitUtils;


namespace wrap_field_source_container{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif


	//constructor for python class wrapping CppFieldSource instance
	//It never will be called directly
	static PyObject* FieldSourceContainer_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  //this is implementation of the __init__ method
  static int FieldSourceContainer_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  self->cpp_obj =  new  FieldSourceContainer();	  
	  ((FieldSourceContainer*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
  static PyObject* FieldSourceContainer_AddFieldSource(PyObject *self, PyObject *args){
	  FieldSourceContainer* cpp_FieldSourceContainer = (FieldSourceContainer*)((pyORBIT_Object*) self)->cpp_obj;
	  BaseFieldSource* fs;
	  PyObject* pyfs;
		if(!PyArg_ParseTuple(	args,"O:",&pyfs))
			error(" AddEffect(external effect) - parameter is needed");
		else {
			fs = (BaseFieldSource*) ((pyORBIT_Object*) pyfs)->cpp_obj;
			cpp_FieldSourceContainer->AddFieldSource(fs);
		}
		Py_INCREF(Py_None);
		return Py_None;
  }
  
  static PyObject* FieldSourceContainer_getFields(PyObject *self, PyObject *args){
	  FieldSourceContainer* cpp_fields = (FieldSourceContainer*)((pyORBIT_Object*) self)->cpp_obj;
		
		double x,y,z,t,Ex,Ey,Ez,Bx,By,Bz;
		if(!PyArg_ParseTuple(	args,"dddd:",&x,&y,&z,&t))
			error(" getFields(x,y,z,t) - parameters are needed");
		else
			cpp_fields->getElectricMagneticField(x,y,z,t,Ex,Ey,Ez,Bx,By,Bz);
		
		return Py_BuildValue("dddddd",Ex,Ey,Ez,Bx,By,Bz);
  }
  
  //-----------------------------------------------------
  //destructor for python FieldSourceContainer class (__del__ method).
  //-----------------------------------------------------
  static void FieldSourceContainer_del(pyORBIT_Object* self){
		//std::cerr<<"The FieldSourceContainer __del__ has been called!"<<std::endl;
		delete ((FieldSourceContainer*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python FieldSourceContainer wrapper class
	// they will be vailable from python level
  static PyMethodDef FieldSourceContainerClassMethods[] = {
		{ "AddFieldSource",				 FieldSourceContainer_AddFieldSource,    	METH_VARARGS,"Adds the field source to the container."},
		{ "getFields",				 		FieldSourceContainer_getFields,    		METH_VARARGS,"Gets superposition of fields in the container."},
    {NULL}
  };
	
	// defenition of the memebers of the python FieldSourceContainer wrapper class
	// they will be vailable from python level
	static PyMemberDef FieldSourceContainerClassMembers [] = {
		{NULL}
	};
	
	//new python FieldSourceContainer wrapper type definition
	static PyTypeObject pyORBIT_FieldSourceContainer_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"FieldSourceContainer", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) FieldSourceContainer_del , /*tp_dealloc*/
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
		"The FieldSourceContainer python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		FieldSourceContainerClassMethods, /* tp_methods */
		FieldSourceContainerClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) FieldSourceContainer_init, /* tp_init */
		0, /* tp_alloc */
		FieldSourceContainer_new, /* tp_new */
	};	
	
	
	
	//--------------------------------------------------
	//Initialization function of the pyFieldSourceContainer class
	//It will be called from utils wrapper initialization
	//--------------------------------------------------
  void initFieldSourceContainer(PyObject* module){
		if (PyType_Ready(&pyORBIT_FieldSourceContainer_Type) < 0) return;
		Py_INCREF(&pyORBIT_FieldSourceContainer_Type);
		PyModule_AddObject(module, "FieldSourceContainer", (PyObject *)&pyORBIT_FieldSourceContainer_Type);
		
	}

#ifdef __cplusplus
}
#endif


}
