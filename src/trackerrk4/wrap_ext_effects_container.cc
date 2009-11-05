#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_ext_effects_container.hh"

#include <iostream>
#include <string>

#include "ExtEffectsContainer.hh"
#include "ExternalEffects.hh"


using namespace OrbitUtils;


namespace wrap_ext_effects_container{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif


	//constructor for python class wrapping CppExternalEffects instance
	//It never will be called directly
	static PyObject* ExtEffectsContainer_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}


  //this is implementation of the __init__ method
  static int ExtEffectsContainer_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  

	  self->cpp_obj =  new  ExtEffectsContainer();
	  ((ExtEffectsContainer*) self->cpp_obj)->setPyWrapper((PyObject*) self);

	
    return 0;
  }
  
		

  static PyObject* ExtEffectsContainer_AddEffect(PyObject *self, PyObject *args){
	  
	  ExtEffectsContainer* cpp_ExtEffectsContainer = (ExtEffectsContainer*)((pyORBIT_Object*) self)->cpp_obj;
	 
	  ExternalEffects* extEf;
	  PyObject* pyExtEffects;

       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"O:",&pyExtEffects))
             error(" AddEffect(external effect) - parameter is needed");
           else {
           extEf = (ExternalEffects*) ((pyORBIT_Object*) pyExtEffects)->cpp_obj;
           cpp_ExtEffectsContainer->AddEffect(extEf);
           }
           
  		    Py_INCREF(Py_None);
  		    return Py_None;

  }
  
  
  
  
  
	
  //-----------------------------------------------------
  //destructor for python ExtEffectsContainer class (__del__ method).
  //-----------------------------------------------------
  static void ExtEffectsContainer_del(pyORBIT_Object* self){
		//std::cerr<<"The LasStripExternalEffects __del__ has been called!"<<std::endl;
		delete ((ExtEffectsContainer*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python ExtEffectsContainer wrapper class
	// they will be vailable from python level
  static PyMethodDef ExtEffectsContainerClassMethods[] = {
		{ "AddEffect",				 ExtEffectsContainer_AddEffect,    	METH_VARARGS,"Adds external effect to container."},

    {NULL}
  };

	// defenition of the memebers of the python ExtEffectsContainer wrapper class
	// they will be vailable from python level
	static PyMemberDef ExtEffectsContainerClassMembers [] = {
		{NULL}
	};

	//new python ExtEffectsContainer wrapper type definition
	static PyTypeObject pyORBIT_ExtEffectsContainer_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"ExtEffectsContainer", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) ExtEffectsContainer_del , /*tp_dealloc*/
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
		"The ExtEffectsContainer python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		ExtEffectsContainerClassMethods, /* tp_methods */
		ExtEffectsContainerClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) ExtEffectsContainer_init, /* tp_init */
		0, /* tp_alloc */
		ExtEffectsContainer_new, /* tp_new */
	};	


		
	//--------------------------------------------------
	//Initialization function of the pyExtEffectsContainer class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initExtEffectsContainer(PyObject* module){
		if (PyType_Ready(&pyORBIT_ExtEffectsContainer_Type) < 0) return;
		Py_INCREF(&pyORBIT_ExtEffectsContainer_Type);
		PyModule_AddObject(module, "ExtEffectsContainer", (PyObject *)&pyORBIT_ExtEffectsContainer_Type);
				
	}

#ifdef __cplusplus
}
#endif


}
