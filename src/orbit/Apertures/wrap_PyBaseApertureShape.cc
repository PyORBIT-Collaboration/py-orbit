#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_bunch.hh"

#include <iostream>

#include "PyBaseApertureShape.hh"

namespace wrap_py_base_aperture_shape{

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	Constructor for python class wrapping c++ PyBaseApertureShape instance.
      It never will be called directly.
	*/
	static PyObject* PyBaseApertureShape_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int PyBaseApertureShape_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  self->cpp_obj =  new PyBaseApertureShape();
	  ((PyBaseApertureShape*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }

	// Sets or returns the center of the shape in X-direction 
  static PyObject* PyBaseApertureShape_centerX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPyBaseApertureShape= (pyORBIT_Object*) self;
		PyBaseApertureShape* cpp_PyBaseApertureShape = (PyBaseApertureShape*) pyPyBaseApertureShape->cpp_obj;	
		int nVars = PyTuple_Size(args);
		if(nVars >= 1){
			double center = 0.;
			if(!PyArg_ParseTuple(	args,"d:centerX",&center)){
				ORBIT_MPI_Finalize("PyBaseApertureShape.centerX(...) - call should be - centerX([center]). Stop.");
			}
			cpp_PyBaseApertureShape->setCenterX(center);
    }
		return Py_BuildValue("d",cpp_PyBaseApertureShape->getCenterX());
  }
  
	// Sets or returns the center of the shape in Y-direction 
  static PyObject* PyBaseApertureShape_centerY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPyBaseApertureShape= (pyORBIT_Object*) self;
		PyBaseApertureShape* cpp_PyBaseApertureShape = (PyBaseApertureShape*) pyPyBaseApertureShape->cpp_obj;	
		int nVars = PyTuple_Size(args);
		if(nVars >= 1){
			double center = 0.;
			if(!PyArg_ParseTuple(	args,"d:centerY",&center)){
				ORBIT_MPI_Finalize("PyBaseApertureShape.centerY(...) - call should be - centerY([center]). Stop.");
			}
			cpp_PyBaseApertureShape->setCenterY(center);
    }
		return Py_BuildValue("d",cpp_PyBaseApertureShape->getCenterY());
  }	  

	// name([name]) - sets or returns the name of the shape 
  static PyObject* PyBaseApertureShape_name(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPyBaseApertureShape= (pyORBIT_Object*) self;
		PyBaseApertureShape* cpp_PyBaseApertureShape = (PyBaseApertureShape*) pyPyBaseApertureShape->cpp_obj;		
    const char* name = NULL;
    if(!PyArg_ParseTuple(	args,"|s:name",&name)){
      ORBIT_MPI_Finalize("PyBaseApertureShape.name(...) - call should be - name([name]). Stop.");
    }
		if(name != NULL){
      std::string name_str(name);
      cpp_PyBaseApertureShape->setName(name_str);
		}
		return Py_BuildValue("s",cpp_PyBaseApertureShape->getName().c_str());
  }
  
	// typeName() - returns the type name of the shape 
  static PyObject* PyBaseApertureShape_typeName(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPyBaseApertureShape= (pyORBIT_Object*) self;
		PyBaseApertureShape* cpp_PyBaseApertureShape = (PyBaseApertureShape*) pyPyBaseApertureShape->cpp_obj;
		int nVars = PyTuple_Size(args);
		if(nVars != 0){
			ORBIT_MPI_Finalize("PyBaseApertureShape.typeName() - has no parameters. Stop.");
    }
		return Py_BuildValue("s",cpp_PyBaseApertureShape->getTypeName().c_str());
  }  
  
	// Returns empty parameters dictionary - could be overloaded in Python subclass
  static PyObject* PyBaseApertureShape_getParamsDict(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPyBaseApertureShape= (pyORBIT_Object*) self;
		PyBaseApertureShape* cpp_PyBaseApertureShape = (PyBaseApertureShape*) pyPyBaseApertureShape->cpp_obj;	
		PyObject* paramsDict = PyDict_New();
		return paramsDict;
  }  
  
  //-----------------------------------------------------
  //destructor for python PyBaseApertureShape class (__del__ method).
  //-----------------------------------------------------
  static void PyBaseApertureShape_del(pyORBIT_Object* self){
		//std::cerr<<"debug PyBaseApertureShape __del__ has been called!"<<std::endl;
		delete ((PyBaseApertureShape*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// definition of the methods of the python PyBaseApertureShape wrapper class
	// they will be vailable from python level
	static PyMethodDef PyBaseApertureShapeClassMethods[] = {
		{ "centerX",      PyBaseApertureShape_centerX,      METH_VARARGS,"Sets or returns the X-shift of the shape center."},	
		{ "centerY",      PyBaseApertureShape_centerY,      METH_VARARGS,"Sets or returns the Y-shift of the shape center."},	
		{ "name",         PyBaseApertureShape_name,         METH_VARARGS,"Sets or returns the name of the shape."},
		{ "typeName",     PyBaseApertureShape_typeName,     METH_VARARGS,"Returns the type of the shape."},
		{ "getParamsDict",PyBaseApertureShape_getParamsDict,METH_VARARGS,"Returns empty parameters dictionary."},			
		{NULL}
  };

	static PyMemberDef PyBaseApertureShapeClassMembers [] = {
		{NULL}
	};
	
	
	//new python PyBaseApertureShape wrapper type definition
	static PyTypeObject pyORBIT_PyBaseApertureShape_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"PyBaseApertureShape", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) PyBaseApertureShape_del , /*tp_dealloc*/
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
		"The PyBaseApertureShape python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		PyBaseApertureShapeClassMethods, /* tp_methods */
		PyBaseApertureShapeClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) PyBaseApertureShape_init, /* tp_init */
		0, /* tp_alloc */
		PyBaseApertureShape_new, /* tp_new */
	};	
	
	//--------------------------------------------------
	//Initialization PyBaseApertureShape class
	//--------------------------------------------------

	void initPyBaseApertureShape(PyObject* module){
		//check that the PyBaseApertureShape wrapper is ready
		if (PyType_Ready(&pyORBIT_PyBaseApertureShape_Type) < 0) return;
		Py_INCREF(&pyORBIT_PyBaseApertureShape_Type);
		PyModule_AddObject(module, "PyBaseApertureShape", (PyObject *)&pyORBIT_PyBaseApertureShape_Type);			
	}

#ifdef __cplusplus
}
#endif


}
