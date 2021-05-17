#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#include "wrap_spacecharge.hh"

#include "wrap_utils.hh"
#include "QuadFieldSource.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_quad_field_source{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python QuadFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping QuadFieldSource instance
	//It never will be called directly
	static PyObject* QuadFieldSource_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  QuadFieldSource class
  //this is implementation of the __init__ method
  static int QuadFieldSource_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new QuadFieldSource();
		((QuadFieldSource*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }

  /** Sets / Returns the length of quad in [m]*/
  static PyObject* QuadFieldSource_length(PyObject *self, PyObject *args){
	  QuadFieldSource* cpp_fieldSource = (QuadFieldSource*)((pyORBIT_Object*) self)->cpp_obj;
	  int nArgs = PyTuple_Size(args);
	  double length;  
	  if(nArgs == 1){
	  	if(!PyArg_ParseTuple(args,"d:length",&length)){
	  		error("QuadFieldSource.length(length) - parameter is needed.");
	  	}
	  	cpp_fieldSource->setLength(length);
	  }
	  length = cpp_fieldSource->getLength();
	  return Py_BuildValue("d",length);
  }	 
  
  /** Sets - Returns field gradient of the quad in [T/m] */
  static PyObject* QuadFieldSource_gradient(PyObject *self, PyObject *args){
  	QuadFieldSource* cpp_fieldSource = (QuadFieldSource*)((pyORBIT_Object*) self)->cpp_obj;
  	int nArgs = PyTuple_Size(args);
  	double gradient;
	  if(nArgs == 1){
	  	if(!PyArg_ParseTuple(args,"d:gardient",&gradient)){
	  		error("QuadFieldSource.gradient(gradient) - parameter is needed.");
	  	}
	  	cpp_fieldSource->setGradient(gradient);
	  }
	  gradient = cpp_fieldSource->getGradient();
	  return Py_BuildValue("d",gradient); 	
  }
  
   /** Sets or returns X,Y,Z axis symmetries */
  static PyObject* QuadFieldSource_getFields(PyObject *self, PyObject *args){
  	QuadFieldSource* cpp_fieldSource = (QuadFieldSource*)((pyORBIT_Object*) self)->cpp_obj;
  	double x,y,z;
  	if(!PyArg_ParseTuple(args,"ddd:getFields",&x,&y,&z)){
  		ORBIT_MPI_Finalize("QuadFieldSource.getFields(x,y,z) - params needed.");
  	}
  	double fe_x; double fe_y; double fe_z;
  	double fm_x; double fm_y; double fm_z;
  	double t = 0.;
  	cpp_fieldSource->getElectricMagneticField(x,y,z,t,fe_x,fe_y,fe_z,fm_x,fm_y,fm_z);
  	return Py_BuildValue("(dddddd)",fe_x,fe_y,fe_z,fm_x,fm_y,fm_z);
  }
  
  /** Sets / Returns the coordinates transformation matrix 4x4 from external to inner system */
  static PyObject* QuadFieldSource_transormfMatrix(PyObject *self, PyObject *args){
	  QuadFieldSource* cpp_fieldSource = (QuadFieldSource*)((pyORBIT_Object*) self)->cpp_obj;
	  int nArgs = PyTuple_Size(args);
	  PyObject* pyMatrix;
	  Matrix* cpp_matrix;
	  if(nArgs == 1){
	  	if(!PyArg_ParseTuple(args,"O:transormfMatrix",&pyMatrix)){
	  		error("QuadFieldSource.transormfMatrix(Matrix) - parameter is needed.");
	  	}
	  	PyObject* pyORBIT_Matrix_Type = getOrbitUtilsType("Matrix");
	  	if(!PyObject_IsInstance(pyMatrix,pyORBIT_Matrix_Type)){
	  		error("QuadFieldSource.transormfMatrix(Matrix) - parameter is not Matrix.");
	  	}
	  	cpp_matrix = (Matrix*) ((pyORBIT_Object*) pyMatrix)->cpp_obj;
	  	if(cpp_matrix->rows() != 4 || cpp_matrix->columns() != 4){
	  		error("QuadFieldSource.transormfMatrix(Matrix) - Matrix is not 4x4.");
	  	}
	  	// the Py_INCREF(pyMatrix) call will be performed inside setCoordsTransformMatrix(...) method 	  	
	  	cpp_fieldSource->setCoordsTransformMatrix(cpp_matrix);
	  	Py_INCREF(Py_None);
	  	return Py_None;			  	
	  }
	  cpp_matrix = cpp_fieldSource->getCoordsTransformMatrix();
	  pyMatrix = (PyObject*) ((pyORBIT_Object*) cpp_matrix->getPyWrapper());
	  if(pyMatrix == NULL){
	  	error("QuadFieldSource.transormfMatrix() - cannot return Matrix 4x4. You have to assign it first.");
	  }
	  Py_INCREF(pyMatrix);
	  return pyMatrix;
  }   
  
  //-----------------------------------------------------
  //destructor for python QuadFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void QuadFieldSource_del(pyORBIT_Object* self){
		delete ((QuadFieldSource*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python QuadFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef QuadFieldSourceClassMethods[] = {
    { "length",   QuadFieldSource_length    ,METH_VARARGS,"Sets or returns length of quad in [m]"},
    { "gradient", QuadFieldSource_gradient  ,METH_VARARGS, "Sets or returns gradient in quad in [T/m]"},   
    { "getFields",QuadFieldSource_getFields ,METH_VARARGS,"Returns E and B fields (Ex,Ey,Ez,Bx,By,Bz) for (x,y,z) point"},
    { "transormfMatrix", QuadFieldSource_transormfMatrix,METH_VARARGS, "Sets or returns the coordinates transformation matrix 4x4 from external to inner system"},    
    {NULL}
  };

	// defenition of the memebers of the python QuadFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef QuadFieldSourceClassMembers [] = {
		{NULL}
	};

	//new python QuadFieldSource wrapper type definition
	static PyTypeObject pyORBIT_QuadFieldSource_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"QuadFieldSource", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) QuadFieldSource_del , /*tp_dealloc*/
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
		"The QuadFieldSource python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		QuadFieldSourceClassMethods, /* tp_methods */
		QuadFieldSourceClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) QuadFieldSource_init, /* tp_init */
		0, /* tp_alloc */
		QuadFieldSource_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the pyQuadFieldSource class
	//It will be called from wrap_field_sources_module
	//--------------------------------------------------
  void initQuadFieldSource(PyObject* module){
		if (PyType_Ready(&pyORBIT_QuadFieldSource_Type) < 0) return;
		Py_INCREF(&pyORBIT_QuadFieldSource_Type);
		PyModule_AddObject(module, const_cast<char*>("QuadFieldSource"), (PyObject *)&pyORBIT_QuadFieldSource_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_quad_field_source
}
