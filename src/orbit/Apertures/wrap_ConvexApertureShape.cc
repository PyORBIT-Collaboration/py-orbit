#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_bunch.hh"

#include <iostream>

#include "ConvexApertureShape.hh"

namespace wrap_convex_aperture_shape{

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	Constructor for python class wrapping c++ Circle, Ellipse, and Rectangular ApertureShape instances.
	*/
	static PyObject* ConvexApertureShape_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int ConvexApertureShape_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
  	self->cpp_obj = new ConvexApertureShape();
 	  ((BaseApertureShape*) self->cpp_obj)->setPyWrapper((PyObject*) self);	  
    return 0;
  }

	// Sets or returns the center of the shape in X-direction 
  static PyObject* ConvexApertureShape_centerX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyConvexApertureShape = (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyConvexApertureShape->cpp_obj;	
		int nVars = PyTuple_Size(args);
		if(nVars >= 1){
			double center = 0.;
			if(!PyArg_ParseTuple(	args,"d:centerX",&center)){
				ORBIT_MPI_Finalize("ConvexApertureShape.centerX(...) - call should be - centerX([center]). Stop.");
			}
			cpp_BaseApertureShape->setCenterX(center);
    }
		return Py_BuildValue("d",cpp_BaseApertureShape->getCenterX());
  }
  
	// Sets or returns the center of the shape in Y-direction 
  static PyObject* ConvexApertureShape_centerY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyConvexApertureShape= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyConvexApertureShape->cpp_obj;	
		int nVars = PyTuple_Size(args);
		if(nVars >= 1){
			double center = 0.;
			if(!PyArg_ParseTuple(	args,"d:centerY",&center)){
				ORBIT_MPI_Finalize("ConvexApertureShape.centerY(...) - call should be - centerY([center]). Stop.");
			}
			cpp_BaseApertureShape->setCenterY(center);
    }
		return Py_BuildValue("d",cpp_BaseApertureShape->getCenterY());
  }	  

	// name([name]) - sets or returns the name of the shape 
  static PyObject* ConvexApertureShape_name(PyObject *self, PyObject *args){
    pyORBIT_Object* pyConvexApertureShape= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyConvexApertureShape->cpp_obj;		
    const char* name = NULL;
    if(!PyArg_ParseTuple(	args,"|s:name",&name)){
      ORBIT_MPI_Finalize("ConvexApertureShape.name(...) - call should be - name([name]). Stop.");
    }
		if(name != NULL){
      std::string name_str(name);
      cpp_BaseApertureShape->setName(name_str);
		}
		return Py_BuildValue("s",cpp_BaseApertureShape->getName().c_str());
  }
  
	// typeName() - returns the type name of the shape 
  static PyObject* ConvexApertureShape_typeName(PyObject *self, PyObject *args){
    pyORBIT_Object* pyConvexApertureShape= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyConvexApertureShape->cpp_obj;
		int nVars = PyTuple_Size(args);
		if(nVars != 0){
			ORBIT_MPI_Finalize("ConvexApertureShape.typeName() - has no parameters. Stop.");
    }
		return Py_BuildValue("s",cpp_BaseApertureShape->getTypeName().c_str());
  }  
  
	// Sets the convex shape points 
  static PyObject* ConvexApertureShape_setPoints(PyObject *self, PyObject *args){
    pyORBIT_Object* pyConvexApertureShape= (pyORBIT_Object*) self;
		ConvexApertureShape* cpp_ConvexApertureShape = (ConvexApertureShape*) pyConvexApertureShape->cpp_obj;
		PyObject* pyPointList = NULL;
  	if(!PyArg_ParseTuple( args,"O:points",&pyPointList)){
  		ORBIT_MPI_Finalize("ConvexApertureShape.setPoints([[x,y],...] - no right parameters. Stop.");
  	}
  	if(!PyList_CheckExact(pyPointList)){
  		ORBIT_MPI_Finalize("ConvexApertureShape.setPoints([[x,y],...] - no right parameters. Stop.");
  	}
  	cpp_ConvexApertureShape->removeAllPoints();
  	int nPoints = PyList_Size(pyPointList);
  	for(int ind = 0; ind < nPoints; ind++){
  		PyObject* pXY_List = PyList_GetItem(pyPointList,ind);
  		double x = PyFloat_AS_DOUBLE(PyList_GetItem(pXY_List,0));
  		double y = PyFloat_AS_DOUBLE(PyList_GetItem(pXY_List,1));
  		cpp_ConvexApertureShape->addPoint(x,y);	
  	}
		cpp_ConvexApertureShape->checkAllPoints();
    Py_INCREF(Py_None);
		return Py_None;		
  }

	// Returns the convex shape points as array [[x,y],...]
  static PyObject* ConvexApertureShape_getPoints(PyObject *self, PyObject *args){
    pyORBIT_Object* pyConvexApertureShape= (pyORBIT_Object*) self;
		ConvexApertureShape* cpp_ConvexApertureShape = (ConvexApertureShape*) pyConvexApertureShape->cpp_obj;
		int nPoints = cpp_ConvexApertureShape->getPointsX().size();
		PyObject* pyPointList = PyList_New(nPoints);
		for(int ind = 0; ind < nPoints; ind++){
			double x = cpp_ConvexApertureShape->getPointsX()[ind];
			double y = cpp_ConvexApertureShape->getPointsY()[ind];
			PyObject* pyXY_List = PyList_New(2);
			PyList_SET_ITEM(pyXY_List,0,PyFloat_FromDouble(x));
			PyList_SET_ITEM(pyXY_List,1,PyFloat_FromDouble(y));
			PyList_SET_ITEM(pyPointList,ind,pyXY_List);
		}
    Py_INCREF(pyPointList);
		return pyPointList;		
  }  

  //-----------------------------------------------------
  //destructor for python ConvexApertureShape class (__del__ method).
  //-----------------------------------------------------
  static void ConvexApertureShape_del(pyORBIT_Object* self){
		delete ((ConvexApertureShape*) self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// definition of the methods of the python ConvexApertureShape wrapper class
	// they will be vailable from python level
	static PyMethodDef ConvexApertureShapeClassMethods[] = {
		{ "centerX",  ConvexApertureShape_centerX,  METH_VARARGS,"Sets or returns the X-shift of the shape center."},	
		{ "centerY",  ConvexApertureShape_centerY,  METH_VARARGS,"Sets or returns the Y-shift of the shape center."},	
		{ "name",     ConvexApertureShape_name,     METH_VARARGS,"Sets or returns the name of the shape."},
		{ "typeName", ConvexApertureShape_typeName, METH_VARARGS,"Returns the type of the shape."},
		{ "setPoints",ConvexApertureShape_setPoints,METH_VARARGS,"Sets the convex shape points as [[x,y],...]"},
		{ "getPoints",ConvexApertureShape_getPoints,METH_VARARGS,"Returns the convex shape points as [[x,y],...]"},
		{NULL}
  };

	static PyMemberDef ConvexApertureShapeClassMembers [] = {
		{NULL}
	};
	
	
	//new python ConvexApertureShape wrapper type definition
	static PyTypeObject pyORBIT_ConvexApertureShape_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"ConvexApertureShape", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) ConvexApertureShape_del , /*tp_dealloc*/
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
		"The ConvexApertureShape python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		ConvexApertureShapeClassMethods, /* tp_methods */
		ConvexApertureShapeClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) ConvexApertureShape_init, /* tp_init */
		0, /* tp_alloc */
		ConvexApertureShape_new, /* tp_new */
	};	
	
	//--------------------------------------------------
	//Initialization ConvexApertureShape class
	//--------------------------------------------------

	void initConvexApertureShape(PyObject* module){
		//check that the ConvexApertureShape wrapper is ready
		if (PyType_Ready(&pyORBIT_ConvexApertureShape_Type) < 0) return;
		Py_INCREF(&pyORBIT_ConvexApertureShape_Type);
		PyModule_AddObject(module, "ConvexApertureShape", (PyObject *)&pyORBIT_ConvexApertureShape_Type);			
	}

#ifdef __cplusplus
}
#endif


}
