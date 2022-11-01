#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_bunch.hh"

#include <iostream>

#include "CircleApertureShape.hh"
#include "EllipseApertureShape.hh"
#include "RectangularApertureShape.hh"

namespace wrap_py_base_aperture_shape{

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	Constructor for python class wrapping c++ Circle, Ellipse, and Rectangular ApertureShape instances.
	*/
	static PyObject* PrimitiveApertureShape_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int PrimitiveApertureShape_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
  	const char* name = NULL;
  	double param1 =  0.0;
  	double param2 = -1.0;
  	if(!PyArg_ParseTuple(	args,"sd|d:parameters",&name,&param1,&param2)){
  		ORBIT_MPI_Finalize("PrimitiveApertureShape(shapeType, par1[,par2]) - no right parameters. Stop.");
  	}
  	std::string typeName(name);
  	self->cpp_obj = NULL;
  	if(typeName == "circle"){
  		self->cpp_obj =  new CircleApertureShape();
  		((CircleApertureShape*) self->cpp_obj)->setRadius(param1);
  	} 
  	if(typeName == "ellipse"){
  		self->cpp_obj =  new EllipseApertureShape();
  		if(param2 < 0.){
  			ORBIT_MPI_Finalize("PrimitiveApertureShape(ellipse, par1,par2) - par2 is not there or < 0. Stop.");
  		}
  		((EllipseApertureShape*) self->cpp_obj)->setHalfAxisX(param1);
  		((EllipseApertureShape*) self->cpp_obj)->setHalfAxisY(param2);
  	}   	
  	if(typeName == "rectangular"){
  		self->cpp_obj =  new RectangularApertureShape();
  		if(param2 < 0.){
  			ORBIT_MPI_Finalize("PrimitiveApertureShape(rectangular, par1,par2) - par2 is not there or < 0. Stop.");
  		}
  		((RectangularApertureShape*) self->cpp_obj)->setHalfX(param1);
  		((RectangularApertureShape*) self->cpp_obj)->setHalfY(param2);  		
  	}   	
  	if(self->cpp_obj == NULL){
  		ORBIT_MPI_Finalize("PrimitiveApertureShape(shapeType, par1[,par2]) - shapeType should be circle,ellipse, or rectangular. Stop.");
  	}
	  ((BaseApertureShape*) self->cpp_obj)->setPyWrapper((PyObject*) self);	  
    return 0;
  }

	// Sets or returns the center of the shape in X-direction 
  static PyObject* PrimitiveApertureShape_centerX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPrimitiveApertureShape = (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyPrimitiveApertureShape->cpp_obj;	
		int nVars = PyTuple_Size(args);
		if(nVars >= 1){
			double center = 0.;
			if(!PyArg_ParseTuple(	args,"d:centerX",&center)){
				ORBIT_MPI_Finalize("PrimitiveApertureShape.centerX(...) - call should be - centerX([center]). Stop.");
			}
			cpp_BaseApertureShape->setCenterX(center);
    }
		return Py_BuildValue("d",cpp_BaseApertureShape->getCenterX());
  }
  
	// Sets or returns the center of the shape in Y-direction 
  static PyObject* PrimitiveApertureShape_centerY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPrimitiveApertureShape= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyPrimitiveApertureShape->cpp_obj;	
		int nVars = PyTuple_Size(args);
		if(nVars >= 1){
			double center = 0.;
			if(!PyArg_ParseTuple(	args,"d:centerY",&center)){
				ORBIT_MPI_Finalize("PrimitiveApertureShape.centerY(...) - call should be - centerY([center]). Stop.");
			}
			cpp_BaseApertureShape->setCenterY(center);
    }
		return Py_BuildValue("d",cpp_BaseApertureShape->getCenterY());
  }	  

	// name([name]) - sets or returns the name of the shape 
  static PyObject* PrimitiveApertureShape_name(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPrimitiveApertureShape= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyPrimitiveApertureShape->cpp_obj;		
    const char* name = NULL;
    if(!PyArg_ParseTuple(	args,"|s:name",&name)){
      ORBIT_MPI_Finalize("PrimitiveApertureShape.name(...) - call should be - name([name]). Stop.");
    }
		if(name != NULL){
      std::string name_str(name);
      cpp_BaseApertureShape->setName(name_str);
		}
		return Py_BuildValue("s",cpp_BaseApertureShape->getName().c_str());
  }
  
	// typeName() - returns the type name of the shape 
  static PyObject* PrimitiveApertureShape_typeName(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPrimitiveApertureShape= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyPrimitiveApertureShape->cpp_obj;
		int nVars = PyTuple_Size(args);
		if(nVars != 0){
			ORBIT_MPI_Finalize("PrimitiveApertureShape.typeName() - has no parameters. Stop.");
    }
		return Py_BuildValue("s",cpp_BaseApertureShape->getTypeName().c_str());
  }  
  
	// Returns the parameters of the shape like radius of the circle 
  static PyObject* PrimitiveApertureShape_getParamsDict(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPrimitiveApertureShape= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyPrimitiveApertureShape->cpp_obj;	
		PyObject* paramsDict = PyDict_New();
		std::string typeName = cpp_BaseApertureShape->getTypeName();
		double param1 =  0.0;
		double param2 = -1.0;
  	if(typeName == "circle"){
  		param1 = ((CircleApertureShape*) cpp_BaseApertureShape)->getRadius();
  		PyObject* pyRadius = Py_BuildValue("d",param1);
  		PyDict_SetItemString(paramsDict,"radius",pyRadius);
  		Py_DECREF(pyRadius);
  	} 
  	if(typeName == "ellipse"){
  		param1 = ((EllipseApertureShape*) cpp_BaseApertureShape)->getHalfAxisX();
  		param2 = ((EllipseApertureShape*) cpp_BaseApertureShape)->getHalfAxisY();
  		PyObject* pyHalfX = Py_BuildValue("d",param1);
  		PyDict_SetItemString(paramsDict,"halfAxisX",pyHalfX);
  		Py_DECREF(pyHalfX);
  		PyObject* pyHalfY = Py_BuildValue("d",param2);
  		PyDict_SetItemString(paramsDict,"halfAxisY",pyHalfY);
  		Py_DECREF(pyHalfY);
  	}   	
  	if(typeName == "rectangular"){
  		param1 = ((RectangularApertureShape*) cpp_BaseApertureShape)->getHalfX();
  		param2 = ((RectangularApertureShape*) cpp_BaseApertureShape)->getHalfY();
  		PyObject* pyHalfX = Py_BuildValue("d",param1);
  		PyDict_SetItemString(paramsDict,"halfSizeX",pyHalfX);
  		Py_DECREF(pyHalfX);
  		PyObject* pyHalfY = Py_BuildValue("d",param2);
  		PyDict_SetItemString(paramsDict,"halfSizeY",pyHalfY);
  		Py_DECREF(pyHalfY);  		
  	}
		return paramsDict;
  }

	// Sets the parameters of the shape like radius of the circle etc.
  static PyObject* PrimitiveApertureShape_setParams(PyObject *self, PyObject *args){
  	double param1 =  0.0;
  	double param2 = -1.0;
  	if(!PyArg_ParseTuple(	args,"d|d:parameters",&param1,&param2)){
  		ORBIT_MPI_Finalize("PrimitiveApertureShape.setParams(par1[,par2]) - no right parameters. Stop.");
  	}
    pyORBIT_Object* pyPrimitiveApertureShape= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyPrimitiveApertureShape->cpp_obj;	  	
  	std::string typeName = cpp_BaseApertureShape->getTypeName();
  	int settings_done = 0;
  	if(typeName == "circle"){
  		((CircleApertureShape*) cpp_BaseApertureShape)->setRadius(param1);
  		settings_done = 1;
  	} 
  	if(typeName == "ellipse"){
  		if(param2 < 0.){
  			ORBIT_MPI_Finalize("ellipse PrimitiveApertureShape.setParams(par1,par2) - par2 is not there or < 0. Stop.");
  		}
  		((EllipseApertureShape*) cpp_BaseApertureShape)->setHalfAxisX(param1);
  		((EllipseApertureShape*) cpp_BaseApertureShape)->setHalfAxisY(param2);
  		settings_done = 1;
  	}   	
  	if(typeName == "rectangular"){
  		if(param2 < 0.){
  			ORBIT_MPI_Finalize("rectangular PrimitiveApertureShape.setParams(par1,par2) - par2 is not there or < 0. Stop.");
  		}
  		((RectangularApertureShape*) cpp_BaseApertureShape)->setHalfX(param1);
  		((RectangularApertureShape*) cpp_BaseApertureShape)->setHalfY(param2);
  		settings_done = 1;
  	}
  	if(settings_done != 1){
  		ORBIT_MPI_Finalize("PrimitiveApertureShape.setParams(par1,par2) - shape does not defined. Stop.");
  	}
    Py_INCREF(Py_None);
		return Py_None;
  }

  //-----------------------------------------------------
  //destructor for python PrimitiveApertureShape class (__del__ method).
  //-----------------------------------------------------
  static void PrimitiveApertureShape_del(pyORBIT_Object* self){
		//std::cerr<<"debug PrimitiveApertureShape __del__ has been called!"<<std::endl;
		delete ((BaseApertureShape*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// definition of the methods of the python PrimitiveApertureShape wrapper class
	// they will be vailable from python level
	static PyMethodDef PrimitiveApertureShapeClassMethods[] = {
		{ "centerX",      PrimitiveApertureShape_centerX,      METH_VARARGS,"Sets or returns the X-shift of the shape center."},	
		{ "centerY",      PrimitiveApertureShape_centerY,      METH_VARARGS,"Sets or returns the Y-shift of the shape center."},	
		{ "name",         PrimitiveApertureShape_name,         METH_VARARGS,"Sets or returns the name of the shape."},
		{ "typeName",     PrimitiveApertureShape_typeName,     METH_VARARGS,"Returns the type of the shape."},
		{ "getParamsDict",PrimitiveApertureShape_getParamsDict,METH_VARARGS,"Returns params. dict. like radius of the circle."},
		{ "setParams",    PrimitiveApertureShape_setParams,    METH_VARARGS,"Sets params. like radius of the circle, or ellipse half axes."},
		{NULL}
  };

	static PyMemberDef PrimitiveApertureShapeClassMembers [] = {
		{NULL}
	};
	
	
	//new python PrimitiveApertureShape wrapper type definition
	static PyTypeObject pyORBIT_PrimitiveApertureShape_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"PrimitiveApertureShape", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) PrimitiveApertureShape_del , /*tp_dealloc*/
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
		"The PrimitiveApertureShape python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		PrimitiveApertureShapeClassMethods, /* tp_methods */
		PrimitiveApertureShapeClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) PrimitiveApertureShape_init, /* tp_init */
		0, /* tp_alloc */
		PrimitiveApertureShape_new, /* tp_new */
	};	
	
	//--------------------------------------------------
	//Initialization PrimitiveApertureShape class
	//--------------------------------------------------

	void initPrimitiveApertureShape(PyObject* module){
		//check that the PrimitiveApertureShape wrapper is ready
		if (PyType_Ready(&pyORBIT_PrimitiveApertureShape_Type) < 0) return;
		Py_INCREF(&pyORBIT_PrimitiveApertureShape_Type);
		PyModule_AddObject(module, "PrimitiveApertureShape", (PyObject *)&pyORBIT_PrimitiveApertureShape_Type);			
	}

#ifdef __cplusplus
}
#endif


}
