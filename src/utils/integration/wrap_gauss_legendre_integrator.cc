#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"

#include <iostream>
#include <string>

#include "GaussLegendreIntegrator.hh"
#include "wrap_gauss_legendre_integrator.hh"
#include "OU_SplineCH.hh"
#include "OU_Function.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_gl_integrator{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	    Constructor for python class wrapping c++ GaussLegendreIntegrator instance.
      It never will be called directly.
	*/
	static PyObject* GaussLegendreIntegrator_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int GaussLegendreIntegrator_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
    //if nVars == 0 call constructor without arguments
    //if nVars == 1 call constructor with nPoints - integration points
		//if nVars == 3 call constructor with nPoints and x_min x_max - integration points and limits
    int nVars = PyTuple_Size(args);		
		if(nVars == 0){
			self->cpp_obj =  new GaussLegendreIntegrator();	
		}
		if(nVars == 1){
			int nPoints = -1;
			if(!PyArg_ParseTuple(args,"i:",&nPoints)){
				error("GaussLegendreIntegrator(nPoints) - parameter is needed");
			}				
			self->cpp_obj =  new GaussLegendreIntegrator(nPoints);	
		}		
		if(nVars == 3){
			int nPoints = -1;
			double x_min = 0.;
			double x_max = 1.;
			if(!PyArg_ParseTuple(args,"idd:",&nPoints,&x_min,&x_max)){
				error("GaussLegendreIntegrator(nPoints,x_min,x_max) - parameters are needed");
			}				
			self->cpp_obj =  new GaussLegendreIntegrator(nPoints,x_min,x_max);	
		}		
		if(self->cpp_obj == NULL){
			error("GaussLegendreIntegrator([nPoints[,x_min,x_max]]) - constructor signature.");
		}
	  ((GaussLegendreIntegrator*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
	/** It will set the number of integration points */
  static PyObject* GaussLegendreIntegrator_setnPoints(PyObject *self, PyObject *args){
	  GaussLegendreIntegrator* cpp_GaussLegendreIntegrator = (GaussLegendreIntegrator*)((pyORBIT_Object*) self)->cpp_obj;
		int nPoints = -1;
		if(!PyArg_ParseTuple(args,"i:",&nPoints)){
			error("gaussLegendreIntegrator.setnPoints(nPoints) - parameter is needed");
		}	
		cpp_GaussLegendreIntegrator->setnPoints(nPoints);
    Py_INCREF(Py_None);
    return Py_None;		
  }
	
 	/** It will set the  limits of integration */
  static PyObject* GaussLegendreIntegrator_setLimits(PyObject *self, PyObject *args){
	  GaussLegendreIntegrator* cpp_GaussLegendreIntegrator = (GaussLegendreIntegrator*)((pyORBIT_Object*) self)->cpp_obj;
		double x_min = 0.;
		double x_max = 1.;
		if(!PyArg_ParseTuple(args,"dd:",&x_min,&x_max)){
			error("gaussLegendreIntegrator.setLimits(x_min,x_max) - parameters are needed");
		}				
		cpp_GaussLegendreIntegrator->setLimits(x_min,x_max);
    Py_INCREF(Py_None);
    return Py_None;		 
  }
	
 	/** It returns the number of integration points */
  static PyObject* GaussLegendreIntegrator_getnPoints(PyObject *self, PyObject *args){
	  GaussLegendreIntegrator* cpp_GaussLegendreIntegrator = (GaussLegendreIntegrator*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("i",cpp_GaussLegendreIntegrator->getnPoints());
  }
	
 	/** It returns the  limits of integration as tuple */
  static PyObject* GaussLegendreIntegrator_getLimits(PyObject *self, PyObject *args){
	  GaussLegendreIntegrator* cpp_GaussLegendreIntegrator = (GaussLegendreIntegrator*)((pyORBIT_Object*) self)->cpp_obj;
		double x_min = cpp_GaussLegendreIntegrator->getPointsAndWeightFunc()->getMinX();
		double x_max = cpp_GaussLegendreIntegrator->getPointsAndWeightFunc()->getMaxX();		
		return Py_BuildValue("(dd)",x_min,x_max);
  }	
	
 	/** It return the list of [x_point,weight] pairs for integration */
  static PyObject* GaussLegendreIntegrator_getPointsAndWeights(PyObject *self, PyObject *args){
	  GaussLegendreIntegrator* cpp_GaussLegendreIntegrator = (GaussLegendreIntegrator*)((pyORBIT_Object*) self)->cpp_obj;
		Function* f = cpp_GaussLegendreIntegrator->getPointsAndWeightFunc();
		int size = f->getSize();
		PyObject* pyRes = PyTuple_New(size);
		double x = 0., y = 0.;
		for(int i = 0; i < size; i++){
			x = f->x(i);
			y = f->y(i);
			if(PyTuple_SetItem(pyRes,i,Py_BuildValue("(dd)",x,y)) != 0){
					error("gaussLegendreIntegrator.getPointsAndWeights()  cannot create a resulting tuple.");
			}			
		}
		return pyRes;
  }		
	
 	/** It return the integral for the Function or SplineCH instance */
  static PyObject* GaussLegendreIntegrator_integral(PyObject *self, PyObject *args){
	  GaussLegendreIntegrator* cpp_GaussLegendreIntegrator = (GaussLegendreIntegrator*)((pyORBIT_Object*) self)->cpp_obj;
		double sum = 0.;
	  PyObject* pyF;
		Function* f = NULL;
		SplineCH* spline = NULL;
		if(!PyArg_ParseTuple(	args,"O:",&pyF)){
			error("gaussLegendreIntegrator.integral(Function F or SplineCH Spline) - parameter is needed");
		}
		else {
			PyObject* pyORBIT_Function_Type = getOrbitUtilsType("Function");
			PyObject* pyORBIT_SplineCH_Type = getOrbitUtilsType("SplineCH");
			if(PyObject_IsInstance(pyF,pyORBIT_Function_Type)){
				f = (Function*) ((pyORBIT_Object*) pyF)->cpp_obj;
				sum = cpp_GaussLegendreIntegrator->integral(f);
				return Py_BuildValue("d",sum);
			}
			if(PyObject_IsInstance(pyF,pyORBIT_SplineCH_Type)){
				spline = (SplineCH*) ((pyORBIT_Object*) pyF)->cpp_obj;
				sum = cpp_GaussLegendreIntegrator->integral(spline);
				return Py_BuildValue("d",sum);
			}			
		}	
		error("gaussLegendreIntegrator.integral(Function F or SplineCH Spline) - parameter is needed");	
    Py_INCREF(Py_None);
    return Py_None;			
  }

  //-----------------------------------------------------
  //destructor for python GaussLegendreIntegrator class (__del__ method).
  //-----------------------------------------------------
  static void GaussLegendreIntegrator_del(pyORBIT_Object* self){
		//std::cerr<<"The GaussLegendreIntegrator __del__ has been called!"<<std::endl;
		delete ((GaussLegendreIntegrator*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python GaussLegendreIntegrator wrapper class
	// they will be vailable from python level
  static PyMethodDef GaussLegendreIntegratorClassMethods[] = {
		{ "setnPoints",			       GaussLegendreIntegrator_setnPoints,    	    METH_VARARGS,"Sets the number of integration points."},
		{ "setLimits",		 	       GaussLegendreIntegrator_setLimits,    	      METH_VARARGS,"Sets the (min,max) limits for the Integrator"},
		{ "getnPoints",				     GaussLegendreIntegrator_getnPoints,          METH_VARARGS,"Returns the number of integration points."},
 		{ "getLimits",				     GaussLegendreIntegrator_getLimits,    	      METH_VARARGS,"Returns a tuple (min,max) with integration limits"},
 		{ "getPointsAndWeights",	 GaussLegendreIntegrator_getPointsAndWeights, METH_VARARGS,"Returns a list with (x,weight) pairs."},
 		{ "integral",				       GaussLegendreIntegrator_integral,    	      METH_VARARGS,"Returns  the integral for Function or SplineCH instance"},
    {NULL}
  };
	
	// defenition of the memebers of the python GaussLegendreIntegrator wrapper class
	// they will be vailable from python level
	static PyMemberDef GaussLegendreIntegratorClassMembers [] = {
		{NULL}
	};
	
	//new python GaussLegendreIntegrator wrapper type definition
	static PyTypeObject pyORBIT_GaussLegendreIntegrator_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"GaussLegendreIntegrator", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) GaussLegendreIntegrator_del , /*tp_dealloc*/
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
		"The GaussLegendreIntegrator python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		GaussLegendreIntegratorClassMethods, /* tp_methods */
		GaussLegendreIntegratorClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) GaussLegendreIntegrator_init, /* tp_init */
		0, /* tp_alloc */
		GaussLegendreIntegrator_new, /* tp_new */
	};	
	
	
	//--------------------------------------------------
	//Initialization function of the pyGaussLegendreIntegrator class
	//--------------------------------------------------
  void initGLIntegrator(PyObject* module){
		if (PyType_Ready(&pyORBIT_GaussLegendreIntegrator_Type) < 0) return;
		Py_INCREF(&pyORBIT_GaussLegendreIntegrator_Type);
		PyModule_AddObject(module, "GaussLegendreIntegrator", (PyObject *)&pyORBIT_GaussLegendreIntegrator_Type);
	}

#ifdef __cplusplus
}
#endif


}
