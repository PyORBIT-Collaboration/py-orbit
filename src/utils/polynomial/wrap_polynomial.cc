#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_polynomial.hh"

#include <iostream>
#include <string>

#include "OU_Polynomial.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_plynomial{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	    Constructor for python class wrapping c++ Polynomial instance.
      It never will be called directly.
	*/
	static PyObject* Polynomial_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int Polynomial_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  self->cpp_obj =  new Polynomial(0);	
	  ((Polynomial*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    int nVars = PyTuple_Size(args);
		int order = -1; 
		if(nVars == 1){
			if(!PyArg_ParseTuple(args,"i:",&order)){
				error("pyPolynomial.Polynomial([order]) - constructor parameter is needed");
			}			
			((Polynomial*)self->cpp_obj)->setOrder(order);
		}		
    return 0;
  }
  
	/** It will set or return order of the Polynomial instance */
  static PyObject* Polynomial_order(PyObject *self, PyObject *args){
	  Polynomial* cpp_Polynomial = (Polynomial*)((pyORBIT_Object*) self)->cpp_obj;
	  int order = -1;
    int nVars = PyTuple_Size(args);
		if(nVars == 0){
			return Py_BuildValue("i",cpp_Polynomial->getOrder());
		}
		if(!PyArg_ParseTuple(args,"i:",&order)){
			error("pyPolynomial.order(order) - parameter is needed");
		}
		cpp_Polynomial->setOrder(order);
		return Py_BuildValue("i",order);
  }
	
 	/** It will set or return the coefficient of the Polynomial instance with index=index*/
  static PyObject* Polynomial_coefficient(PyObject *self, PyObject *args){
	  Polynomial* cpp_Polynomial = (Polynomial*)((pyORBIT_Object*) self)->cpp_obj;
		int index = -1;
		double val = 0.;
    int nVars = PyTuple_Size(args);		
		if(nVars == 1){
			if(!PyArg_ParseTuple(args,"i:",&index)){
				error("pyPolynomial.coefficient(index) - a parameter is needed");
			}
			return Py_BuildValue("d",cpp_Polynomial->getCoef(index));		
		}
		if(nVars == 2){
			if(!PyArg_ParseTuple(args,"id:",&index,&val)){
				error("pyPolynomial.coefficient(index,val) - parameters are needed");
			}	
			cpp_Polynomial->setCoef(index,val);
			Py_INCREF(Py_None);
			return Py_None;
		}
		error("pyPolynomial.coef(index[,val]) - parameters are needed");
  }
	
 	/** It will return the polynomial value for x */
  static PyObject* Polynomial_value(PyObject *self, PyObject *args){
	  Polynomial* cpp_Polynomial = (Polynomial*)((pyORBIT_Object*) self)->cpp_obj;
		double x;
		if(!PyArg_ParseTuple(	args,"d:",&x)){
			error("pyPolynomial.value(x) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_Polynomial->value(x));
  }
	
 	/** It will put the derivative into the other polynomial */
  static PyObject* Polynomial_derivative(PyObject *self, PyObject *args){
	  Polynomial* cpp_Polynomial = (Polynomial*)((pyORBIT_Object*) self)->cpp_obj;
	  PyObject* pyP;
		Polynomial* p = NULL;
		if(!PyArg_ParseTuple(	args,"O:",&pyP))
			error("pyPolynomial.derivative(polinomial)- parameter is needed");
		else {
			PyObject* pyORBIT_Polynomial_Type = getOrbitUtilsType("Polynomial");
			if(!PyObject_IsInstance(pyP,pyORBIT_Polynomial_Type)){
				error("pyPolynomial.derivative(polinomial)- parameter is needed");
			}			
			p = (Polynomial*) ((pyORBIT_Object*) pyP)->cpp_obj;
			cpp_Polynomial->derivative(p);
			int order = p->getOrder();
			for(int i = 0; i < (order+1); i++){
			}
		}		
		Py_INCREF(pyP);
		return pyP;
  }
	
 	/** It will put the copy into the other polynomial */
  static PyObject* Polynomial_copyTo(PyObject *self, PyObject *args){
	  Polynomial* cpp_Polynomial = (Polynomial*)((pyORBIT_Object*) self)->cpp_obj;
	  PyObject* pyP;
		Polynomial* p = NULL;
		if(!PyArg_ParseTuple(	args,"O:",&pyP))
			error("pyPolynomial.copyTo(polinomial)- parameter is needed");
		else {
			PyObject* pyORBIT_Polynomial_Type = getOrbitUtilsType("Polynomial");
			if(!PyObject_IsInstance(pyP,pyORBIT_Polynomial_Type)){
				error("pyPolynomial.copyTo(polinomial)- parameter is needed");
			}			
			p = (Polynomial*) ((pyORBIT_Object*) pyP)->cpp_obj;
			cpp_Polynomial->copyTo(p);
		}		
		Py_INCREF(pyP);
		return pyP;
  }
	
  //-----------------------------------------------------
  //destructor for python Polynomial class (__del__ method).
  //-----------------------------------------------------
  static void Polynomial_del(pyORBIT_Object* self){
		//std::cerr<<"The Polynomial __del__ has been called!"<<std::endl;
		delete ((Polynomial*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python Polynomial wrapper class
	// they will be vailable from python level
  static PyMethodDef PolynomialClassMethods[] = {
		{ "order",				 Polynomial_order,    	  METH_VARARGS,"Sets or returns order of the Polynomial."},
		{ "coefficient",	 Polynomial_coefficient, METH_VARARGS,"Sets or gets the cofficient with index - coefficient(index[,val])"},
		{ "value",		 	   Polynomial_value,       METH_VARARGS,"Returns the value of the polynomial."},
		{ "derivative",		 Polynomial_derivative,  METH_VARARGS,"derivative(polynomial) - initializes another polynomial as derivative"},
		{ "copyTo",		 	   Polynomial_copyTo,    	METH_VARARGS,"copyTo(polynomial) - initializes another polynomial as a copy"},
    {NULL}
  };
	
	// defenition of the memebers of the python Polynomial wrapper class
	// they will be vailable from python level
	static PyMemberDef PolynomialClassMembers [] = {
		{NULL}
	};
	
	//new python Polynomial wrapper type definition
	static PyTypeObject pyORBIT_Polynomial_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Polynomial", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Polynomial_del , /*tp_dealloc*/
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
		"The Polynomial python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		PolynomialClassMethods, /* tp_methods */
		PolynomialClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Polynomial_init, /* tp_init */
		0, /* tp_alloc */
		Polynomial_new, /* tp_new */
	};	
	
	//--------------------------------------------------
	//Initialization plynomial of the pyPolynomial class
	//--------------------------------------------------
  void initPolynomial(PyObject* module){
		if (PyType_Ready(&pyORBIT_Polynomial_Type) < 0) return;
		Py_INCREF(&pyORBIT_Polynomial_Type);
		PyModule_AddObject(module, "Polynomial", (PyObject *)&pyORBIT_Polynomial_Type);
		
	}

#ifdef __cplusplus
}
#endif


}
