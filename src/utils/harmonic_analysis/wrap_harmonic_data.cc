#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_harmonic_data.hh"

#include <iostream>
#include <string>
#include <cfloat>

#include "HarmonicData.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_harmonicdata{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	    Constructor for python class wrapping c++ HarmonicData instance.
      It never will be called directly.
	*/
	static PyObject* HarmonicData_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int HarmonicData_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
    int nVars = PyTuple_Size(args);
		if(nVars == 2){
			int order = -1; 
			PyObject* pyF;
			if(!PyArg_ParseTuple(args,"iO:",&order,&pyF)){
				error("HarmonicData(order,pyFunction) - constructor parameters are needed");
			}
			PyObject* pyORBIT_Function_Type = getOrbitUtilsType("Function");
			if(!PyObject_IsInstance(pyF,pyORBIT_Function_Type)){
				error("HarmonicData(order,pyFunction) - pyFunction parameter is needed");
			}
			Function* f_in = (Function*) ((pyORBIT_Object*) pyF)->cpp_obj;
			self->cpp_obj =  new HarmonicData(order,f_in);
		}
		else if(nVars == 1){
			PyObject* pyHarmData;
			if(!PyArg_ParseTuple(args,"O:",&pyHarmData)){
				error("HarmonicData(harmonicData) - constructor parameter is needed");
			}
			PyObject* pyORBIT_HarmonicData_Type = getOrbitUtilsType("HarmonicData");
			if(!PyObject_IsInstance(pyHarmData,pyORBIT_HarmonicData_Type)){
				error("HarmonicData(harmonicData) - constructor parameter is needed");
			}
			self->cpp_obj =  new HarmonicData( (HarmonicData*) ((pyORBIT_Object *) pyHarmData)->cpp_obj);
		}
		else{
			error("HarmonicData(order,pyFunction) or HarmonicData(harmonicData)- constructor parameters are needed");
		}
		((HarmonicData*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
	/** It will set or return order of the HarmonicData instance */
  static PyObject* HarmonicData_order(PyObject *self, PyObject *args){
	  HarmonicData* cpp_HarmonicData = (HarmonicData*)((pyORBIT_Object*) self)->cpp_obj;
	  int order = -1;
    int nVars = PyTuple_Size(args);
		if(nVars == 0){
			return Py_BuildValue("i",cpp_HarmonicData->getOrder());
		}
		if(!PyArg_ParseTuple(args,"i:",&order)){
			error("pyHarmonicData.order(order) - parameter is needed");
		}
		cpp_HarmonicData->setOrder(order);
		return Py_BuildValue("i",order);
  }
	
 	/** It will return number of point in the initial data y(x) */
  static PyObject* HarmonicData_dataSize(PyObject *self, PyObject *args){
	  HarmonicData* cpp_HarmonicData = (HarmonicData*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("i",cpp_HarmonicData->dataSize());
  } 
  
	/** It will set the new initial data for fitting in HarmonicData instance */
  static PyObject* HarmonicData_setDataFunction(PyObject *self, PyObject *args){
	  HarmonicData* cpp_HarmonicData = (HarmonicData*)((pyORBIT_Object*) self)->cpp_obj;
		PyObject* pyF;
		if(!PyArg_ParseTuple(args,"O:",&pyF)){
				error("HarmonicData.setDataFunction(pyFunction) - pyFunction parameter is needed");
		}
		PyObject* pyORBIT_Function_Type = getOrbitUtilsType("Function");
		if(!PyObject_IsInstance(pyF,pyORBIT_Function_Type)){
			error("HarmonicData.setDataFunction(pyFunction) - pyFunction parameter is needed");
		}
		Function* f_in = (Function*) ((pyORBIT_Object*) pyF)->cpp_obj;
		cpp_HarmonicData->setDataFunction(f_in);
    Py_INCREF(Py_None);
    return Py_None;
  }
  
 	/** It will set or return the parameters of the HarmonicData instance with index=index*/
  static PyObject* HarmonicData_parameter(PyObject *self, PyObject *args){
	  HarmonicData* cpp_HarmonicData = (HarmonicData*)((pyORBIT_Object*) self)->cpp_obj;
		int index = -1;
		double val = 0.;
    int nVars = PyTuple_Size(args);		
		if(nVars == 1){
			if(!PyArg_ParseTuple(args,"i:",&index)){
				error("pyHarmonicData.parameter(index) - a parameter is needed");
			}
			return Py_BuildValue("d",cpp_HarmonicData->getParam(index));		
		}
		if(nVars == 2){
			if(!PyArg_ParseTuple(args,"id:",&index,&val)){
				error("pyHarmonicData.parameter(index,val) - parameters are needed");
			}	
			cpp_HarmonicData->setParam(index,val);
			Py_INCREF(Py_None);
			return Py_None;
		}
		error("pyHarmonicData.parameter(index[,val]) - parameters are needed");
  }
	
 	/** It will return the y data value for given x index */
  static PyObject* HarmonicData_valueY(PyObject *self, PyObject *args){
	  HarmonicData* cpp_HarmonicData = (HarmonicData*)((pyORBIT_Object*) self)->cpp_obj;
		int ind;
		if(!PyArg_ParseTuple(	args,"i:",&ind)){
			error("pyHarmonicData.valueY(indX) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_HarmonicData->valueY(ind));
  }
  
 	/** It will return the error of y data value for given x index */
  static PyObject* HarmonicData_valueErr(PyObject *self, PyObject *args){
	  HarmonicData* cpp_HarmonicData = (HarmonicData*)((pyORBIT_Object*) self)->cpp_obj;
		int ind;
		if(!PyArg_ParseTuple(	args,"i:",&ind)){
			error("pyHarmonicData.valueErr(indX) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_HarmonicData->valueErr(ind));
  }  
  
 	/** It will return the x data value for given x index */
  static PyObject* HarmonicData_valueX(PyObject *self, PyObject *args){
	  HarmonicData* cpp_HarmonicData = (HarmonicData*)((pyORBIT_Object*) self)->cpp_obj;
		int ind;
		if(!PyArg_ParseTuple(	args,"i:",&ind)){
			error("pyHarmonicData.valueX(indX) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_HarmonicData->valueX(ind));
  }  
  
 	/** It will return the fit Y hramonic value for given x value */
  static PyObject* HarmonicData_fitValueY(PyObject *self, PyObject *args){
	  HarmonicData* cpp_HarmonicData = (HarmonicData*)((pyORBIT_Object*) self)->cpp_obj;
		double x;
		if(!PyArg_ParseTuple(	args,"d:",&x)){
			error("pyHarmonicData.fitValueY(x) - parameter is needed");
		}	
		return Py_BuildValue("d",cpp_HarmonicData->fitValueY(x));
  }  
	
 	/** It clean all data and fitted parameters in the HarmonicData */
  static PyObject* HarmonicData_clean(PyObject *self, PyObject *args){
	  HarmonicData* cpp_HarmonicData = (HarmonicData*)((pyORBIT_Object*) self)->cpp_obj;
    cpp_HarmonicData->clean();
    Py_INCREF(Py_None);
    return Py_None;
  }
  
 	/** It calculates the sum of squares of diff between fit and initial data  */
  static PyObject* HarmonicData_sumDiff2(PyObject *self, PyObject *args){
	  HarmonicData* cpp_HarmonicData = (HarmonicData*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_HarmonicData->sumDiff2());
  }  
	
  //-----------------------------------------------------
  //destructor for python HarmonicData class (__del__ method).
  //-----------------------------------------------------
  static void HarmonicData_del(pyORBIT_Object* self){
		//std::cerr<<"The HarmonicData __del__ has been called! order="<< ((HarmonicData*)self->cpp_obj)->getOrder()<<std::endl;
		delete ((HarmonicData*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python HarmonicData wrapper class
	// they will be vailable from python level
  static PyMethodDef HarmonicDataClassMethods[] = {
		{ "order",				  HarmonicData_order,    	  METH_VARARGS,"Sets or returns order of the HarmonicData."},
		{ "dataSize",       HarmonicData_dataSize,    METH_VARARGS,"Returns the number of points in the initial data y(x)."},
		{ "setDataFunction",HarmonicData_setDataFunction,    METH_VARARGS,"Sets the new function with data for fitting."},		
		{ "parameter",	    HarmonicData_parameter,   METH_VARARGS,"Sets or gets the harmonic parameter with index - parameter(index[,val])"},	
		{ "valueY",		 	    HarmonicData_valueY,      METH_VARARGS,"Returns the Y value of the initial data point with indexX."},
		{ "valueErr",		 	  HarmonicData_valueErr,    METH_VARARGS,"Returns the error of Y value of the initial data point with indexX."},
		{ "valueX",		 	    HarmonicData_valueX,      METH_VARARGS,"Returns the X value of the initial data point with indexX."},
		{ "fitValueY",		  HarmonicData_fitValueY,   METH_VARARGS,"Returns the Y value of the fit function."},
		{ "clean",		 	    HarmonicData_clean,    	  METH_VARARGS,"It cleans all data and fitted parameters"},
		{ "sumDiff2",		 	  HarmonicData_sumDiff2,    METH_VARARGS,"Calculates the sum of squares of diff between fit and initial data"},
    {NULL}
  };
	
	// defenition of the memebers of the python HarmonicData wrapper class
	// they will be vailable from python level
	static PyMemberDef HarmonicDataClassMembers [] = {
		{NULL}
	};
	
	//new python HarmonicData wrapper type definition
	static PyTypeObject pyORBIT_HarmonicData_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"HarmonicData", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) HarmonicData_del , /*tp_dealloc*/
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
		"The HarmonicData python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		HarmonicDataClassMethods, /* tp_methods */
		HarmonicDataClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) HarmonicData_init, /* tp_init */
		0, /* tp_alloc */
		HarmonicData_new, /* tp_new */
	};	
	
	//--------------------------------------------------
	//Initialization plynomial of the pyHarmonicData class
	//--------------------------------------------------
  void initHarmonicData(PyObject* module){
		if (PyType_Ready(&pyORBIT_HarmonicData_Type) < 0) return;
		Py_INCREF(&pyORBIT_HarmonicData_Type);
		PyModule_AddObject(module, "HarmonicData", (PyObject *)&pyORBIT_HarmonicData_Type);
	}

#ifdef __cplusplus
}
#endif


}
