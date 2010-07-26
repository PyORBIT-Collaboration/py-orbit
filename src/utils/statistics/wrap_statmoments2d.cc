#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_statmoments2d.hh"
#include "wrap_mpi_comm.hh"

#include <iostream>
#include <string>

#include "StatMoments2D.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_statmoments2d{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	    Constructor for python class wrapping c++ StatMoments2D instance.
      It never will be called directly.
	*/
	static PyObject* StatMoments2D_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int StatMoments2D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
    //if nVars == 0 constructor without parameters
    //if nVars == 1 constructor with param. - max order
    int nVars = PyTuple_Size(args);		
		if(nVars == 0){
			self->cpp_obj =  new StatMoments2D();
		} else {
			int max_order = -2;
			if(!PyArg_ParseTuple(	args,"i:",&max_order))
				error("pyStatMoments2D(max_order) - constructor - parameter is needed");
			else {
				self->cpp_obj =  new StatMoments2D(max_order);
			}			
		}
	  ((StatMoments2D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
	/** It will set the max order of the moments to the StatMoments2D instance */
  static PyObject* StatMoments2D_setMaxOrder(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
	  int max_order = -2;
		if(!PyArg_ParseTuple(	args,"i:",&max_order))
			error("pyStatMoments2D.setMaxOrder(max_order) - parameter is needed");
		else {
			cpp_StatMoments2D->setMaxOrder(max_order);
		}
		Py_INCREF(Py_None);
		return Py_None;
  }
		
	/** It will account for pair (u,up) inside the StatMoments2D instance */
  static PyObject* StatMoments2D_account(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
	  double u,up;
		if(!PyArg_ParseTuple(	args,"dd:",&u,&up))
			error("pyStatMoments2D.account(u,up) - parameters are needed");
		else {
			cpp_StatMoments2D->account(u,up);
		}
		Py_INCREF(Py_None);
		return Py_None;
  }	

	/** It will return the statistical moments for (u,up) in the StatMoments2D instance */
  static PyObject* StatMoments2D_getStatMoment(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
	  int order_u, order_up;
		double res = 0.;
		if(!PyArg_ParseTuple(	args,"ii:",&order_u,&order_up))
			error("pyStatMoments2D.getStatMoment(order_u,order_up) - parameters are needed");
		else {
			res = cpp_StatMoments2D->getStatMoment(order_u,order_up);
		}
		return Py_BuildValue("d",res);
  }	
	
	/** It will return the statistical moment for u in the StatMoments2D instance*/
  static PyObject* StatMoments2D_getStatMomentU(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
	  int order_u;
		double res = 0.;
		if(!PyArg_ParseTuple(	args,"i:",&order_u))
			error("pyStatMoments2D.getStatMomentU(order_u) - parameter is needed");
		else {
			res = cpp_StatMoments2D->getStatMomentU(order_u);
		}
		return Py_BuildValue("d",res);
  }	
	
	/** It will return the statistical moment for up in the StatMoments2D instance*/
  static PyObject* StatMoments2D_getStatMomentUP(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
	  int order_up;
		double res = 0.;
		if(!PyArg_ParseTuple(	args,"i:",&order_up))
			error("pyStatMoments2D.getStatMomentUP(order_up) - parameter is needed");
		else {
			res = cpp_StatMoments2D->getStatMomentUP(order_up);
		}
		return Py_BuildValue("d",res);
  }	
	
 	/** It returns the minimal value of u in the StatMoments2D instance */
  static PyObject* StatMoments2D_getMinU(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_StatMoments2D->getMinU());
  }	
	
 	/** It returns the minimal value of up in the StatMoments2D instance */
  static PyObject* StatMoments2D_getMinUP(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_StatMoments2D->getMinUP());
  }	
	
 	/** It returns the maximal value of u in the StatMoments2D instance */
  static PyObject* StatMoments2D_getMaxU(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_StatMoments2D->getMaxU());
  }	
	
 	/** It returns the maximal value of up in the StatMoments2D instance */
  static PyObject* StatMoments2D_getMaxUP(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_StatMoments2D->getMaxUP());
  }	
	
 	/** It will return the number of pair (u,up) accounted for in the StatMoments2D instance */
  static PyObject* StatMoments2D_getCount(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("i",cpp_StatMoments2D->getCount());
  }
	
 	/** It will return the maximal order in analysis in the StatMoments2D instance */
  static PyObject* StatMoments2D_getMaxOrder(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("i",cpp_StatMoments2D->getMaxOrder());
  }
	
 	/** It will remove all points and set to zero all variables in the StatMoments2D */
  static PyObject* StatMoments2D_clean(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		cpp_StatMoments2D->clean();
	 	Py_INCREF(Py_None);
		return Py_None; 
  }	
	
	/** It will synchronize the moments through the MPI communicator */ 	
  static PyObject* StatMoments2D_synchronizeMPI(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		int nVars = PyTuple_Size(args);
		if(nVars == 0){
			cpp_StatMoments2D->synchronizeMPI(NULL);
		} else {
			PyObject* py_mpi_comm_type = wrap_orbit_mpi_comm::getMPI_CommType("MPI_Comm");
			PyObject* pyMPIComm = PyTuple_GetItem(args,0);			
			if((!PyObject_IsInstance(pyMPIComm,py_mpi_comm_type))){
				error("StatMoments2D.synchronizeMPI(MPI_Comm) - input parameter is not MPI_Comm");
			}					
			cpp_StatMoments2D->synchronizeMPI((pyORBIT_MPI_Comm*) pyMPIComm);
		}
	 	Py_INCREF(Py_None);
		return Py_None; 
  }		

	//------------------------------------------------------------
	//Twiss functions
	//------------------------------------------------------------
 	/** It returns the emittance */
  static PyObject* StatMoments2D_getEmittance(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_StatMoments2D->getEmittance());
  }		
	
 	/** It returns the Twiss alpha */
  static PyObject* StatMoments2D_getAlpha(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_StatMoments2D->getAlpha());
  }		
	
 	/** It returns the Twiss beta */
  static PyObject* StatMoments2D_getBeta(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_StatMoments2D->getBeta());
  }		
	
 	/** It returns the Twiss gamma */
  static PyObject* StatMoments2D_getGamma(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_StatMoments2D->getGamma());
  }		
	
 	/** It returns the rms value of u */
  static PyObject* StatMoments2D_getRmsU(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_StatMoments2D->getRmsU());
  }		
	
 	/** It returns the rms value of up */
  static PyObject* StatMoments2D_getRmsUP(PyObject *self, PyObject *args){
	  StatMoments2D* cpp_StatMoments2D = (StatMoments2D*)((pyORBIT_Object*) self)->cpp_obj;
		return Py_BuildValue("d",cpp_StatMoments2D->getRmsUP());
  }		
	
	
  //-----------------------------------------------------
  //destructor for python StatMoments2D class (__del__ method).
  //-----------------------------------------------------
  static void StatMoments2D_del(pyORBIT_Object* self){
		//std::cerr<<"The StatMoments2D __del__ has been called!"<<std::endl;
		delete ((StatMoments2D*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python StatMoments2D wrapper class
	// they will be vailable from python level
  static PyMethodDef StatMoments2DClassMethods[] = {
		{ "setMaxOrder",				 StatMoments2D_setMaxOrder,    	METH_VARARGS,"setMaxOrder(max order) sets the maximal order to the StatMoments2D."},
 		{ "getMaxOrder",			 	 StatMoments2D_getMaxOrder,      METH_VARARGS,"Returns the maximal order of statiscis in the StatMoments2D"},		
		{ "account",				     StatMoments2D_account,          METH_VARARGS,"account for (u,up) pair in StatMoments2D"},
 		{ "getStatMoment",			 StatMoments2D_getStatMoment,    METH_VARARGS,"Returns the stat. moment for (order_u, order_up)"},
 		{ "getStatMomentU",			 StatMoments2D_getStatMomentU,   METH_VARARGS,"Returns the stat. moment of order_u for u-variable"},
 		{ "getStatMomentUP",		 StatMoments2D_getStatMomentUP,  METH_VARARGS,"Returns the stat. moment of order_up for up-variable"},
 		{ "getMinU",		 	       StatMoments2D_getMinU,    	    METH_VARARGS,"Returns the minimal U value in the StatMoments2D"},
 		{ "getMaxU",		 	       StatMoments2D_getMaxU,    	    METH_VARARGS,"Returns the maximal U value in the StatMoments2D"},
 		{ "getMinUP",		      	 StatMoments2D_getMinUP,    	    METH_VARARGS,"Returns the minimal UP value in the StatMoments2D"},
 		{ "getMaxUP",		 	       StatMoments2D_getMaxUP,    	    METH_VARARGS,"Returns the maximal UP value in the StatMoments2D"},
 		{ "getCount",			 	     StatMoments2D_getCount,    	    METH_VARARGS,"Returns the number of (u,up) points in the StatMoments2D"},
  	{ "clean",			 	       StatMoments2D_clean,    	      METH_VARARGS,"It will remove all points in the StatMoments2D"},		
  	{ "synchronizeMPI",			 StatMoments2D_synchronizeMPI,   METH_VARARGS,"It will synchronize data over CPUs in specified communicator the StatMoments2D"},
  	{ "getEmittance",			 	 StatMoments2D_getEmittance,    	METH_VARARGS,"Returns the emittance"},		
  	{ "getAlpha",			 	     StatMoments2D_getAlpha,    	    METH_VARARGS,"Returns Twiss alpha"},		
  	{ "getBeta",			 	     StatMoments2D_getBeta,    	    METH_VARARGS,"Returns Twiss beta"},		
  	{ "getGamma",			 	     StatMoments2D_getGamma,    	    METH_VARARGS,"Returns Twiss gamma"},		
  	{ "getRmsU",			 	     StatMoments2D_getRmsU,    	    METH_VARARGS,"Returns the rms value of u"},		
  	{ "getRmsUP",			 	     StatMoments2D_getRmsUP,    	    METH_VARARGS,"Returns the rms value of up"},			
   {NULL}
  };
	
	// defenition of the memebers of the python StatMoments2D wrapper class
	// they will be vailable from python level
	static PyMemberDef StatMoments2DClassMembers [] = {
		{NULL}
	};
	
	//new python StatMoments2D wrapper type definition
	static PyTypeObject pyORBIT_StatMoments2D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"StatMoments2D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) StatMoments2D_del , /*tp_dealloc*/
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
		"The StatMoments2D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		StatMoments2DClassMethods, /* tp_methods */
		StatMoments2DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) StatMoments2D_init, /* tp_init */
		0, /* tp_alloc */
		StatMoments2D_new, /* tp_new */
	};	
	
	
	
	//--------------------------------------------------
	//Initialization StatMoments2D of the pyStatMoments2D class
	//--------------------------------------------------
  void initstatmoments2d(PyObject* module){
		if (PyType_Ready(&pyORBIT_StatMoments2D_Type) < 0) return;
		Py_INCREF(&pyORBIT_StatMoments2D_Type);
		PyModule_AddObject(module, "StatMoments2D", (PyObject *)&pyORBIT_StatMoments2D_Type);
	}

#ifdef __cplusplus
}
#endif


}
