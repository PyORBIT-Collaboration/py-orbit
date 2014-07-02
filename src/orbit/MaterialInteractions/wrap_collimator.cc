#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_collimator.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "Collimator.hh"

namespace wrap_collimator{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	Constructor for python class wrapping c++ Collimator instance.
      It never will be called directly.
	*/
	static PyObject* Collimator_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int Collimator_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){

	  double length = 0.;  int ma = 0.; double density_fac = 1; int shape = 0.;
	  double a = 0.; double b = 0.;  double c = 0.; double d = 0.; double angle = 0.;
	  double pos = 0.;
	  
	  //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
	  if(!PyArg_ParseTuple(	args,"dididddddd:arguments",&length,&ma,&density_fac,&shape,&a,&b,&c,&d,&angle,&pos)){
		  error("PyBunch - addParticle - cannot parse arguments! It should be (length,ma,density_fac,shape,a,b,c,d,angle,pos)");
	  }
		self->cpp_obj =  new Collimator(length,ma,density_fac,shape,a,b,c,d,angle,pos);
	  ((Collimator*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
  /** Performs the collimation tracking of the bunch */
  static PyObject* Collimator_collimateBunch(PyObject *self, PyObject *args){
	  Collimator* cpp_Collimator = (Collimator*)((pyORBIT_Object*) self)->cpp_obj;
		PyObject* pyBunch;
		PyObject* pyLostBunch;
		if(!PyArg_ParseTuple(args,"OO:collimateBunch",&pyBunch, &pyLostBunch)){
			ORBIT_MPI_Finalize("Collimator - collimateBunch(Bunch* bunch, Bunch* bunch) - parameter are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type) || !PyObject_IsInstance(pyLostBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("Collimator - collimateBunch(Bunch* bunch, Bunch* bunch) - method needs a Bunch.");
		}
	  
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		Bunch* cpp_lostbunch = (Bunch*) ((pyORBIT_Object*)pyLostBunch)->cpp_obj;
		cpp_Collimator->collimateBunch(cpp_bunch, cpp_lostbunch);
		Py_INCREF(Py_None);
		return Py_None;
  }
	
	
	/** Sets the position of the element in the lattice */
	static PyObject* Collimator_setPosition(PyObject *self, PyObject *args){
		Collimator* cpp_Collimator = (Collimator*)((pyORBIT_Object*) self)->cpp_obj;
		double position = 0;
		if(!PyArg_ParseTuple(	args,"d:arguments",&position)){
			error("PyBunch - setPosition - cannot parse arguments! It should be (position)");
		}
		cpp_Collimator->setPosition(position);
		Py_INCREF(Py_None);
		return Py_None;
	}
	
		
	
  //-----------------------------------------------------
  //destructor for python Collimator class (__del__ method).
  //-----------------------------------------------------
  static void Collimator_del(pyORBIT_Object* self){
		//std::cerr<<"The Collimator __del__ has been called!"<<std::endl;
		delete ((Collimator*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// definition of the methods of the python Collimator wrapper class
	// they will be vailable from python level
	static PyMethodDef CollimatorClassMethods[] = {
		{ "collimateBunch",				 Collimator_collimateBunch,    	METH_VARARGS,"Performs the collimation of the bunch."},
		{ "setPosition",				 Collimator_setPosition,        METH_VARARGS,"Sets the start position of the collimator."},
   {NULL}
  };

	static PyMemberDef CollimatorClassMembers [] = {
		{NULL}
	};
	
	
	//new python Collimator wrapper type definition
	static PyTypeObject pyORBIT_Collimator_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Collimator", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Collimator_del , /*tp_dealloc*/
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
		"The Collimator python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		CollimatorClassMethods, /* tp_methods */
		CollimatorClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Collimator_init, /* tp_init */
		0, /* tp_alloc */
		Collimator_new, /* tp_new */
	};	
	
	
	
	//--------------------------------------------------
	//Initialization Collimator of the pyBunchCollimator class
	//--------------------------------------------------

	void initcollimator(){
		//check that the Collimator wrapper is ready
		if (PyType_Ready(&pyORBIT_Collimator_Type) < 0) return;
		Py_INCREF(&pyORBIT_Collimator_Type);
		//create new module
		PyObject* module = Py_InitModule("collimator",CollimatorClassMethods);
		PyModule_AddObject(module, "Collimator", (PyObject *)&pyORBIT_Collimator_Type);			
	}

#ifdef __cplusplus
}
#endif


}
