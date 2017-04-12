#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_aperture.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "PhaseAperture.hh"

namespace wrap_phase_aperture{

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	Constructor for python class wrapping c++ PhaseAperture instance.
      It never will be called directly.
	*/
	static PyObject* PhaseAperture_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int PhaseAperture_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){

	  double frequency = 402.5e+6;
	  
	  //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
	  if(!PyArg_ParseTuple(	args,"d:arguments",&frequency)){
	  	ORBIT_MPI_Finalize("PhaseAperture class constructor - cannot parse arguments! It should be (frequency)");
	  }
	  self->cpp_obj =  new PhaseAperture(frequency);
	  ((PhaseAperture*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
  /** Performs the collimation tracking of the bunch */
  static PyObject* PhaseAperture_checkBunch(PyObject *self, PyObject *args){
	  PhaseAperture* cpp_PhaseAperture = (PhaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		PyObject* pyBunch;
		PyObject* pyLostBunch;
		int nVars = PyTuple_Size(args);
		if(nVars == 2){
			if(!PyArg_ParseTuple(args,"OO:checkBunch",&pyBunch, &pyLostBunch)){
				ORBIT_MPI_Finalize("PhaseAperture - checkBunch(Bunch* bunch, Bunch* bunch) - parameters are needed.");
			}
			PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
			if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type) || !PyObject_IsInstance(pyLostBunch,pyORBIT_Bunch_Type)){
				ORBIT_MPI_Finalize("PhaseAperture - checkBunch(Bunch* bunch, Bunch* bunch) - method needs a Bunch.");
			}
			Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
			Bunch* cpp_lostbunch = (Bunch*) ((pyORBIT_Object*)pyLostBunch)->cpp_obj;
			cpp_PhaseAperture->checkBunch(cpp_bunch, cpp_lostbunch);
		}
		else{
			if(!PyArg_ParseTuple(args,"O:checkBunch",&pyBunch)){
				ORBIT_MPI_Finalize("PhaseAperture - checkBunch(Bunch* bunch) - parameter is needed.");
			}
			PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
			if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
				ORBIT_MPI_Finalize("PhaseAperture - checkBunch(Bunch* bunch) - method needs a Bunch.");
			}
			Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
			cpp_PhaseAperture->checkBunch(cpp_bunch, NULL);			
		}
		Py_INCREF(Py_None);
		return Py_None;
  }
		
  /** Sets the min and max phases of the phase aperture class */
	static PyObject* PhaseAperture_setMinMaxPhase(PyObject *self, PyObject *args){
		PhaseAperture* cpp_PhaseAperture = (PhaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double minPhase = 0.;
		double maxPhase = 0.;
		if(!PyArg_ParseTuple(	args,"dd:arguments",&minPhase,&maxPhase)){
			ORBIT_MPI_Finalize("PyBunch - setMinMaxPhase - cannot parse arguments! It should be (minPhase,maxPhase)");
		}
		cpp_PhaseAperture->setPhaseLimits(minPhase,maxPhase);
		Py_INCREF(Py_None);
		return Py_None;
	}  
  
  /** Returns the min and max phases of the phase aperture class */
	static PyObject* PhaseAperture_getMinMaxPhase(PyObject *self, PyObject *args){
		PhaseAperture* cpp_PhaseAperture = (PhaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double minPhase = cpp_PhaseAperture->getMinPhase();
		double maxPhase = cpp_PhaseAperture->getMaxPhase();
		return Py_BuildValue("(dd)",minPhase,maxPhase);
	}  
  
  /** Returns the RF frequency of the phase aperture class */
	static PyObject* PhaseAperture_getRfFrequency(PyObject *self, PyObject *args){
		PhaseAperture* cpp_PhaseAperture = (PhaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double frequency = cpp_PhaseAperture->getRfFrequency();
		return Py_BuildValue("d",frequency);
	} 
  
 	/** Sets the RF frequency of the phase aperture class */
	static PyObject* PhaseAperture_setRfFrequency(PyObject *self, PyObject *args){
		PhaseAperture* cpp_PhaseAperture = (PhaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double frequency = 0.;
		if(!PyArg_ParseTuple(	args,"d:arguments",&frequency)){
			ORBIT_MPI_Finalize("PyBunch - setRfFrequency - cannot parse arguments! It should be (frequency)");
		}
		cpp_PhaseAperture->setRfFrequency(frequency);
		Py_INCREF(Py_None);
		return Py_None;
	} 
  
  /** Returns the position of the element in the lattice */
	static PyObject* PhaseAperture_getPosition(PyObject *self, PyObject *args){
		PhaseAperture* cpp_PhaseAperture = (PhaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double position = cpp_PhaseAperture->getPosition();
		return Py_BuildValue("d",position);
	} 	
	
	/** Sets the position of the element in the lattice */
	static PyObject* PhaseAperture_setPosition(PyObject *self, PyObject *args){
		PhaseAperture* cpp_PhaseAperture = (PhaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double position = 0;
		if(!PyArg_ParseTuple(	args,"d:arguments",&position)){
			ORBIT_MPI_Finalize("PyBunch - setPosition - cannot parse arguments! It should be (position)");
		}
		cpp_PhaseAperture->setPosition(position);
		Py_INCREF(Py_None);
		return Py_None;
	}
	
  //-----------------------------------------------------
  //destructor for python PhaseAperture class (__del__ method).
  //-----------------------------------------------------
  static void PhaseAperture_del(pyORBIT_Object* self){
		//std::cerr<<"The PhaseAperture __del__ has been called!"<<std::endl;
		delete ((PhaseAperture*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// definition of the methods of the python PhaseAperture wrapper class
	// they will be vailable from python level
	static PyMethodDef PhaseApertureClassMethods[] = {
		{ "checkBunch",				PhaseAperture_checkBunch,    	 METH_VARARGS,"Performs the phase aperture check of the bunch."},
		{ "setPosition",			PhaseAperture_setPosition,		 METH_VARARGS,"Sets the position of the element in lattice."},
		{ "getPosition",			PhaseAperture_getPosition,		 METH_VARARGS,"Returns the position of the element in lattice."},
		{ "setMinMaxPhase",  PhaseAperture_setMinMaxPhase, METH_VARARGS,"Sets the min and max phases of the phase aperture"},
		{ "getMinMaxPhase",	PhaseAperture_getMinMaxPhase, METH_VARARGS,"Returns the min and max phases of the phase aperture"},
		{ "getRfFrequency",	  PhaseAperture_getRfFrequency,	 METH_VARARGS,"Returns the RF frequency of the phase aperture"},
		{ "setRfFrequency",		PhaseAperture_setRfFrequency,	 METH_VARARGS,"Sets the RF frequency of the phase aperture"},	
   {NULL}
  };

	static PyMemberDef PhaseApertureClassMembers [] = {
		{NULL}
	};
	
	
	//new python PhaseAperture wrapper type definition
	static PyTypeObject pyORBIT_PhaseAperture_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"PhaseAperture", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) PhaseAperture_del , /*tp_dealloc*/
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
		"The PhaseAperture python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		PhaseApertureClassMethods, /* tp_methods */
		PhaseApertureClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) PhaseAperture_init, /* tp_init */
		0, /* tp_alloc */
		PhaseAperture_new, /* tp_new */
	};	
	
	//--------------------------------------------------
	//Initialization PhaseAperture class
	//--------------------------------------------------

	void initPhaseAperture(PyObject* module){
		//check that the PhaseAperture wrapper is ready
		if (PyType_Ready(&pyORBIT_PhaseAperture_Type) < 0) return;
		Py_INCREF(&pyORBIT_PhaseAperture_Type);
		PyModule_AddObject(module, "PhaseAperture", (PyObject *)&pyORBIT_PhaseAperture_Type);			
	}

#ifdef __cplusplus
}
#endif


}
