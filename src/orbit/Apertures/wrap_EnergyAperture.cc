#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_aperture.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "EnergyAperture.hh"

namespace wrap_energy_aperture{

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	Constructor for python class wrapping c++ EnergyAperture instance.
      It never will be called directly.
	*/
	static PyObject* EnergyAperture_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int EnergyAperture_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){	  
	  self->cpp_obj =  new EnergyAperture();
	  ((EnergyAperture*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
  /** Performs the collimation tracking of the bunch */
  static PyObject* EnergyAperture_checkBunch(PyObject *self, PyObject *args){
	  EnergyAperture* cpp_EnergyAperture = (EnergyAperture*)((pyORBIT_Object*) self)->cpp_obj;
		PyObject* pyBunch;
		PyObject* pyLostBunch;
		int nVars = PyTuple_Size(args);
		if(nVars == 2){
			if(!PyArg_ParseTuple(args,"OO:checkBunch",&pyBunch, &pyLostBunch)){
				ORBIT_MPI_Finalize("EnergyAperture - checkBunch(Bunch* bunch, Bunch* bunch) - parameters are needed.");
			}
			PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
			if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type) || !PyObject_IsInstance(pyLostBunch,pyORBIT_Bunch_Type)){
				ORBIT_MPI_Finalize("EnergyAperture - checkBunch(Bunch* bunch, Bunch* bunch) - method needs a Bunch.");
			}
			Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
			Bunch* cpp_lostbunch = (Bunch*) ((pyORBIT_Object*)pyLostBunch)->cpp_obj;
			cpp_EnergyAperture->checkBunch(cpp_bunch, cpp_lostbunch);
		}
		else{
			if(!PyArg_ParseTuple(args,"O:checkBunch",&pyBunch)){
				ORBIT_MPI_Finalize("EnergyAperture - checkBunch(Bunch* bunch) - parameter is needed.");
			}
			PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
			if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
				ORBIT_MPI_Finalize("EnergyAperture - checkBunch(Bunch* bunch) - method needs a Bunch.");
			}
			Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
			cpp_EnergyAperture->checkBunch(cpp_bunch, NULL);			
		}
		Py_INCREF(Py_None);
		return Py_None;
  }
		
  /** Sets the min and max energy of the energy aperture class */
	static PyObject* EnergyAperture_setMinMaxEnergy(PyObject *self, PyObject *args){
		EnergyAperture* cpp_EnergyAperture = (EnergyAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double minEnergy = 0.;
		double maxEnergy = 0.;
		if(!PyArg_ParseTuple(	args,"dd:arguments",&minEnergy,&maxEnergy)){
			ORBIT_MPI_Finalize("PyBunch - setMinMaxEnergy - cannot parse arguments! It should be (minEnergy,maxEnergy)");
		}
		cpp_EnergyAperture->setEnergyLimits(minEnergy,maxEnergy);
		Py_INCREF(Py_None);
		return Py_None;
	}  
  
  /** Returns the min and max energy of the energy aperture class */
	static PyObject* EnergyAperture_getMinMaxEnergy(PyObject *self, PyObject *args){
		EnergyAperture* cpp_EnergyAperture = (EnergyAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double minEnergy = cpp_EnergyAperture->getMinEnergy();
		double maxEnergy = cpp_EnergyAperture->getMaxEnergy();
		return Py_BuildValue("(dd)",minEnergy,maxEnergy);
	}  

  /** Returns the position of the element in the lattice */
	static PyObject* EnergyAperture_getPosition(PyObject *self, PyObject *args){
		EnergyAperture* cpp_EnergyAperture = (EnergyAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double position = cpp_EnergyAperture->getPosition();
		return Py_BuildValue("d",position);
	} 	
	
	/** Sets the position of the element in the lattice */
	static PyObject* EnergyAperture_setPosition(PyObject *self, PyObject *args){
		EnergyAperture* cpp_EnergyAperture = (EnergyAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double position = 0;
		if(!PyArg_ParseTuple(	args,"d:arguments",&position)){
			ORBIT_MPI_Finalize("PyBunch - setPosition - cannot parse arguments! It should be (position)");
		}
		cpp_EnergyAperture->setPosition(position);
		Py_INCREF(Py_None);
		return Py_None;
	}
	
  //-----------------------------------------------------
  //destructor for python EnergyAperture class (__del__ method).
  //-----------------------------------------------------
  static void EnergyAperture_del(pyORBIT_Object* self){
		//std::cerr<<"The EnergyAperture __del__ has been called!"<<std::endl;
		delete ((EnergyAperture*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// definition of the methods of the python EnergyAperture wrapper class
	// they will be vailable from python level
	static PyMethodDef EnergyApertureClassMethods[] = {
		{ "checkBunch",				EnergyAperture_checkBunch,    	METH_VARARGS,"Performs the phase aperture check of the bunch."},
		{ "setPosition",			EnergyAperture_setPosition,		  METH_VARARGS,"Sets the position of the element in lattice."},
		{ "getPosition",			EnergyAperture_getPosition,		  METH_VARARGS,"Returns the position of the element in lattice."},
		{ "setMinMaxEnergy",  EnergyAperture_setMinMaxEnergy, METH_VARARGS,"Sets the min and max energy of the energy aperture"},
		{ "getMinMaxEnergy",	EnergyAperture_getMinMaxEnergy, METH_VARARGS,"Returns the min and max energy of the energy aperture"},
   {NULL}
  };

	static PyMemberDef EnergyApertureClassMembers [] = {
		{NULL}
	};
	
	
	//new python EnergyAperture wrapper type definition
	static PyTypeObject pyORBIT_EnergyAperture_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"EnergyAperture", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) EnergyAperture_del , /*tp_dealloc*/
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
		"The EnergyAperture python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		EnergyApertureClassMethods, /* tp_methods */
		EnergyApertureClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) EnergyAperture_init, /* tp_init */
		0, /* tp_alloc */
		EnergyAperture_new, /* tp_new */
	};	
	
	//--------------------------------------------------
	//Initialization EnergyAperture class
	//--------------------------------------------------

	void initEnergyAperture(PyObject* module){
		//check that the EnergyAperture wrapper is ready
		if (PyType_Ready(&pyORBIT_EnergyAperture_Type) < 0) return;
		Py_INCREF(&pyORBIT_EnergyAperture_Type);
		PyModule_AddObject(module, "EnergyAperture", (PyObject *)&pyORBIT_EnergyAperture_Type);			
	}

#ifdef __cplusplus
}
#endif


}
