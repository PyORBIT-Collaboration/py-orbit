#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_bunch.hh"

#include <iostream>

#include "BaseAperture.hh"
#include "BaseApertureShape.hh"

namespace wrap_base_aperture{

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	Constructor for python class wrapping c++ BaseAperture instance.
      It never will be called directly.
	*/
	static PyObject* BaseAperture_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int BaseAperture_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){  	
	  self->cpp_obj =  new BaseAperture();
	  ((BaseAperture*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
  /** Sets the pyBaseApertureShape object for inside(...) method */
  static PyObject* BaseAperture_setApertureShape(PyObject *self, PyObject *args){
	  BaseAperture* cpp_BaseAperture = (BaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
	  PyObject* pyBaseApertureShape;
		if(!PyArg_ParseTuple(args,"O:setApertureShape",&pyBaseApertureShape)){
				ORBIT_MPI_Finalize("BaseAperture - setApertureShape(BaseApertureShape) - parameter is needed. Stop.");
		}
		cpp_BaseAperture->setApertureShape((BaseApertureShape*) ((pyORBIT_Object*) pyBaseApertureShape)->cpp_obj);
		Py_INCREF(Py_None);
		return Py_None;	  
	}
	
  /** Returns the pyBaseApertureShape object for inside(...) method */
  static PyObject* BaseAperture_getApertureShape(PyObject *self, PyObject *args){
	  BaseAperture* cpp_BaseAperture = (BaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
	  BaseApertureShape* baseApertureShape = cpp_BaseAperture->getApertureShape();
	  if(baseApertureShape == NULL){
	  	Py_INCREF(Py_None);
	  	return Py_None;
	  }
	  PyObject* pyBaseApertureShape = baseApertureShape->getPyWrapper();
		Py_INCREF(pyBaseApertureShape);
		return pyBaseApertureShape;	  
	}	
  
  /** Performs the collimation tracking of the bunch */
  static PyObject* BaseAperture_checkBunch(PyObject *self, PyObject *args){
	  BaseAperture* cpp_BaseAperture = (BaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		PyObject* pyBunch;
		PyObject* pyLostBunch;
		int nVars = PyTuple_Size(args);
		if(nVars == 2){
			if(!PyArg_ParseTuple(args,"OO:checkBunch",&pyBunch, &pyLostBunch)){
				ORBIT_MPI_Finalize("BaseAperture - checkBunch(Bunch* bunch, Bunch* lostbunch) - parameters are needed. Stop.");
			}
			PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
			if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type) || !PyObject_IsInstance(pyLostBunch,pyORBIT_Bunch_Type)){
				ORBIT_MPI_Finalize("BaseAperture - checkBunch(Bunch* bunch, Bunch* lostbunch) - method needs a Bunch. Stop.");
			}
			Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
			Bunch* cpp_lostbunch = (Bunch*) ((pyORBIT_Object*)pyLostBunch)->cpp_obj;
			cpp_BaseAperture->checkBunch(cpp_bunch, cpp_lostbunch);
		}
		else{
			if(!PyArg_ParseTuple(args,"O:checkBunch",&pyBunch)){
				ORBIT_MPI_Finalize("BaseAperture - checkBunch(Bunch* bunch) - parameter is needed. Stop.");
			}
			PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
			if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
				ORBIT_MPI_Finalize("BaseAperture - checkBunch(Bunch* bunch) - method needs a Bunch. Stop.");
			}
			Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
			cpp_BaseAperture->checkBunch(cpp_bunch, NULL);			
		}
		Py_INCREF(Py_None);
		return Py_None;
  }
  
	// name([name]) - sets or returns the name of the aperture
  static PyObject* BaseAperture_name(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBaseAperture= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyBaseAperture->cpp_obj;		
    const char* name = NULL;
    if(!PyArg_ParseTuple(	args,"|s:name",&name)){
      ORBIT_MPI_Finalize("BaseAperture.name(...) - call should be - name([name]). Stop.");
    }
		if(name != NULL){
      std::string name_str(name);
      cpp_BaseApertureShape->setName(name_str);
		}
		return Py_BuildValue("s",cpp_BaseApertureShape->getName().c_str());
  }  
		
	/** Sets/Returns the position of the element in the lattice */
	static PyObject* BaseAperture_position(PyObject *self, PyObject *args){
		BaseAperture* cpp_BaseAperture = (BaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		double position = 0;
		int nVars = PyTuple_Size(args);
		if(nVars == 1){
			if(!PyArg_ParseTuple(	args,"d:arguments",&position)){
				ORBIT_MPI_Finalize("BaseAperture - setPosition - cannot parse arguments! It should be (position). Stop.");
			}
			cpp_BaseAperture->setPosition(position);
		}
		position = cpp_BaseAperture->getPosition();
		return Py_BuildValue("d",position);
	}
	
	/** Returns the number of lost particles at this aperture across all CPUs */
	static PyObject* BaseAperture_getNumberOfLost(PyObject *self, PyObject *args){
		BaseAperture* cpp_BaseAperture = (BaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		int nLostParts = cpp_BaseAperture->getNumberOfLost();
		return Py_BuildValue("i",nLostParts);
	}	
	
	/** Sets or returns the aperture state 1 - is active, 0 - switched off*/
	static PyObject* BaseAperture_setOnOff(PyObject *self, PyObject *args){
		BaseAperture* cpp_BaseAperture = (BaseAperture*)((pyORBIT_Object*) self)->cpp_obj;
		int nVars = PyTuple_Size(args);
		int isActive = 1;
		if(nVars == 1){
			if(!PyArg_ParseTuple(	args,"i:onOff",&isActive)){
				ORBIT_MPI_Finalize("BaseAperture.onOff(isActive)  - cannot parse arguments! Stop.");
			}
			if(isActive != 1 && isActive != 0){
				ORBIT_MPI_Finalize("BaseAperture.onOff(isActive)  - isActive should be 1 or 0! Stop.");
			}
			cpp_BaseAperture->setOnOff(isActive);
		}
		isActive = cpp_BaseAperture->getOnOff();
		return Py_BuildValue("b",isActive);	
	}	
	
	
  //-----------------------------------------------------
  //destructor for python BaseAperture class (__del__ method).
  //-----------------------------------------------------
  static void BaseAperture_del(pyORBIT_Object* self){
		//std::cerr<<"The BaseAperture __del__ has been called!"<<std::endl;
		delete ((BaseAperture*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// definition of the methods of the python BaseAperture wrapper class
	// they will be vailable from python level
	static PyMethodDef BaseApertureClassMethods[] = {
		{ "setApertureShape", BaseAperture_setApertureShape, METH_VARARGS,"Sets the pyBaseApertureShape object for inside(...) method."},
		{ "getApertureShape", BaseAperture_getApertureShape, METH_VARARGS,"Returns the pyBaseApertureShape object for inside(...) method."},
		{ "checkBunch",				BaseAperture_checkBunch,       METH_VARARGS,"Performs the aperture check of the bunch."},
		{ "name",             BaseAperture_name,             METH_VARARGS,"Sets or returns the name of the aperture."},
		{ "position",			    BaseAperture_position,         METH_VARARGS,"Sets/Returns the position of the element in the lattice."},
		{ "getNumberOfLost",  BaseAperture_getNumberOfLost,  METH_VARARGS,"Returns the number of lost particles at this aperture across all CPUs."},
		{ "onOff",            BaseAperture_setOnOff,         METH_VARARGS,"Sets/Returns the state of aperture 0 or 1."},
   {NULL}
  };

	static PyMemberDef BaseApertureClassMembers [] = {
		{NULL}
	};
	
	
	//new python BaseAperture wrapper type definition
	static PyTypeObject pyORBIT_BaseAperture_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"BaseAperture", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) BaseAperture_del , /*tp_dealloc*/
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
		"The BaseAperture python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		BaseApertureClassMethods, /* tp_methods */
		BaseApertureClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) BaseAperture_init, /* tp_init */
		0, /* tp_alloc */
		BaseAperture_new, /* tp_new */
	};	
	
	//--------------------------------------------------
	//Initialization BaseAperture class
	//--------------------------------------------------

	void initBaseAperture(PyObject* module){
		//check that the BaseAperture wrapper is ready
		if (PyType_Ready(&pyORBIT_BaseAperture_Type) < 0) return;
		Py_INCREF(&pyORBIT_BaseAperture_Type);
		PyModule_AddObject(module, "BaseAperture", (PyObject *)&pyORBIT_BaseAperture_Type);			
	}

#ifdef __cplusplus
}
#endif


}
