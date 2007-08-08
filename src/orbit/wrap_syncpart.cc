///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "wrap_syncpart.hh"

#include "pyORBIT_Object.hh"

#include "Bunch.hh"
#include "SyncPart.hh"

namespace wrap_orbit_syncpart{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python SyncParticle class definition
	//---------------------------------------------------------

	//constructor for python class wrapping SyncPart instance
	//It never will be called directly
	static PyObject* SyncPart_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  SyncParticle class
  //this is implementation of the __init__ method
  static int SyncPart_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		/**
    int nArgs = PyTuple_Size(args);
    if(nArgs != 1){
      error("SyncParticle constructor needs pyBunch instance as parameter.");
    }
		pyORBIT_Object* pyBunch = (pyORBIT_Object *) PyTuple_GetItem(args,0);
		Bunch* cpp_bunch = (Bunch*) pyBunch->cpp_obj;
		if(cpp_bunch->getSyncPart()->getPyWrapper() != NULL){
			error("You should not create SyncParticle class instance directly!");
		}
    self->cpp_obj = (void*) cpp_bunch->getSyncPart();
		cpp_bunch->getSyncPart()->setPyWrapper((PyObject *) self);
		Py_INCREF((PyObject *) self);
    return 0;
		*/

	  int nArgs = PyTuple_Size(args);
    if(nArgs != 1){
      error("SyncParticle constructor needs pyBunch instance as parameter.");
    }
		PyObject* pyBunch = PyTuple_GetItem(args,0);
		PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
		Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

		if(cpp_bunch->getSyncPart()->getPyWrapper() != NULL){
			error("You should not create SyncParticle class instance directly!");
		}

    self->cpp_obj = (void*) cpp_bunch->getSyncPart();
		cpp_bunch->getSyncPart()->setPyWrapper((PyObject *) self);

		Py_DECREF(py_bunch_ref);

		Py_INCREF((PyObject *) self);
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python SyncParticle class.
	//The destructor of cpp_obj is in the bunch class
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static void SyncPart_del(pyORBIT_Object* self){
		self->ob_type->tp_free((PyObject*)self);
  }

 //Sets or returns the kinEnergy for the SyncPart object
  //  the action is depended on the number of arguments
  //  kinEnergy() - returns kinEnergy GeV
  //  kinEnergy(value) - sets the new value for pz and px=0,py=0
  static PyObject* SyncPart_kinEnergy(PyObject *self, PyObject *args){
    //if nVars == 0 this is get kinEnergy
    //if nVars == 1 this is set kinEnergy
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->getEnergy();
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:kinEnergy",&val)){
          error("PySyncPart - kinEnergy(value) - a new value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
				val = cpp_SyncPart->energyToMomentum(val);
        cpp_SyncPart->setPXYZ(0.,0.,val);
				val = cpp_SyncPart->momentumToEnergy(val);
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PySyncPart. You should call kinEnergy() or kinEnergy(value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

  //  beta() - returns beta
  static PyObject* SyncPart_beta(PyObject *self, PyObject *args){
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
		SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
		val = cpp_SyncPart->getBeta();
		return Py_BuildValue("d",val);
  }

  //  gamma() - returns gamma
  static PyObject* SyncPart_gamma(PyObject *self, PyObject *args){
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
		SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
		val = cpp_SyncPart->getGamma();
		return Py_BuildValue("d",val);
  }

  //  mass() - returns mass in GeV
  static PyObject* SyncPart_mass(PyObject *self, PyObject *args){
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
		SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
		val = cpp_SyncPart->getMass();
		return Py_BuildValue("d",val);
  }

 //Sets or returns the momentum for the SyncPart object
  //  the action is depended on the number of arguments
  //  momentum() - returns momentum
  //  momentum(value) - sets the new value for pz and px=0,py=0
  static PyObject* SyncPart_momentum(PyObject *self, PyObject *args){
    //if nVars == 0 this is get momentum
    //if nVars == 1 this is set momentum
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->getMomentum();
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:momentum",&val)){
          error("PySyncPart - momentum(value) - a new value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        cpp_SyncPart->setPXYZ(0.,0.,val);
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PySyncPart. You should call momentum() or momentum(value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }


 //Sets or returns the time in seconds for the SyncPart object
  //  the action is depended on the number of arguments
  //  time() - returns time in seconds
  //  time(value) - sets the new value in seconds
  static PyObject* SyncPart_time(PyObject *self, PyObject *args){
    //if nVars == 0 this is get time
    //if nVars == 1 this is set time
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->getTime();
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:time",&val)){
          error("PySyncPart - time(value) - pySyncPart object and new value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        cpp_SyncPart->setTime(val);
      }
    }
    else{
      error("PySyncPart. You should call time() or time(value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }


 //Sets or returns the x for the SyncPart object
  //  the action is depended on the number of arguments
  //  x() - returns x
  //  x(value) - sets the new value
  static PyObject* SyncPart_x(PyObject *self, PyObject *args){
    //if nVars == 0 this is get x
    //if nVars == 1 this is set x
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->getX();
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:x",&val)){
          error("PySyncPart - x(value) - pySyncPart object and new value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        cpp_SyncPart->setX(val);
      }
    }
    else{
      error("PySyncPart. You should call x() or x(value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

 //Sets or returns the y for the SyncPart object
  //  the action is depended on the number of arguments
  //  y() - returns y
  //  y(value) - sets the new value
  static PyObject* SyncPart_y(PyObject *self, PyObject *args){
    //if nVars == 0 this is get y
    //if nVars == 1 this is set y
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->getY();
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:y",&val)){
          error("PySyncPart - y(value) - a new value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        cpp_SyncPart->setY(val);
      }
    }
    else{
      error("PySyncPart. You should call y() or y(value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

 //Sets or returns the z for the SyncPart object
  //  the action is depended on the number of arguments
  //  z() - returns z
  //  z(value) - sets the new value
  static PyObject* SyncPart_z(PyObject *self, PyObject *args){
    //if nVars == 0 this is get z
    //if nVars == 1 this is set z
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->getZ();
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:z",&val)){
          error("PySyncPart - z(value) - a new value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        cpp_SyncPart->setZ(val);
      }
    }
    else{
      error("PySyncPart. You should call z() or z(value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns the px for the SyncPart object
  //  the action is depended on the number of arguments
  //  px() - returns px
  //  px(value) - sets the new value
  static PyObject* SyncPart_px(PyObject *self, PyObject *args){
    //if nVars == 0 this is get px
    //if nVars == 1 this is set px
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->getPX();
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:px",&val)){
          error("PySyncPart - px(value) - a new value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        cpp_SyncPart->setPX(val);
      }
    }
    else{
      error("PySyncPart. You should call px() or px(value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns the py for the SyncPart object
  //  the action is depended on the number of arguments
  //  py() - returns py
  //  py(value) - sets the new value
  static PyObject* SyncPart_py(PyObject *self, PyObject *args){
    //if nVars == 0 this is get py
    //if nVars == 1 this is set py
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->getPY();
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:py",&val)){
          error("PySyncPart - py(value) - pySyncPart object and new value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        cpp_SyncPart->setPY(val);
      }
    }
    else{
      error("PySyncPart. You should call py() or py(value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns the pz for the SyncPart object
  //  the action is depended on the number of arguments
  //  pz() - returns pz
  //  pz(value) - sets the new value
  static PyObject* SyncPart_pz(PyObject *self, PyObject *args){
    //if nVars == 0 this is get pz
    //if nVars == 1 this is set pz
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->getPZ();
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:pz",&val)){
          error("PySyncPart - pz(value) - pySyncPart object and new value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        cpp_SyncPart->setPZ(val);
      }

    }
    else{
      error("PySyncPart. You should call pz() or pz(value)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Transforms momentum into kinetic energy
  static PyObject* SyncPart_pToE(PyObject *self, PyObject *args){
    //should be nVars == 1
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 1){
        if(!PyArg_ParseTuple(	args,"d:momentumToEnergy",&val)){
          error("PySyncPart - momentumToEnergy(value) - pySyncPart object and input value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->momentumToEnergy(val);
    }
    else{
      error("PySyncPart. You should call momentumToEnergy(p)");
    }
		return Py_BuildValue("d",val);
  }

  //Transforms kinetic energy into momentum
  static PyObject* SyncPart_eToP(PyObject *self, PyObject *args){
    //should be nVars == 1
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pySyncPart = (pyORBIT_Object*) self;
    double val = 0.;
    if(nVars == 1){
        if(!PyArg_ParseTuple(	args,"d:energyToMomentum",&val)){
          error("PySyncPart - energyToMomentum(value) - pySyncPart object and input value are needed");
        }
        SyncPart* cpp_SyncPart = (SyncPart*) pySyncPart->cpp_obj;
        val = cpp_SyncPart->energyToMomentum(val);
    }
    else{
      error("PySyncPart. You should call energyToMomentum(p)");
    }
		return Py_BuildValue("d",val);
  }

	// defenition of the methods of the python SyncPart wrapper class
	// they will be vailable from python level
  static PyMethodDef SyncPartClassMethods[] = {
    { (char *)"mass",             (PyCFunction) SyncPart_mass         ,METH_VARARGS,"Returns mass in GeV"},
    { (char *)"momentum",         (PyCFunction) SyncPart_momentum     ,METH_VARARGS,"Returns or sets momentum in GeV/c. If setting px=0,py=0,pz=val."},
    { (char *)"beta",             (PyCFunction) SyncPart_beta         ,METH_VARARGS,"Returns beta=v/c"},
		{ (char *)"gamma",            (PyCFunction) SyncPart_gamma        ,METH_VARARGS,"Returns gamma=1/sqrt(1-(v/c)**2)"},
		{ (char *)"kinEnergy",        (PyCFunction) SyncPart_kinEnergy    ,METH_VARARGS,"Returns or sets kinetic energy of the synchronous particle in MeV"},
		{ (char *)"time",		          (PyCFunction) SyncPart_time         ,METH_VARARGS,"Sets or returns time in sec"},
		{ (char *)"x",		            (PyCFunction) SyncPart_x            ,METH_VARARGS,"Sets or returns the x-coordinate"},
		{ (char *)"y",		            (PyCFunction) SyncPart_y            ,METH_VARARGS,"Sets or returns the y-coordinate"},
		{ (char *)"z",		            (PyCFunction) SyncPart_z            ,METH_VARARGS,"Sets or returns the z-coordinate"},
		{ (char *)"px",		            (PyCFunction) SyncPart_px           ,METH_VARARGS,"Sets or returns the x-momentum"},
		{ (char *)"py",		            (PyCFunction) SyncPart_py           ,METH_VARARGS,"Sets or returns the y-momentum"},
		{ (char *)"pz",		            (PyCFunction) SyncPart_pz           ,METH_VARARGS,"Sets or returns the z-momentum"},
		{ (char *)"energyToMomentum", (PyCFunction) SyncPart_eToP         ,METH_VARARGS,"Transforms the energy to momentum"},
		{ (char *)"momentumToEnergy", (PyCFunction) SyncPart_pToE         ,METH_VARARGS,"Transforms the momentum to energy"},
    {NULL,NULL}
  };

	// defenition of the memebers of the python SyncPart wrapper class
	// they will be vailable from python level
	static PyMemberDef SyncPartClassMembers [] = {
		{NULL}
	};

	//new python SyncPart wrapper type definition
	static PyTypeObject pyORBIT_SyncPart_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SyncParticle", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SyncPart_del , /*tp_dealloc*/
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
		"The SyncPart python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SyncPartClassMethods, /* tp_methods */
		SyncPartClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SyncPart_init, /* tp_init */
		0, /* tp_alloc */
		SyncPart_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the pySyncPart class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initsyncpart(PyObject* module){
		if (PyType_Ready(&pyORBIT_SyncPart_Type) < 0) return;
		Py_INCREF(&pyORBIT_SyncPart_Type);
		PyModule_AddObject(module, "SyncParticle", (PyObject *)&pyORBIT_SyncPart_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_syncpart
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
