//////////////////////////////// -*- C++ -*- //////////////////////////////
//
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "wrap_bunch.hh"
#include "wrap_syncpart.hh"

#include "pyORBIT_Object.hh"

#include "Bunch.hh"
#include "ParticleAttributesFactory.hh"

namespace wrap_orbit_bunch{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

	//---------------------------------------------------------
	//Python Bunch class definition
	//---------------------------------------------------------

	//constructor for python class wrapping Bunch instance
	//It never will be called directly
	static PyObject* Bunch_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python Bunch class
  //this is implementation of the __init__ method
  static int Bunch_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
    //std::cerr<<"The Bunch __init__ has been called!"<<std::endl;

		//instantiation of a new c++ Bunch
		self->cpp_obj = (void*) new Bunch();
		//This is the way to create new class instance from the C-level
		// Template: PyObject* PyObject_CallMethod(	PyObject *o, char *method, char *format, ...)
		//see Python/C API documentation
		//It will create a SyncParticle object and set the reference to it from pyBunch
		PyObject* mod = PyImport_ImportModule("bunch");
		PyObject* pySyncPart = PyObject_CallMethod(mod,"SyncParticle","O",self);

		//the references should be decreased because they were created as "new reference"
		Py_DECREF(pySyncPart);
		Py_DECREF(mod);
    return 0;		
  }

  //---------------------------------------------------------------
  //
  // methods related to synchronous particle etc.
  //
  //----------------------------------------------------------------

  //returns the SyncPart python class wrapper instance
	static PyObject* Bunch_getSyncParticle(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
		PyObject* pySyncPart = cpp_bunch->getSyncPart()->getPyWrapper();
		Py_INCREF(pySyncPart);
    return pySyncPart;
  }

  //---------------------------------------------------------------
  //
  // add and remove particles, compress etc.
  //
  //----------------------------------------------------------------

  //adds a particle to the Bunch object
  //this is implementation of the  addParticle(...) method
  static PyObject* Bunch_addParticle(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;

    double x = 0.;  double xp = 0.; double y = 0.;
    double yp = 0.; double z = 0.;  double zp = 0.;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"dddddd:coordinates",&x,&xp,&y,&yp,&z,&zp)){
      error("PyBunch - addParticle - cannot parse arguments! It should be (x,xp,y,yp,z,zp)");
    }
    int ind = cpp_bunch->addParticle(x,xp,y,yp,z,zp);
    return Py_BuildValue("i",ind);
  }


  //removes a particle to the Bunch object
  //returns the number of particles in the bunch
  //this is implementation of the deleteParticle(int index)  method
  static PyObject* Bunch_deleteParticle(PyObject *self, PyObject *args){
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    int ind;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"i:deleteParticle",&ind)){
      error("PyBunch - deleteParticle - needs index of particle for deleting");
    }

    cpp_bunch->deleteParticle(ind);
    int size = cpp_bunch->getSize();
		
    return Py_BuildValue("i",size);
  }

  //removes a particle to the Bunch object
  //returns the index of removed macro-particle
  //this is implementation of the deleteParticleFast(int index)  method
  static PyObject* Bunch_deleteParticleFast(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    int ind;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"i:deleteParticleFast",&ind)){
      error("PyBunch - deleteParticleFast - needs index of particle for deleting");
    }

    cpp_bunch->deleteParticleFast(ind);
    return Py_BuildValue("i",ind);
  }

  //removes all particles from the Bunch object
  //this is implementation of the deleteAllParticles()  method
  static PyObject* Bunch_deleteAllParticles(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    cpp_bunch->deleteAllParticles();
    Py_INCREF(Py_None);
    return Py_None;
  }

  //compress the bunch. This method should be called after deleting one
  //  or more macro-particles
  //this is implementation of the compress()  method
  static PyObject* Bunch_compress(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    cpp_bunch->compress();
    Py_INCREF(Py_None);
    return Py_None;
  }

  //---------------------------------------------------------------
  //
  // related to the macro-particles' coordinates
  //
  //----------------------------------------------------------------

  //Sets or returns x coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns x-coordinate
  //  (index, value) - sets the new value to the x-coordinate
  //this is implementation of the x(int index)  method
  static PyObject* Bunch_x(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"i:x",&index)){
          error("PyBunch - x(index) - index is needed");
        }
        val = cpp_bunch->x(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"id:x",&index,&val)){
          error("PyBunch - x(index, value) - index and value are needed");
        }
        cpp_bunch->x(index) = val;
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call x(index) or x(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns y coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns y-coordinate
  //  (index, value) - sets the new value to the y-coordinate
  //this is implementation of the y(int index)  method
  static PyObject* Bunch_y(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"i:y",&index)){
          error("PyBunch - y(index) - index is needed");
        }
        val = cpp_bunch->y(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"id:y",&index,&val)){
          error("PyBunch - y(index, value) - index and value are needed");
        }
        cpp_bunch->y(index) = val;
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call y(index) or y(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns z coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns z(phi)-coordinate
  //  (index, value) - sets the new value to the z(phi)-coordinate
  //this is implementation of the (z or phi)(int index)  method
  static PyObject* Bunch_z(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"i:z",&index)){
          error("PyBunch - z(index) - index is needed");
        }
        val = cpp_bunch->z(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"id:z",&index,&val)){
          error("PyBunch - z(index, value) - index and value are needed");
        }
        cpp_bunch->z(index) = val;
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call z(index) or z(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }


  //Sets or returns px coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns px-coordinate
  //  (index, value) - sets the new value to the px-coordinate
  //this is implementation of the px(int index)  method
  static PyObject* Bunch_px(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"i:px",&index)){
          error("PyBunch - px(index) - index is needed");
        }
        val = cpp_bunch->px(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"id:px",&index,&val)){
          error("PyBunch - px(index, value) - index and value are needed");
        }
        cpp_bunch->px(index) = val;
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call px(index) or px(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns py coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns y-coordinate
  //  (index, value) - sets the new value to the py-coordinate
  //this is implementation of the py(int index)  method
  static PyObject* Bunch_py(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"i:py",&index)){
          error("PyBunch - py(index) - index is needed");
        }
        val = cpp_bunch->py(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"id:py",&index,&val)){
          error("PyBunch - py(index, value) - index and value are needed");
        }
        cpp_bunch->py(index) = val;
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call py(index) or py(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns pz or dE coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns pz(dE)-coordinate
  //  (index, value) - sets the new value to the pz(dE)-coordinate
  //this is implementation of the (pz or dE)(int index)  method
  static PyObject* Bunch_pz(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"i:pz",&index)){
          error("PyBunch - pz(index) - index is needed");
        }
        val = cpp_bunch->pz(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"id:pz",&index,&val)){
          error("PyBunch - pz(index, value) - index and value are needed");
        }
        cpp_bunch->pz(index) = val;
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call pz(index) or pz(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns flag of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns flag
  //  (index, value) - sets the new flag
  //this is implementation of the flag(int index)  method
  static PyObject* Bunch_flag(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 1 get flag
    //if nVars == 2 set flag
    int nVars = PyTuple_Size(args);

    int index = 0;
    int flag = 0;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"i:flag",&index)){
          error("PyBunch - flag(index) - index is needed");
        }
        flag = cpp_bunch->flag(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"ii:flag",&index,&flag)){
          error("PyBunch - flag(index, flag) - index and value are needed");
        }
        cpp_bunch->flag(index) = flag;
      }
			return Py_BuildValue("i",flag);
    }
    else{
      error("PyBunch. You should call flag(index) or flag(index,flag)");
    }

    Py_INCREF(Py_None);
    return Py_None;		
  }

	//Wraps long. coords in the bunch
	//ringwrap(ring_length)
  static PyObject* Bunch_ringwrap(PyObject *self, PyObject *args) {
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;	
    double ring_length = 0.;
		
		//NO NEW OBJECT CREATED BY PyArg_ParseTuple! //NO NEED OF Py_DECREF()
		if(!PyArg_ParseTuple(	args,"d:py",&ring_length)){
			error("PyBunch - ringwrap(ring_length) - pyBunch object needed");
		}
		
		cpp_bunch->ringwrap(ring_length);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //---------------------------------------------------------------
  //
  // related to the bunch predefined attributes
  //
  //----------------------------------------------------------------

  //Sets or returns mass of the macro-particle in MeV
  //  the action is depended on the number of arguments
  //  mass() - returns mass
  //  mass(value) - sets the new value
  //this is implementation of the getMass() and setMass  methods of the Bunch class
  static PyObject* Bunch_mass(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 0 get mass
    //if nVars == 1 set mass
    int nVars = PyTuple_Size(args);

    double val = 0.;

    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        val = cpp_bunch->getMass();
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:mass",&val)){
          error("PyBunch - mass(value) - value is needed");
        }
        cpp_bunch->setMass(val);
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call mass() or mass(value)");
    }

    Py_INCREF(Py_None);
    return Py_None;		
  }

  //Sets or returns classicalRadius of the macro-particle in MeV
  //  the action is depended on the number of arguments
  //  classicalRadius() - returns classicalR
  //  classicalRadius(value) - sets the new value
  //this is implementation of the getClassicalRadius() and
  //setClassicalRadius() methods of the Bunch class
  static PyObject* Bunch_classicalRadius(PyObject *self, PyObject *args){	
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 0 this is get classicalRadius
    //if nVars == 1 this is set classicalRadius
    int nVars = PyTuple_Size(args);

    double val = 0.;

    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        val = cpp_bunch->getClassicalRadius();
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:classicalRadius",&val)){
          error("PyBunch - classicalRadius(value) - value is needed");
        }
        cpp_bunch->setClassicalRadius(val);
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call classicalRadius() or classicalRadius(value)");
    }

    Py_INCREF(Py_None);
    return Py_None;		
  }

  //Sets or returns charge of the macro-particle in e-charge
  //  the action is depended on the number of arguments
  //  charge() - returns charge
  //  charge(value) - sets the new value
  //this is implementation of the getCharge() and setCharge  methods of the Bunch class
  static PyObject* Bunch_charge(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 0 get charge
    //if nVars == 1 set charge
    int nVars = PyTuple_Size(args);
		
    double val = 0.;
		
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        val = cpp_bunch->getCharge();
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:charge",&val)){
          error("PyBunch - charge(value) - value is needed");
        }
        cpp_bunch->setCharge(val);
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call charge() or charge(value)");
    }
		
    Py_INCREF(Py_None);
    return Py_None;			
	}

  //Sets or returns macroSize of the macro-particle
  //  the action is depended on the number of arguments
  //  macroSize() - returns macroSize
  //  macroSize(value) - sets the new value
  //this is implementation of the getMacroSize() and setMacroSize  methods of the Bunch class
  static PyObject* Bunch_macroSize(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 0 get macroSize
    //if nVars == 1 set macroSize
    int nVars = PyTuple_Size(args);

    double val = 0.;

    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        val = cpp_bunch->getMacroSize();
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"d:macroSize",&val)){
          error("PyBunch - macroSize(value) - value is needed");
        }
        cpp_bunch->setMacroSize(val);
      }
			return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call macroSize() or macroSize(value)");
    }

    Py_INCREF(Py_None);
    return Py_None;			
  }

  //---------------------------------------------------------------
  //
  // related to the bunch attributes
  //
  //----------------------------------------------------------------

  //initilizes bunch attributes from the bunch file
  static PyObject* Bunch_initBunchAttr(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* file_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"s:initBunchAttr",&file_name)){
      error("PyBunch - initBunchAttr(fileName) - the file name are needed");
    }
    cpp_bunch->initBunchAttributes(file_name);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns a double bunch attribute
  //  the action is depended on the number of arguments
  //  (attr_name) - returns double-value
  //  (index, value) - sets the new value to the attribute
  //this is implementation of
  // getBunchAttributeDouble(name)
  // setBunchAttributeDouble(name,value) Bunch methods
  static PyObject* Bunch_bunchAttrDouble(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    //if nVars == 2 this is get attribute
    //if nVars == 3 this is set attribute
    int nVars = PyTuple_Size(args);

    const char* attr_name = NULL;
    double val = 0.;

    if(nVars == 2 ||  nVars == 3){
      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"s:bunchAttrDouble",&attr_name)){
          error("PyBunch - bunchAttrDouble(name) - name are needed");
        }
        std::string attr_name_str(attr_name);
        val = cpp_bunch->getBunchAttributeDouble(attr_name_str);
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"sd:bunchAttrDouble",&attr_name,&val)){
          error("PyBunch - bunchAttrDouble(name,value) - name and double value are needed");
        }

        std::string attr_name_str(attr_name);
        cpp_bunch->setBunchAttribute( attr_name_str, val);
				return Py_BuildValue("d",val);
      }
    }
    else{
      error("PyBunch. You should call bunchAttrDouble(name) or bunchAttrDouble(name,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns a integer bunch attribute
  //  the action is depended on the number of arguments
  //  (attr_name) - returns int-value
  //  (index, value) - sets the new value to the attribute
  //this is implementation of
  // getBunchAttributeInt(name)
  // setBunchAttributeInt(name,value) Bunch methods
  static PyObject* Bunch_bunchAttrInt(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 2 this is get attribute
    //if nVars == 3 this is set attribute
    int nVars = PyTuple_Size(args);

    const char* attr_name = NULL;
    int val = 0;

    if(nVars == 2 ||  nVars == 3){
      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"s:bunchAttrInt",&attr_name)){
          error("PyBunch - bunchAttrInt(name) - pyBunch object and name are needed");
        }
        std::string attr_name_str(attr_name);
        val = cpp_bunch->getBunchAttributeInt(attr_name_str);
        return Py_BuildValue("i",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"si:bunchAttrInt",&attr_name,&val)){
          error("PyBunch - bunchAttrInt(name,value) - name, and double value are needed");
        }

        std::string attr_name_str(attr_name);
        cpp_bunch->setBunchAttribute( attr_name_str, val);
				return Py_BuildValue("i",val);
      }
    }
    else{
      error("PyBunch. You should call bunchAttrInt(name) or bunchAttrInt(name,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns a list (tuple) of  ther double bunch attribute names
  static PyObject* Bunch_bunchAttrDoubleNames(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    std::vector<std::string> names;
    cpp_bunch->getDoubleBunchAttributeNames(names);
		//create tuple with names
		PyObject* resTuple = PyTuple_New(names.size());
		for(int i = 0, n = names.size(); i < n; i++){
			PyObject* py_nm = PyString_FromString(names[i].c_str());
			if(PyTuple_SetItem(resTuple,i,py_nm)){
				error("PyBunch - bunchAttrDoubleNames - cannot create tuple with bunch attr names");
			}
		}
    return resTuple;
  }

  //Returns a list (tuple) of  ther integer bunch attribute names
  static PyObject* Bunch_bunchAttrIntNames(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    std::vector<std::string> names;
    cpp_bunch->getIntBunchAttributeNames(names);
		//create tuple with names
		PyObject* resTuple = PyTuple_New(names.size());
		for(int i = 0, n = names.size(); i < n; i++){
			PyObject* py_nm = PyString_FromString(names[i].c_str());
			if(PyTuple_SetItem(resTuple,i,py_nm)){
				error("PyBunch - bunchAttrIntNames() - cannot create tuple with bunch attr names");
			}
		}
    return resTuple;		
  }

  //Returns 0 or 1. The result is 1 if the bunch has an attribute with a particular name
  static PyObject* Bunch_hasBunchAttrDouble(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"s:hasBunchAttrDouble",&attr_name)){
      error("PyBunch - hasBunchAttrDouble(name) - a bunch attr. name are needed");
    }
    std::string attr_name_str(attr_name);
    int res = cpp_bunch->getBunchAttributes()->hasDoubleAttribute(attr_name_str);
    return Py_BuildValue("i",res);
  }

  //Returns 0 or 1. The result is 1 if the bunch has an attribute with a particular name
  static PyObject* Bunch_hasBunchAttrInt(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"s:hasBunchAttrInt",&attr_name)){
      error("PyBunch - hasBunchAttrInt(name) - a bunch attr. name are needed");
    }
    std::string attr_name_str(attr_name);
    int res = cpp_bunch->getBunchAttributes()->hasIntAttribute(attr_name_str);
    return Py_BuildValue("i",res);
  }

  //---------------------------------------------------------------
  //
  // related to particles' attributes
  //
  //----------------------------------------------------------------

  //Adds a particles' attributes with a particular name to the bunch
  static PyObject* Bunch_addPartAttr(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"s:addPartAttr",&attr_name)){
      error("PyBunch - addPartAttr(name) - a particle attr. name are needed");
    }
    std::string attr_name_str(attr_name);
    cpp_bunch->addParticleAttributes(attr_name_str);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Removes a particles' attributes with a particular name from the bunch
  static PyObject* Bunch_removePartAttr(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"s:removePartAttr",&attr_name)){
      error("PyBunch - removePartAttr(name) - pyBunch object and a particle attr. name are needed");
    }
    std::string attr_name_str(attr_name);
    cpp_bunch->removeParticleAttributes(attr_name_str);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Removes all particles' attributes from the bunch
  static PyObject* Bunch_removeAllPartAttr(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    cpp_bunch->removeAllParticleAttributes();
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns a list (tuple) of  the particles' attributes names
  static PyObject* Bunch_getPartAttrNames(PyObject *self, PyObject *args){
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    std::vector<std::string> names;
    cpp_bunch->getParticleAttributesNames(names);
   //create tuple with names
   PyObject* resTuple = PyTuple_New(names.size());
   for(int i = 0, n = names.size(); i < n; i++){
     PyObject* py_nm = PyString_FromString(names[i].c_str());
     if(PyTuple_SetItem(resTuple,i,py_nm)){
       error("PyBunch - getPartAttrNames() - cannot create tuple with bunch attr names");
     }
   }
    return resTuple;
  }

   //Returns a list (tuple) of the possible particles' attributes names
	 static PyObject* Bunch_getPossiblePartAttrNames(PyObject *self, PyObject *args){
		 std::vector<std::string> names;
		 ParticleAttributesFactory::getParticleAttributesNames(names);
		 //create tuple with names
		 PyObject* resTuple = PyTuple_New(names.size());
		 for(int i = 0, n = names.size(); i < n; i++){
			 PyObject* py_nm = PyString_FromString(names[i].c_str());
			 if(PyTuple_SetItem(resTuple,i,py_nm)){
				 error("PyBunch - getPossiblePartAttrNames - cannot create tuple with bunch attr names");
			 }
		 }
		 return resTuple;
	 }

  //temporary removes and memorizes all particles' attributes names
  static PyObject* Bunch_clearAllPartAttrAndMemorize(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    cpp_bunch->clearAllParticleAttributesAndMemorize();
    Py_INCREF(Py_None);
    return Py_None;
  }

  //restores all particles' attributes names from memory
  static PyObject* Bunch_restoreAllPartAttrFromMemory(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    cpp_bunch->restoreAllParticleAttributesFromMemory();
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns 0 or 1. The result is 1 if the bunch has a particles' attributes with a particular name
  static PyObject* Bunch_hasPartAttr(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"s:hasPartAttr",&attr_name)){
      error("PyBunch - hasPartAttr(name) - a particles' attr. name are needed");
    }
    std::string attr_name_str(attr_name);
    int res = cpp_bunch->hasParticleAttributes(attr_name_str);
    return Py_BuildValue("i",res);
  }

  //Returns a list (tuple) of  their bunch particles attribute names specified in the bunch file
  static PyObject* Bunch_readPartAttrNames(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* file_name = NULL;
    std::vector<std::string> names;
    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"s:readPartAttrNames",&file_name)){
      error("PyBunch - readPartAttrNames(fileName) - a file name are needed");
    }
    cpp_bunch->readParticleAttributesNames(file_name,names);
		//create tuple with names
		PyObject* resTuple = PyTuple_New(names.size());
		for(int i = 0, n = names.size(); i < n; i++){
			PyObject* py_nm = PyString_FromString(names[i].c_str());
			if(PyTuple_SetItem(resTuple,i,py_nm)){
				error("PyBunch - readPartAttrNames(fileName) - cannot create tuple with particles attr. names");
			}
		}
    return resTuple;
  }

  //initilizes particles' attributes from the bunch file
  static PyObject* Bunch_initPartAttr(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* file_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"s:initPartAttr",&file_name)){
      error("PyBunch - initPartAttr(fileName) - pyBunch object and file name are needed");
    }
    cpp_bunch->initParticleAttributes(file_name);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns the number of variables in the particles' attributes with a particular name
  static PyObject* Bunch_getPartAttrSize(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"s:getPartAttrSize",&attr_name)){
      error("PyBunch - getPartAttrSize(name) - a particles' attr. name are needed");
    }
    std::string attr_name_str(attr_name);
    int size = cpp_bunch->getParticleAttributes(attr_name_str)->getAttSize();
    return Py_BuildValue("i",size);
  }


  //Sets or returns a particles' attributes' value
  //  the action is depended on the number of arguments
  //  (attr_name,part_index, attr_index) - returns double-value
  //  (attr_name, part_index, attr_index, value) - sets the new value to the attribute
  //This is slow. In the C++ code you have to get reference to
  //particles' attributes object and operate through it
  static PyObject* Bunch_partAttrValue(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 3 this is get attribute
    //if nVars == 4 this is set attribute
    int nVars = PyTuple_Size(args);

    const char* attr_name = NULL;
    int part_index = 0;
    int attr_index = 0;
    double val = 0.;

    if(nVars == 4 ||  nVars == 5){
      if(nVars == 4){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"sii:partAttrValue",&attr_name,&part_index ,&attr_index)){
          error("PyBunch - partAttrValue(attr_name,part_index,atr_index) - params. are needed");
        }
        std::string attr_name_str(attr_name);
        val = cpp_bunch->getParticleAttributes(attr_name_str)->attValue(part_index,attr_index);
        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"siid:",&attr_name,&part_index ,&attr_index,&val)){
          error("PyBunch - partAttrValue(attr_name,part_index,atr_index,value) - params. are needed");
        }
        std::string attr_name_str(attr_name);
        cpp_bunch->getParticleAttributes(attr_name_str)->attValue(part_index,attr_index) = val;
				return Py_BuildValue("d",val);
      }
    }
    else{
      error("PyBunch. You should call partAttrValue(attr_name,part_ind,attr_ind) or partAttrValue(attr_name,part_ind,attr_ind,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //---------------------------------------------------------------
  //
  // getSize, getSizeGlobal, getSizeGlobalFromMemory, setTotalCount
  // getCapacity
  //
  //----------------------------------------------------------------

  //returns the number of macro-particles in the bunch
  //this is implementation of the "getSize()" method
  static PyObject* Bunch_getSize(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;		
    return Py_BuildValue("i",cpp_bunch->getSize());
  }

  //returns the number of macro-particles in the bunch in all CPUs
  //this is implementation of the "getSizeGlobal()" method
  static PyObject* Bunch_getSizeGlobal(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;	
    return Py_BuildValue("i",cpp_bunch->getSizeGlobal());
  }

  //returns the number of macro-particles in the bunch in all CPUs
  //    that was calculated in the previous call of getSizeGlobal()
  //this is implementation of the "getSizeGlobalFromMemory()" method
  static PyObject* Bunch_getSizeGlobalFromMemory(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;	
    return Py_BuildValue("i",cpp_bunch->getSizeGlobalFromMemory());
  }

  //returns the number of all macro-particles - alive, dead, new
  //this is implementation of the "getTotalCount()" method
  static PyObject* Bunch_getTotalCount(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    return Py_BuildValue("i",cpp_bunch->getTotalCount());
  }

  //returns the capacity of the bunch-container. It could be changed.
  //this is implementation of the "getCapacity()" method
  static PyObject* Bunch_getCapacity(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    return Py_BuildValue("i",cpp_bunch->getCapacity());
  }

  //---------------------------------------------------------------
  //
  // write into file or print Bunch
  //
  //----------------------------------------------------------------

  //Prints bunch into the std::cout stream
  static PyObject* Bunch_dumpBunch(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 0 dumpBunchs into std::cout
    //if nVars == 1 dumpBunchs into the file
    int nVars = PyTuple_Size(args);
    const char* file_name = NULL;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        cpp_bunch->print(std::cout);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"s:dumpBunch",&file_name)){
          error("PyBunch - dumpBunch(fileName) - a new value are needed");
        }
        cpp_bunch->print(file_name);
      }
    }
    else{
      error("PyBunch. You should call dumpBunch() or dumpBunch(file_name)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Reads bunch info from the file
  static PyObject* Bunch_readBunch(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 1 reads all macro-particles
    //if nVars == 2 reads only specified number of macro-particles
    int nVars = PyTuple_Size(args);
    const char* file_name = NULL;
    int nParts = 0;
    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"s:read",&file_name)){
          error("PyBunch - readBunch(fileName) - a file name are needed");
        }
        cpp_bunch->readBunch(file_name);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"si:read",&file_name,&nParts)){
          error("PyBunch - readBunch(fileName,nParts) - file name, and number of particles are needed");
        }
        cpp_bunch->readBunch(file_name,nParts);
      }
    }
    else{
      error("PyBunch. You should call readBunch(file_name) or readBunch(file_name,nParts)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Copy bunch attrubutes and structure to another bunch
  static PyObject* Bunch_copyEmptyBunchTo(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    PyObject* pyBunch_Target;
		//NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
		if(!PyArg_ParseTuple(	args,"O:copyEmptyBunchTo",&pyBunch_Target)){
			error("PyBunch - copyEmptyBunchTo(pyBunch) - target pyBunch object is needed");
		}
		Bunch* cpp_target_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch_Target)->cpp_obj;
		cpp_bunch->copyEmptyBunchTo(cpp_target_bunch);
    Py_INCREF(Py_None);
    return Py_None;
  }
	
  //Copy bunch all info including particles coordinates and attributes to another bunch
  static PyObject* Bunch_copyBunchTo(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    PyObject* pyBunch_Target;
		//NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
		if(!PyArg_ParseTuple(	args,"O:copyBunchTo",&pyBunch_Target)){
			error("PyBunch - copyBunchTo(pyBunch) - target pyBunch object is needed");
		}
		Bunch* cpp_target_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch_Target)->cpp_obj;
		cpp_bunch->copyBunchTo(cpp_target_bunch);
		Py_INCREF(Py_None);
    return Py_None;
  }

  //Copy particles coordinates from one bunch to another
  static PyObject* Bunch_addParticlesTo(PyObject *self, PyObject *args){
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    PyObject* pyBunch_Target;
		
		//NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
		if(!PyArg_ParseTuple(	args,"O:addParticlesTo",&pyBunch_Target)){
			error("PyBunch - addParticlesTo(pyBunch) - target pyBunch object is needed");
		}
		Bunch* cpp_target_bunch =(Bunch*) ((pyORBIT_Object *) pyBunch_Target)->cpp_obj ;
		cpp_bunch->addParticlesTo(cpp_target_bunch);
		Py_INCREF(Py_None);
    return Py_None;
  }

  //-----------------------------------------------------
  //destructor for python Bunch class
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static void Bunch_del(pyORBIT_Object* self){
		Bunch* cpp_bunch = (Bunch*) self->cpp_obj;
		delete cpp_bunch;
		self->ob_type->tp_free((PyObject*)self);
  }	

  static PyMethodDef BunchClassMethods[] = {
    //--------------------------------------------------------
    // class Bunch wrapper                        START
    //--------------------------------------------------------
    { "getSyncParticle",                Bunch_getSyncParticle               ,METH_VARARGS,"Returns syncParticle class instance"},
    { "addParticle",                    Bunch_addParticle                   ,METH_VARARGS,"Adds a macro-particle to the bunch"},
    { "deleteParticle",                 Bunch_deleteParticle                ,METH_VARARGS,"Removes macro-particle from the bunch and call compress inside"},
    { "deleteParticleFast",             Bunch_deleteParticleFast            ,METH_VARARGS,"Removes macro-particle from the bunch very fast"},
    { "deleteAllParticles",             Bunch_deleteAllParticles            ,METH_VARARGS,"Removes all macro-particles from the bunch"},
    { "compress",                       Bunch_compress                      ,METH_VARARGS,"Compress the bunch"},
    { "x",                              Bunch_x                             ,METH_VARARGS,"Set x(index,value) or get x(index) coordinate"},
    { "y",                              Bunch_y                             ,METH_VARARGS,"Set y(index,value) or get y(index) coordinate"},
    { "z",                              Bunch_z                             ,METH_VARARGS,"Set z(index,value) or get z(index) coordinate"},
    { "px",                             Bunch_px                            ,METH_VARARGS,"Set px(index,value) or get px(index) coordinate"},
    { "py",                             Bunch_py                            ,METH_VARARGS,"Set py(index,value) or get py(index) coordinate"},
    { "pz",                             Bunch_pz                            ,METH_VARARGS,"Set pz(index,value) or get pz(index) coordinate"},
    { "dE",                             Bunch_pz                            ,METH_VARARGS,"Set dE(index,value) or get dE(index) coordinate"},
    { "xp",                             Bunch_px                            ,METH_VARARGS,"Set xp(index,value) or get xp(index) coordinate"},
    { "yp",                             Bunch_py                            ,METH_VARARGS,"Set yp(index,value) or get yp(index) coordinate"},
    { "flag",                           Bunch_flag                          ,METH_VARARGS,"Set flag(index,value) or get flag(index) coordinate"},
    { "ringwrap",                       Bunch_ringwrap                      ,METH_VARARGS,"Perform the ring wrap. Usage: ringwrap(ring_length)"},
    { "mass",                           Bunch_mass                          ,METH_VARARGS,"Set mass(value) or get mass() the mass of particle in MeV"},
    { "classicalRadius",                Bunch_classicalRadius               ,METH_VARARGS,"Set and get a classical radius of particle in [m]"},
    { "charge",                         Bunch_charge                        ,METH_VARARGS,"Set charge(value) or get charge() the charge of particle in e-charge"},
    { "macroSize",                      Bunch_macroSize                     ,METH_VARARGS,"Set macroSize(value) or get macroSize() the charge of particle in e-charge"},
    { "initBunchAttr",                  Bunch_initBunchAttr                 ,METH_VARARGS,"Reads and initilizes bunch attributes from a bunch file"},
    { "bunchAttrDouble",                Bunch_bunchAttrDouble               ,METH_VARARGS,"Returns and sets a double bunch attribute"},
    { "bunchAttrInt",                   Bunch_bunchAttrInt                  ,METH_VARARGS,"Returns and sets an integer bunch attribute"},
    { "bunchAttrDoubleNames",           Bunch_bunchAttrDoubleNames          ,METH_VARARGS,"Returns a list of double bunch attribute names"},
    { "bunchAttrIntNames",              Bunch_bunchAttrIntNames             ,METH_VARARGS,"Returns a list of integer bunch attribute names"},
    { "hasBunchAttrDouble",             Bunch_hasBunchAttrDouble            ,METH_VARARGS,"Returns 1 if there is a double bunch attr. with this name, 0 - otherwise"},
    { "hasBunchAttrInt",                Bunch_hasBunchAttrInt               ,METH_VARARGS,"Returns 1 if there is a int bunch attr. with this name, 0 - otherwise"},
    { "addPartAttr",                    Bunch_addPartAttr                   ,METH_VARARGS,"Adds a particles' attributes to the bunch"},
    { "removePartAttr",                 Bunch_removePartAttr                ,METH_VARARGS,"Removes a particles' attributes from the bunch"},
    { "removeAllPartAttr",              Bunch_removeAllPartAttr             ,METH_VARARGS,"Removes all particles' attributes from the bunch"},
    { "getPartAttrNames",               Bunch_getPartAttrNames              ,METH_VARARGS,"Returns all particles' attributes names in the bunch at this moment"},
    { "getPossiblePartAttrNames",       Bunch_getPossiblePartAttrNames      ,METH_VARARGS,"Returns all possible particles' attributes names"},
    { "clearAllPartAttrAndMemorize",    Bunch_clearAllPartAttrAndMemorize   ,METH_VARARGS,"Temporary removes and memorizes all particles' attributes names"},
    { "restoreAllPartAttrFromMemory",   Bunch_restoreAllPartAttrFromMemory  ,METH_VARARGS,"Restores all particles' attributes names from memory"},
    { "hasPartAttr",                    Bunch_hasPartAttr                   ,METH_VARARGS,"Returns 1 if there is a particles' attr. with this name, 0 - otherwis"},
    { "readPartAttrNames",              Bunch_readPartAttrNames             ,METH_VARARGS,"Returns a tuple with particles' attr. names in the bunch file"},
    { "initPartAttr",                   Bunch_initPartAttr                  ,METH_VARARGS,"Initializes the particles' attr. from the bunch file"},
    { "getPartAttrSize",                Bunch_getPartAttrSize               ,METH_VARARGS,"Returns the number of variables in the particles' attributes with a particular name"},
    { "partAttrValue",                  Bunch_partAttrValue                 ,METH_VARARGS,"Sets or returns a particles' attribute value"},
    { "getSize",                        Bunch_getSize                       ,METH_VARARGS,"Returns number of macro-particles"},
    { "getSizeGlobal",                  Bunch_getSizeGlobal                 ,METH_VARARGS,"Returns number of macro-particles in all CPUs"},
    { "getSizeGlobalFromMemory",        Bunch_getSizeGlobalFromMemory       ,METH_VARARGS,"Returns number of macro-particles in all CPUs from memory"},
    { "getTotalCount",                  Bunch_getTotalCount                 ,METH_VARARGS,"Returns number of all particles - alive,dead,new"},
    { "getCapacity",                    Bunch_getCapacity                   ,METH_VARARGS,"Returns the capacity of the bunch-contaiter"},
    { "dumpBunch",                      Bunch_dumpBunch                     ,METH_VARARGS,"Prints the bunch info into a standart output stream or file"},
    { "readBunch",                      Bunch_readBunch                     ,METH_VARARGS,"Reads the bunch info from a file"},
    { "copyEmptyBunchTo",               Bunch_copyEmptyBunchTo              ,METH_VARARGS,"Copy bunch attrubutes and structure to another bunch"},
    { "copyBunchTo",                    Bunch_copyBunchTo                   ,METH_VARARGS,"Copy bunch all info including particles coordinates and attributes to another bunch"},
    { "addParticlesTo",                 Bunch_addParticlesTo                ,METH_VARARGS,"Copy particles coordinates from one bunch to another"},
    {NULL,NULL}
    //--------------------------------------------------------
    // class Bunch wrapper                        STOP
    //--------------------------------------------------------
  };
	
	// defenition of the memebers of the python Bunch wrapper class
	// they will be vailable from python level
	static PyMemberDef BunchClassMembers [] = {
		{NULL}
	};
	
	//new python Bunch wrapper type definition
	static PyTypeObject pyORBIT_Bunch_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Bunch", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Bunch_del , /*tp_dealloc*/
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
		"The Bunch python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		BunchClassMethods, /* tp_methods */
		BunchClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Bunch_init, /* tp_init */
		0, /* tp_alloc */
		Bunch_new, /* tp_new */
	};	
	
  static PyMethodDef BunchModuleMethods[] = { {NULL,NULL} };


#ifdef __cplusplus
extern "C" {
#endif

  void initbunch(){
		//check that the Bunch wrapper is ready
		if (PyType_Ready(&pyORBIT_Bunch_Type) < 0) return;
		Py_INCREF(&pyORBIT_Bunch_Type);
    //create new module
    PyObject* module = Py_InitModule("bunch",BunchModuleMethods);
		PyModule_AddObject(module, "Bunch", (PyObject *)&pyORBIT_Bunch_Type);			
		//add the SyncParticle python class
		wrap_orbit_syncpart::initsyncpart(module);
  }


#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_bunch
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
