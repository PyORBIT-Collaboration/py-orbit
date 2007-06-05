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

#include "Bunch.hh"
#include "ParticleAttributesFactory.hh"

namespace wrap_orbit_bunch{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

	//---------------------------------------------------------
	//Python Bunch class definition
	//---------------------------------------------------------

  //constructor for python Bunch class
  //this is implementation of the __init__ method
  static PyObject* Bunch_init(PyObject *self, PyObject *args){
    //std::cerr<<"The Bunch __init__ has been called!"<<std::endl;

    PyObject* pyBunch = PyTuple_GetItem(args,0);
    Bunch* cpp_bunch = new Bunch();
    PyObject* py_bunch_ref = PyCObject_FromVoidPtr((void *) cpp_bunch, NULL);
    if(PyObject_SetAttrString(pyBunch	, "cpp_ptr", py_bunch_ref) < 0){
      error("Bunch wrapper. Can not set c-bunch reference.");
    }
    Py_DECREF(py_bunch_ref);

		//This is the way to create new class instance from the C-level
		// Template: PyObject* PyObject_CallMethod(	PyObject *o, char *method, char *format, ...)
		//see Python/C API documentation
		PyObject* mod = PyImport_ImportModule("bunch");
		PyObject* pySyncPart = PyObject_CallMethod(mod,"SyncParticle","O",pyBunch);

		//the references should be decreased because they were created as "new reference"
		Py_DECREF(pySyncPart);
		Py_DECREF(mod);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //---------------------------------------------------------------
  //
  // methods related to synchronous particle etc.
  //
  //----------------------------------------------------------------

  //returns the SyncPart python class wrapper instance
  static PyObject* Bunch_getSyncParticle(PyObject *self, PyObject *args){
    PyObject* pyBunch;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"O:coordinates",&pyBunch)){
      error("PyBunch - Bunch_getSyncParticle - cannot parse arguments!");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);
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
    PyObject* pyBunch;
    double x = 0.;
    double xp = 0.;
    double y = 0.;
    double yp = 0.;
    double z = 0.;
    double zp = 0.;

    int nArgs = PyTuple_Size(args);
    if(nArgs != 1 && nArgs != 7){
      error("PyBunch - addParticle - needs 6 coordinates (x,px,y,py,z,pz) or nothing!");
    }

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"O|dddddd:coordinates",&pyBunch,&x,&xp,&y,&yp,&z,&zp)){
      error("PyBunch - addParticle - cannot parse arguments!");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    int ind = cpp_bunch->addParticle(x,xp,y,yp,z,zp);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",ind);
  }


  //removes a particle to the Bunch object
  //returns the number of particles in the bunch
  //this is implementation of the deleteParticle(int index)  method
  static PyObject* Bunch_deleteParticle(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    int ind;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Oi:deleteParticle",&pyBunch,&ind)){
      error("PyBunch - deleteParticle - needs index of particle for deleting");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->deleteParticle(ind);
    int size = cpp_bunch->getSize();

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",size);
  }

  //removes a particle to the Bunch object
  //returns the index of removed macro-particle
  //this is implementation of the deleteParticleFast(int index)  method
  static PyObject* Bunch_deleteParticleFast(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    int ind;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Oi:deleteParticleFast",&pyBunch,&ind)){
      error("PyBunch - deleteParticleFast - needs index of particle for deleting");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->deleteParticleFast(ind);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",ind);
  }

  //removes all particles from the Bunch object
  //this is implementation of the deleteAllParticles()  method
  static PyObject* Bunch_deleteAllParticles(PyObject *self, PyObject *args){
    PyObject* pyBunch;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"O:deleteAllParticles",&pyBunch)){
      error("PyBunch - deleteAllParticles- needs index pyBunch object");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->deleteAllParticles();

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //compress the bunch. This method should be called after deleting one
  //  or more macro-particles
  //this is implementation of the compress()  method
  static PyObject* Bunch_compress(PyObject *self, PyObject *args){
    PyObject* pyBunch;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"O:compress",&pyBunch)){
      error("PyBunch - compress - pyBunch object needed");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->compress();

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

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
    //if nVars == 2 this is get coordinate
    //if nVars == 3 this is set coordinate
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    int index = 0;
    double val = 0.;

    if(nVars == 2 ||  nVars == 3){

      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oi:x",&pyBunch,&index)){
          error("PyBunch - x(index) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->x(index);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oid:x",&pyBunch,&index,&val)){
          error("PyBunch - x(index, value) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->x(index) = val;

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

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
    //if nVars == 2 this is get coordinate
    //if nVars == 3 this is set coordinate
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    int index = 0;
    double val = 0.;

    if(nVars == 2 ||  nVars == 3){

      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oi:y",&pyBunch,&index)){
          error("PyBunch - y(index) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->y(index);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oid:y",&pyBunch,&index,&val)){
          error("PyBunch - y(index, value) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->y(index) = val;

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

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
    //if nVars == 2 this is get coordinate
    //if nVars == 3 this is set coordinate
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    int index = 0;
    double val = 0.;

    if(nVars == 2 ||  nVars == 3){

      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oi:z",&pyBunch,&index)){
          error("PyBunch - z or phi(index) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->z(index);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oid:z",&pyBunch,&index,&val)){
          error("PyBunch - z or phi (index, value) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->z(index) = val;

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

    }
    else{
      error("PyBunch. You should call z or phi(index) or z or phi(index,value)");
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
    //if nVars == 2 this is get coordinate
    //if nVars == 3 this is set coordinate
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    int index = 0;
    double val = 0.;

    if(nVars == 2 ||  nVars == 3){

      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oi:px",&pyBunch,&index)){
          error("PyBunch - px(index) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->px(index);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oid:px",&pyBunch,&index,&val)){
          error("PyBunch - px(index, value) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->px(index) = val;

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

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
    //if nVars == 2 this is get coordinate
    //if nVars == 3 this is set coordinate
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    int index = 0;
    double val = 0.;

    if(nVars == 2 ||  nVars == 3){

      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oi:py",&pyBunch,&index)){
          error("PyBunch - py(index) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->py(index);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oid:py",&pyBunch,&index,&val)){
          error("PyBunch - py(index, value) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->py(index) = val;

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

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
    //if nVars == 2 this is get coordinate
    //if nVars == 3 this is set coordinate
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    int index = 0;
    double val = 0.;

    if(nVars == 2 ||  nVars == 3){

      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oi:pz",&pyBunch,&index)){
          error("PyBunch - pz or dE(index) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->pz(index);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oid:pz",&pyBunch,&index,&val)){
          error("PyBunch - pz or dE (index, value) - pyBunch object needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->pz(index) = val;

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

    }
    else{
      error("PyBunch. You should call pz or dE(index) or pz or dE(index,value)");
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
    //if nVars == 2 this is get flag
    //if nVars == 3 this is set flag
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    int index = 0;
    int flag = 0;

    if(nVars == 2 ||  nVars == 3){

      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oi:flag",&pyBunch,&index)){
          error("PyBunch - flag(index) - pyBunch object and index are needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        flag = cpp_bunch->flag(index);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("i",flag);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Oii:flag",&pyBunch,&index,&flag)){
          error("PyBunch - flag(index, value) - pyBunch object,index and new value are needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->flag(index) = flag;

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

    }
    else{
      error("PyBunch. You should call flag(index) or flag(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

	//Wraps long. coords in the bunch
	//ringwrap(ring_length)
  static PyObject* Bunch_ringwrap(PyObject *self, PyObject *args) {

    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    double ring_length = 0.;

    if(nVars == 2){
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple!
			//NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:py",&pyBunch,&ring_length)){
				error("PyBunch - ringwrap(ring_length) - pyBunch object needed");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
			cpp_bunch->ringwrap(ring_length);

			//clear the reference created by
			//PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);
    }
    else{
      error("PyBunch. You should call py(index) or py(index,value)");
    }

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
    //if nVars == 1 this is get mass
    //if nVars == 2 this is set mass
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    double val = 0;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:mass",&pyBunch)){
          error("PyBunch - mass() - pyBunch object is needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->getMass();

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:mass",&pyBunch,&val)){
          error("PyBunch - mass(value) - pyBunch object and new value are needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->setMass(val);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

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
    //if nVars == 1 this is get classicalRadius
    //if nVars == 2 this is set classicalRadius
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    double val = 0;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:classicalRadius",&pyBunch)){
          error("PyBunch - classicalRadius() - pyBunch object is needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->getClassicalRadius();

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:classicalRadius",&pyBunch,&val)){
          error("PyBunch - classicalRadius(value) - pyBunch object and new value are needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->setClassicalRadius(val);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

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
    //if nVars == 1 this is get charge
    //if nVars == 2 this is set charge
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    double val = 0;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:charge",&pyBunch)){
          error("PyBunch - charge() - pyBunch object is needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->getCharge();

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:charge",&pyBunch,&val)){
          error("PyBunch - charge(value) - pyBunch object and new value are needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->setCharge(val);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

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
    //if nVars == 1 this is get macroSize
    //if nVars == 2 this is set macroSize
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    double val = 0;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:macroSize",&pyBunch)){
          error("PyBunch - macroSize() - pyBunch object is needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->getMacroSize();

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:macroSize",&pyBunch,&val)){
          error("PyBunch - macroSize(value) - pyBunch object and new value are needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->setMacroSize(val);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }

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
    PyObject* pyBunch;
    const char* file_name = NULL;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Os:initBunchAttr",&pyBunch,&file_name)){
      error("PyBunch - initBunchAttr - pyBunch object and file name are needed");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->initBunchAttributes(file_name);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

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
    //if nVars == 2 this is get attribute
    //if nVars == 3 this is set attribute
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    const char* attr_name = NULL;
    double val = 0.;

    if(nVars == 2 ||  nVars == 3){

      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Os:bunchAttrDouble",&pyBunch,&attr_name)){
          error("PyBunch - bunchAttrDouble - pyBunch object and name are needed");
        }

        std::string attr_name_str(attr_name);

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->getBunchAttributeDouble(attr_name_str);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Osd:bunchAttrDouble",&pyBunch,&attr_name,&val)){
          error("PyBunch - bunchAttrDouble - pyBunch object, name, and double value are needed");
        }

        std::string attr_name_str(attr_name);

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->setBunchAttribute( attr_name_str, val);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
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
    //if nVars == 2 this is get attribute
    //if nVars == 3 this is set attribute
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    const char* attr_name = NULL;
    int val = 0;

    if(nVars == 2 ||  nVars == 3){

      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Os:bunchAttrInt",&pyBunch,&attr_name)){
          error("PyBunch - bunchAttrInt - pyBunch object and name are needed");
        }

        std::string attr_name_str(attr_name);

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->getBunchAttributeInt(attr_name_str);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("i",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Osi:bunchAttrInt",&pyBunch,&attr_name,&val)){
          error("PyBunch - bunchAttrInt - pyBunch object, name, and double value are needed");
        }

        std::string attr_name_str(attr_name);

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->setBunchAttribute( attr_name_str, val);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
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
    PyObject* pyBunch;
    std::vector<std::string> names;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"O:bunchAttrDoubleNames",&pyBunch)){
      error("PyBunch - bunchAttrDoubleNames - pyBunch object needed");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->getDoubleBunchAttributeNames(names);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

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
    PyObject* pyBunch;
    std::vector<std::string> names;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"O:bunchAttrIntNames",&pyBunch)){
      error("PyBunch - bunchAttrIntNames - pyBunch object needed");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->getIntBunchAttributeNames(names);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

   //create tuple with names
   PyObject* resTuple = PyTuple_New(names.size());

   for(int i = 0, n = names.size(); i < n; i++){
     PyObject* py_nm = PyString_FromString(names[i].c_str());
     if(PyTuple_SetItem(resTuple,i,py_nm)){
       error("PyBunch - bunchAttrIntNames - cannot create tuple with bunch attr names");
     }
   }

    return resTuple;
  }

  //Returns 0 or 1. The result is 1 if the bunch has an attribute with a particular name
  static PyObject* Bunch_hasBunchAttrDouble(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    const char* attr_name = NULL;
    int res = 0;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Os:hasBunchAttrDouble",&pyBunch,&attr_name)){
      error("PyBunch - hasBunchAttrDouble - pyBunch object and a bunch attr. name are needed");
    }

    std::string attr_name_str(attr_name);

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    res = cpp_bunch->getBunchAttributes()->hasDoubleAttribute(attr_name_str);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",res);
  }

  //Returns 0 or 1. The result is 1 if the bunch has an attribute with a particular name
  static PyObject* Bunch_hasBunchAttrInt(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    const char* attr_name = NULL;
    int res = 0;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Os:hasBunchAttrInt",&pyBunch,&attr_name)){
      error("PyBunch - hasBunchAttrInt - pyBunch object and a bunch attr. name are needed");
    }

    std::string attr_name_str(attr_name);

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    res = cpp_bunch->getBunchAttributes()->hasIntAttribute(attr_name_str);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",res);
  }

  //---------------------------------------------------------------
  //
  // related to particles' attributes
  //
  //----------------------------------------------------------------

  //Adds a particles' attributes with a particular name to the bunch
  static PyObject* Bunch_addPartAttr(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    const char* attr_name = NULL;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Os:addPartAttr",&pyBunch,&attr_name)){
      error("PyBunch - addPartAttr - pyBunch object and a particle attr. name are needed");
    }

    std::string attr_name_str(attr_name);

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->addParticleAttributes(attr_name_str);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Removes a particles' attributes with a particular name from the bunch
  static PyObject* Bunch_removePartAttr(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    const char* attr_name = NULL;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Os:removePartAttr",&pyBunch,&attr_name)){
      error("PyBunch - removePartAttr - pyBunch object and a particle attr. name are needed");
    }

    std::string attr_name_str(attr_name);

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->removeParticleAttributes(attr_name_str);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Removes all particles' attributes from the bunch
  static PyObject* Bunch_removeAllPartAttr(PyObject *self, PyObject *args){
    PyObject* pyBunch;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"O:removeAllPartAttr",&pyBunch)){
      error("PyBunch - removeAllPartAttr - pyBunch object needed");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->removeAllParticleAttributes();

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns a list (tuple) of  the particles' attributes names
  static PyObject* Bunch_getPartAttrNames(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    std::vector<std::string> names;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"O:getPartAttrNames",&pyBunch)){
      error("PyBunch -getPartAttrNames- pyBunch object needed");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->getParticleAttributesNames(names);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

   //create tuple with names
   PyObject* resTuple = PyTuple_New(names.size());

   for(int i = 0, n = names.size(); i < n; i++){
     PyObject* py_nm = PyString_FromString(names[i].c_str());
     if(PyTuple_SetItem(resTuple,i,py_nm)){
       error("PyBunch - getPartAttrNames - cannot create tuple with bunch attr names");
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
    PyObject* pyBunch;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"O:clearAllPartAttrAndMemorize",&pyBunch)){
      error("PyBunch -clearAllPartAttrAndMemorize- pyBunch object needed");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->clearAllParticleAttributesAndMemorize();

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //restores all particles' attributes names from memory
  static PyObject* Bunch_restoreAllPartAttrFromMemory(PyObject *self, PyObject *args){
    PyObject* pyBunch;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"O:restoreAllPartAttrFromMemory",&pyBunch)){
      error("PyBunch -restoreAllPartAttrFromMemory- pyBunch object needed");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->restoreAllParticleAttributesFromMemory();

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns 0 or 1. The result is 1 if the bunch has a particles' attributes with a particular name
  static PyObject* Bunch_hasPartAttr(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    const char* attr_name = NULL;
    int res = 0;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Os:hasPartAttr",&pyBunch,&attr_name)){
      error("PyBunch - hasPartAttr - pyBunch object and a particles' attr. name are needed");
    }

    std::string attr_name_str(attr_name);

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    res = cpp_bunch->hasParticleAttributes(attr_name_str);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",res);
  }

  //Returns a list (tuple) of  their bunch particles attribute names
  static PyObject* Bunch_readPartAttrNames(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    const char* file_name = NULL;
    std::vector<std::string> names;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Os:readPartAttrNames",&pyBunch,&file_name)){
      error("PyBunch - readPartAttrNames - pyBunch object and file name are needed");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->readParticleAttributesNames(file_name,names);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

   //create tuple with names
   PyObject* resTuple = PyTuple_New(names.size());

   for(int i = 0, n = names.size(); i < n; i++){
     PyObject* py_nm = PyString_FromString(names[i].c_str());
     if(PyTuple_SetItem(resTuple,i,py_nm)){
       error("PyBunch - readPartAttrNames - cannot create tuple with particles attr. names");
     }
   }

    return resTuple;
  }

  //initilizes particles' attributes from the bunch file
  static PyObject* Bunch_initPartAttr(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    const char* file_name = NULL;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Os:initPartAttr",&pyBunch,&file_name)){
      error("PyBunch - initPartAttr - pyBunch object and file name are needed");
    }

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    cpp_bunch->initParticleAttributes(file_name);

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns the number of variables in the particles' attributes with a particular name
  static PyObject* Bunch_getPartAttrSize(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    const char* attr_name = NULL;
    int size = 0;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
    //NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(	args,"Os:getPartAttrSize",&pyBunch,&attr_name)){
      error("PyBunch - getPartAttrSize - pyBunch object and a particles' attr. name are needed");
    }

    std::string attr_name_str(attr_name);

    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    size = cpp_bunch->getParticleAttributes(attr_name_str)->getAttSize();

    //clear the reference created by
    //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",size);
  }


  //Sets or returns a particles' attributes' value
  //  the action is depended on the number of arguments
  //  (attr_name,part_index, attr_index) - returns double-value
  //  (attr_name, part_index, attr_index, value) - sets the new value to the attribute
  //This is slow. In the C++ code you have to get reference to
  //particles' attributes object and operate through it
  static PyObject* Bunch_partAttrValue(PyObject *self, PyObject *args){
    //if nVars == 4 this is get attribute
    //if nVars == 5 this is set attribute
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    const char* attr_name = NULL;
    int part_index = 0;
    int attr_index = 0;
    double val = 0.;

    if(nVars == 4 ||  nVars == 5){

      if(nVars == 4){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Osii:partAttrValue",&pyBunch,&attr_name,&part_index ,&attr_index)){
          error("PyBunch - partAttrValue - pyBunch object,name of attr., part. and attr. indexes are needed");
        }

        std::string attr_name_str(attr_name);

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        val = cpp_bunch->getParticleAttributes(attr_name_str)->attValue(part_index,attr_index);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Osiid:",&pyBunch,&attr_name,&part_index ,&attr_index,&val)){
          error("PyBunch - partAttrValue - pyBunch object,name of attr., part. and attr. indexes, and value are needed");
        }

        std::string attr_name_str(attr_name);

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->getParticleAttributes(attr_name_str)->attValue(part_index,attr_index) = val;

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
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
    PyObject* pyBunch = PyTuple_GetItem(args,0);
    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    int nMacroParts = cpp_bunch->getSize();

    //clear the reference created by PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",nMacroParts);
  }

  //returns the number of macro-particles in the bunch in all CPUs
  //this is implementation of the "getSizeGlobal()" method
  static PyObject* Bunch_getSizeGlobal(PyObject *self, PyObject *args){
    PyObject* pyBunch = PyTuple_GetItem(args,0);
    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    int nMacroPartsGlobal = cpp_bunch->getSizeGlobal();

    //clear the reference created by PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",nMacroPartsGlobal);
  }

  //returns the number of macro-particles in the bunch in all CPUs
  //    that was calculated in the previous call of getSizeGlobal()
  //this is implementation of the "getSizeGlobalFromMemory()" method
  static PyObject* Bunch_getSizeGlobalFromMemory(PyObject *self, PyObject *args){
    PyObject* pyBunch = PyTuple_GetItem(args,0);
    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    int nMacroPartsGlobalFromMemory = cpp_bunch->getSizeGlobalFromMemory();

    //clear the reference created by PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",nMacroPartsGlobalFromMemory);
  }

  //returns the number of all macro-particles - alive, dead, new
  //this is implementation of the "getTotalCount()" method
  static PyObject* Bunch_getTotalCount(PyObject *self, PyObject *args){
    PyObject* pyBunch = PyTuple_GetItem(args,0);
    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    int nTotalCount = cpp_bunch->getTotalCount();

    //clear the reference created by PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",nTotalCount);
  }

  //returns the capacity of the bunch-container. It could be changed.
  //this is implementation of the "getCapacity()" method
  static PyObject* Bunch_getCapacity(PyObject *self, PyObject *args){
    PyObject* pyBunch = PyTuple_GetItem(args,0);
    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    int nCapacity = cpp_bunch->getCapacity();

    //clear the reference created by PyObject_GetAttrString( pyBunch ,"cpp_ptr")
    Py_DECREF(py_bunch_ref);

    return Py_BuildValue("i",nCapacity);
  }

  //---------------------------------------------------------------
  //
  // write into file or print Bunch
  //
  //----------------------------------------------------------------

  //Prints bunch into the std::cout stream
  static PyObject* Bunch_dumpBunch(PyObject *self, PyObject *args){
    //if nVars == 1 dumpBunchs into std::cout
    //if nVars == 2 dumpBunchs into the file
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    const char* file_name = NULL;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:dumpBunch",&pyBunch)){
          error("PyBunch - dumpBunch() - pyBunch object is needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->print(std::cout);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Os:dumpBunch",&pyBunch,&file_name)){
          error("PyBunch - dumpBunch - pyBunch object and new value are needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->print(file_name);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
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
    //if nVars == 1 reads all macro-particles
    //if nVars == 2 reads only specified number of macro-particles
    int nVars = PyTuple_Size(args);

    PyObject* pyBunch;
    const char* file_name = NULL;
    int nParts = 0;

    if(nVars == 2 ||  nVars == 3){

      if(nVars == 2){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Os:read",&pyBunch,&file_name)){
          error("PyBunch - readBunch - pyBunch object and file name are needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->readBunch(file_name);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Osi:read",&pyBunch,&file_name,&nParts)){
          error("PyBunch - readBunch - pyBunch object, file name, and number of particles are needed");
        }

        PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
        Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
        cpp_bunch->readBunch(file_name,nParts);

        //clear the reference created by
        //PyObject_GetAttrString( pyBunch ,"cpp_ptr")
        Py_DECREF(py_bunch_ref);
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
    PyObject* pyBunch;
    PyObject* pyBunch_Target;

		//NO NEW OBJECT CREATED BY PyArg_ParseTuple!
		//NO NEED OF Py_DECREF()
		if(!PyArg_ParseTuple(	args,"OO:copyEmptyBunchTo",&pyBunch,&pyBunch_Target)){
			error("PyBunch - copyEmptyBunchTo - target pyBunch object is needed");
		}

		PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
		Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

		PyObject* py_bunch_target_ref = PyObject_GetAttrString( pyBunch_Target ,"cpp_ptr");
		Bunch* cpp_target_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_target_ref);

		cpp_bunch->copyEmptyBunchTo(cpp_target_bunch);

		//clear the reference created by
		//PyObject_GetAttrString( pyBunch ,"cpp_ptr")
		Py_DECREF(py_bunch_ref);
		Py_DECREF(py_bunch_target_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Copy bunch all info including particles coordinates and attributes to another bunch
  static PyObject* Bunch_copyBunchTo(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    PyObject* pyBunch_Target;

		//NO NEW OBJECT CREATED BY PyArg_ParseTuple!
		//NO NEED OF Py_DECREF()
		if(!PyArg_ParseTuple(	args,"OO:copyBunchTo",&pyBunch,&pyBunch_Target)){
			error("PyBunch - copyBunchTo - target pyBunch object is needed");
		}

		PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
		Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

		PyObject* py_bunch_target_ref = PyObject_GetAttrString( pyBunch_Target ,"cpp_ptr");
		Bunch* cpp_target_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_target_ref);

		cpp_bunch->copyBunchTo(cpp_target_bunch);

		//clear the reference created by
		//PyObject_GetAttrString( pyBunch ,"cpp_ptr")
		Py_DECREF(py_bunch_ref);
		Py_DECREF(py_bunch_target_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Copy particles coordinates from one bunch to another
  static PyObject* Bunch_addParticlesTo(PyObject *self, PyObject *args){
    PyObject* pyBunch;
    PyObject* pyBunch_Target;

		//NO NEW OBJECT CREATED BY PyArg_ParseTuple!
		//NO NEED OF Py_DECREF()
		if(!PyArg_ParseTuple(	args,"OO:addParticlesTo",&pyBunch,&pyBunch_Target)){
			error("PyBunch - addParticlesTo - target pyBunch object is needed");
		}

		PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
		Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

		PyObject* py_bunch_target_ref = PyObject_GetAttrString( pyBunch_Target ,"cpp_ptr");
		Bunch* cpp_target_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_target_ref);

		cpp_bunch->addParticlesTo(cpp_target_bunch);

		//clear the reference created by
		//PyObject_GetAttrString( pyBunch ,"cpp_ptr")
		Py_DECREF(py_bunch_ref);
		Py_DECREF(py_bunch_target_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  //-----------------------------------------------------
  //destructor for python Bunch class
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static PyObject* Bunch_del(PyObject *self, PyObject *args){
    //std::cerr<<"The Bunch __del__ has been called!"<<std::endl;

    PyObject* pyBunch = PyTuple_GetItem(args,0);
    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

		PyObject* pySyncPart = cpp_bunch->getSyncPart()->getPyWrapper();

    delete cpp_bunch;
		Py_DECREF(pySyncPart);
    Py_DECREF(py_bunch_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }


  static PyMethodDef BunchClassMethods[] = {
    //--------------------------------------------------------
    // class Bunch wrapper                        START
    //--------------------------------------------------------
    { "__init__",                       Bunch_init                          ,METH_VARARGS,"Constructor of Bunch class"},
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
    { "bunchAttrDouble",                Bunch_bunchAttrDouble               ,METH_VARARGS,"Gets and sets a double bunch attribute"},
    { "bunchAttrInt",                   Bunch_bunchAttrInt                  ,METH_VARARGS,"Gets and sets an integer bunch attribute"},
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
    { "__del__",                        Bunch_del                           ,METH_VARARGS,"Destructor of Bunch class"},
    {NULL,NULL}
    //--------------------------------------------------------
    // class Bunch wrapper                        STOP
    //--------------------------------------------------------
  };

  static PyMethodDef BunchModuleMethods[] = { {NULL,NULL} };


#ifdef __cplusplus
extern "C" {
#endif

  void initbunch(){
    //create new module
    PyObject* module = Py_InitModule("bunch",BunchModuleMethods);
    PyObject* moduleDict = PyModule_GetDict(module);

    //create Bunch class object
    PyObject* bunchClassDict = PyDict_New();
    PyObject* bunchClassName = PyString_FromString("Bunch");
    PyObject* bunchClass = PyClass_New(NULL,bunchClassDict,bunchClassName);

    //add Bunch class to the bunch module
    PyDict_SetItemString(moduleDict,"Bunch",bunchClass);

    //clear unnecessary references
    Py_DECREF(bunchClassDict);
    Py_DECREF(bunchClassName);
    Py_DECREF(bunchClass);


    //adds methods to the Bunch class
    PyMethodDef* def;
    for( def = BunchClassMethods; def->ml_name != NULL; def++){
      PyObject* func = PyCFunction_New(def,NULL);
      PyObject* method = PyMethod_New(func,NULL,bunchClass);
      PyDict_SetItemString(bunchClassDict,def->ml_name,method);
      Py_DECREF(func);
      Py_DECREF(method);
    }

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
