///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "wrap_syncpart.hh"

#include "Bunch.hh"
#include "SyncPart.hh"

namespace wrap_orbit_syncpart{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

	//---------------------------------------------------------
	//Python SyncParticle class definition
	//---------------------------------------------------------

  //constructor for python  SyncParticle class
  //this is implementation of the __init__ method
  static PyObject* SyncPart_init(PyObject *self, PyObject *args){
    //std::cerr<<"The SyncParticle __init__ has been called!"<<std::endl;

    int nArgs = PyTuple_Size(args);
    if(nArgs != 2){
      error("SyncParticle constructor needs pyBunch instance as parameter.");
    }

    PyObject* pySyncPart = PyTuple_GetItem(args,0);
		PyObject* pyBunch = PyTuple_GetItem(args,1);
		PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

    SyncPart* cpp_syncPart = cpp_bunch->getSyncPart();
		if(cpp_syncPart->getPyWrapper() != NULL){
			error("You should not create SyncParticle class instance directly!");
		}

		cpp_syncPart->setPyWrapper(pySyncPart);
		Py_INCREF(pySyncPart);

    PyObject* py_syncPart_ref = PyCObject_FromVoidPtr((void *) cpp_syncPart, NULL);
    if(PyObject_SetAttrString(pySyncPart	, "cpp_ptr", py_syncPart_ref) < 0){
      error("SyncParticle  wrapper. Can not set c-SyncParticle reference.");
    }

    //clear the reference created by PyObject_GetAttrString(...)
		// and PyCObject_FromVoidPtr(...)
		// According to Python documentation they create "New Reference"
    Py_DECREF(py_syncPart_ref);
		Py_DECREF(py_bunch_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

 //Sets or returns the kinEnergy for the SyncPart object
  //  the action is depended on the number of arguments
  //  kinEnergy() - returns kinEnergy MeV
  //  kinEnergy(value) - sets the new value for pz and px=0,py=0
  static PyObject* SyncPart_kinEnergy(PyObject *self, PyObject *args){
    //if nVars == 1 this is get kinEnergy
    //if nVars == 2 this is set kinEnergy
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:kinEnergy",&pySyncPart)){
          error("PySyncPart - kinEnergy() - pySyncPart object is needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->getEnergy();

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:kinEnergy",&pySyncPart,&val)){
          error("PySyncPart - kinEnergy(value) - pySyncPart object and new value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
				val = cpp_SyncPart->energyToMomentum(val);
        cpp_SyncPart->setPXYZ(0.,0.,val);
				val = cpp_SyncPart->momentumToEnergy(val);

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);
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
    PyObject* pySyncPart;
    double val = 0.;

		if(!PyArg_ParseTuple(	args,"O:beta",&pySyncPart)){
			error("PySyncPart - beta() - pySyncPart object is needed");
		}

		PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
		SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
		val = cpp_SyncPart->getBeta();

		Py_DECREF(py_SyncPart_ref);

		return Py_BuildValue("d",val);
  }

  //  gamma() - returns gamma
  static PyObject* SyncPart_gamma(PyObject *self, PyObject *args){
    PyObject* pySyncPart;
    double val = 0.;

		if(!PyArg_ParseTuple(	args,"O:gamma",&pySyncPart)){
			error("PySyncPart - gamma() - pySyncPart object is needed");
		}

		PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
		SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
		val = cpp_SyncPart->getGamma();

		Py_DECREF(py_SyncPart_ref);

		return Py_BuildValue("d",val);
  }

  //  mass() - returns mass
  static PyObject* SyncPart_mass(PyObject *self, PyObject *args){
    PyObject* pySyncPart;
    double val = 0.;

		if(!PyArg_ParseTuple(	args,"O:mass",&pySyncPart)){
			error("PySyncPart - mass() - pySyncPart object is needed");
		}

		PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
		SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
		val = cpp_SyncPart->getMass();

		Py_DECREF(py_SyncPart_ref);

		return Py_BuildValue("d",val);
  }

 //Sets or returns the momentum for the SyncPart object
  //  the action is depended on the number of arguments
  //  momentum() - returns momentum
  //  momentum(value) - sets the new value for pz and px=0,py=0
  static PyObject* SyncPart_momentum(PyObject *self, PyObject *args){
    //if nVars == 1 this is get momentum
    //if nVars == 2 this is set momentum
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:momentum",&pySyncPart)){
          error("PySyncPart - momentum() - pySyncPart object is needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->getMomentum();

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:momentum",&pySyncPart,&val)){
          error("PySyncPart - momentum(value) - pySyncPart object and new value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        cpp_SyncPart->setPXYZ(0.,0.,val);

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);
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
  //  time() - returns time
  //  time(value) - sets the new value
  static PyObject* SyncPart_time(PyObject *self, PyObject *args){
    //if nVars == 1 this is get time
    //if nVars == 2 this is set time
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:time",&pySyncPart)){
          error("PySyncPart - time() - pySyncPart object is needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->getTime();

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:time",&pySyncPart,&val)){
          error("PySyncPart - time(value) - pySyncPart object and new value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        cpp_SyncPart->setTime(val);

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);
      }

    }
    else{
      error("PySyncPart. You should call time() or time(value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

 //Sets or returns the rfFrequency in 1/second for the SyncPart object
  //  the action is depended on the number of arguments
  //  rfFrequency() - returns rfFrequency
  //  rfFrequency(value) - sets the new value
  static PyObject* SyncPart_rfFrequency(PyObject *self, PyObject *args){
    //if nVars == 1 this is get rfFrequency
    //if nVars == 2 this is set rfFrequency
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:rfFrequency",&pySyncPart)){
          error("PySyncPart - rfFrequency() - pySyncPart object is needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->getFrequency();

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:rfFrequency",&pySyncPart,&val)){
          error("PySyncPart - rfFrequency(value) - pySyncPart object and new value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        cpp_SyncPart->setFrequency(val);

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);
      }

    }
    else{
      error("PySyncPart. You should call rfFrequency() or rfFrequency(value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }


 //Sets or returns the x for the SyncPart object
  //  the action is depended on the number of arguments
  //  x() - returns x
  //  x(value) - sets the new value
  static PyObject* SyncPart_x(PyObject *self, PyObject *args){
    //if nVars == 1 this is get x
    //if nVars == 2 this is set x
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:x",&pySyncPart)){
          error("PySyncPart - x() - pySyncPart object is needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->getX();

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:x",&pySyncPart,&val)){
          error("PySyncPart - x(value) - pySyncPart object and new value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        cpp_SyncPart->setX(val);

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);
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
    //if nVars == 1 this is get y
    //if nVars == 2 this is set y
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:y",&pySyncPart)){
          error("PySyncPart - y() - pySyncPart object is needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->getY();

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:y",&pySyncPart,&val)){
          error("PySyncPart - y(value) - pySyncPart object and new value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        cpp_SyncPart->setY(val);

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);
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
    //if nVars == 1 this is get z
    //if nVars == 2 this is set z
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:z",&pySyncPart)){
          error("PySyncPart - z() - pySyncPart object is needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->getZ();

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:z",&pySyncPart,&val)){
          error("PySyncPart - z(value) - pySyncPart object and new value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        cpp_SyncPart->setZ(val);

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);
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
    //if nVars == 1 this is get px
    //if nVars == 2 this is set px
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:px",&pySyncPart)){
          error("PySyncPart - px() - pySyncPart object is needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->getPX();

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:px",&pySyncPart,&val)){
          error("PySyncPart - px(value) - pySyncPart object and new value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        cpp_SyncPart->setPX(val);

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);
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
    //if nVars == 1 this is get py
    //if nVars == 2 this is set py
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:py",&pySyncPart)){
          error("PySyncPart - py() - pySyncPart object is needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->getPY();

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:py",&pySyncPart,&val)){
          error("PySyncPart - py(value) - pySyncPart object and new value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        cpp_SyncPart->setPY(val);

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);
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
    //if nVars == 1 this is get pz
    //if nVars == 2 this is set pz
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"O:pz",&pySyncPart)){
          error("PySyncPart - pz() - pySyncPart object is needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->getPZ();

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);

        return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"Od:pz",&pySyncPart,&val)){
          error("PySyncPart - pz(value) - pySyncPart object and new value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        cpp_SyncPart->setPZ(val);

        //clear the reference created by
        //PyObject_GetAttrString( pySyncPart ,"cpp_ptr")
        Py_DECREF(py_SyncPart_ref);
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
    //should be nVars == 2
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 2){
        if(!PyArg_ParseTuple(	args,"Od:momentumToEnergy",&pySyncPart,&val)){
          error("PySyncPart - momentumToEnergy(value) - pySyncPart object and input value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->momentumToEnergy(val);

        Py_DECREF(py_SyncPart_ref);
    }
    else{
      error("PySyncPart. You should call momentumToEnergy(p)");
    }

		return Py_BuildValue("d",val);
  }

  //Transforms kinetic energy into momentum
  static PyObject* SyncPart_eToP(PyObject *self, PyObject *args){
    //should be nVars == 2
    int nVars = PyTuple_Size(args);

    PyObject* pySyncPart;
    double val = 0.;

    if(nVars == 2){
        if(!PyArg_ParseTuple(	args,"Od:energyToMomentum",&pySyncPart,&val)){
          error("PySyncPart - energyToMomentum(value) - pySyncPart object and input value are needed");
        }

        PyObject* py_SyncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
        SyncPart* cpp_SyncPart = (SyncPart*) PyCObject_AsVoidPtr(py_SyncPart_ref);
        val = cpp_SyncPart->energyToMomentum(val);

        Py_DECREF(py_SyncPart_ref);
    }
    else{
      error("PySyncPart. You should call energyToMomentum(p)");
    }

		return Py_BuildValue("d",val);
  }

  //-----------------------------------------------------
  //destructor for python SyncParticle class
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static PyObject* SyncPart_del(PyObject *self, PyObject *args){
    //std::cerr<<"The SyncParticle __del__ has been called!"<<std::endl;

    PyObject* pySyncPart = PyTuple_GetItem(args,0);
    PyObject* py_syncPart_ref = PyObject_GetAttrString( pySyncPart ,"cpp_ptr");
    Py_DECREF(py_syncPart_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  static PyMethodDef SyncPartClassMethods[] = {
    //--------------------------------------------------------
    // class SyncParticle wrapper                        START
    //--------------------------------------------------------
    { "__init__",         SyncPart_init         ,METH_VARARGS,"Constructor of SyncParticle class"},
    { "mass",             SyncPart_mass         ,METH_VARARGS,"Returns mass in MeV"},
    { "momentum",         SyncPart_momentum     ,METH_VARARGS,"Returns or sets momentum in MeV/c. If setting px=0,py=0,pz=val."},
    { "beta",             SyncPart_beta         ,METH_VARARGS,"Returns beta=v/c"},
		{ "gamma",            SyncPart_gamma        ,METH_VARARGS,"Returns gamma=1/sqrt(1-(v/c)**2)"},
		{ "kinEnergy",        SyncPart_kinEnergy    ,METH_VARARGS,"Returns or sets kinetic energy of the synchronous particle in MeV"},
		{ "rfFrequency",      SyncPart_rfFrequency  ,METH_VARARGS,"Returns or sets the rf frequency in Hz"},
		{ "time",		          SyncPart_time         ,METH_VARARGS,"Sets or returns time in sec"},
		{ "x",		            SyncPart_x            ,METH_VARARGS,"Sets or returns the x-coordinate"},
		{ "y",		            SyncPart_y            ,METH_VARARGS,"Sets or returns the y-coordinate"},
		{ "z",		            SyncPart_z            ,METH_VARARGS,"Sets or returns the z-coordinate"},
		{ "px",		            SyncPart_px           ,METH_VARARGS,"Sets or returns the x-momentum"},
		{ "py",		            SyncPart_py           ,METH_VARARGS,"Sets or returns the y-momentum"},
		{ "pz",		            SyncPart_pz           ,METH_VARARGS,"Sets or returns the z-momentum"},
		{ "energyToMomentum", SyncPart_eToP         ,METH_VARARGS,"Transforms the energy to momentum"},
		{ "momentumToEnergy", SyncPart_pToE         ,METH_VARARGS,"Transforms the momentum to energy"},
    { "__del__",          SyncPart_del          ,METH_VARARGS,"Destructor of SyncParticle class"},
    {NULL,NULL}
    //--------------------------------------------------------
    // class SyncParticle wrapper                        STOP
    //--------------------------------------------------------
  };

#ifdef __cplusplus
extern "C" {
#endif

  void initsyncpart(PyObject* module){


    //create new module
    PyObject* moduleDict = PyModule_GetDict(module);

    //create SyncPart class object
    PyObject* syncPartClassDict = PyDict_New();
    PyObject* syncPartClassName = PyString_FromString("SyncParticle");
    PyObject* syncPartClass = PyClass_New(NULL,syncPartClassDict,syncPartClassName);

    //add Bunch class to the bunch module
    PyDict_SetItemString(moduleDict,"SyncParticle",syncPartClass);

    //clear unnecessary references
    Py_DECREF(syncPartClassDict);
    Py_DECREF(syncPartClassName);
    Py_DECREF(syncPartClass);


    //adds methods to the Bunch class
    PyMethodDef* def;
    for( def = SyncPartClassMethods; def->ml_name != NULL; def++){
      PyObject* func = PyCFunction_New(def,NULL);
      PyObject* method = PyMethod_New(func,NULL,syncPartClass);
      PyDict_SetItemString(syncPartClassDict,def->ml_name,method);
      Py_DECREF(func);
      Py_DECREF(method);
    }
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
