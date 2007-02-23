///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "wrap_multexpansion3D.hh"
#include "MultipoleExpansion3D.hh"

#include <iostream>
#include <cstdlib>
#include <string>

namespace wrap_orbit_multexpansion3D{

	//---------------------------------------------------------
	//Python MultipoleExpansion3D class definition
	//---------------------------------------------------------

  //constructor for python  MultipoleExpansion3D class
  //this is implementation of the __init__ method
  static PyObject* MultExp3D_init(PyObject *self, PyObject *args){

    int nArgs = PyTuple_Size(args);

    if(nArgs != 2 && nArgs != 1 ){
      std::cerr<<"Stop. Error in MultipoleExpansion3D class constructor"<<std::endl;
    }

		if(nArgs == 2){
			const char* file_name = NULL;
			PyObject* pyMultExp3D = NULL;

			PyArg_ParseTuple(	args,"Os:pyMultExp3D",&pyMultExp3D,&file_name);

			std::string file_name_str(file_name);
			MultipoleExpansion3D* cpp_multExp3D = new MultipoleExpansion3D(file_name_str);

			PyObject* py_multExp3D_ref = PyCObject_FromVoidPtr((void *) cpp_multExp3D, NULL);
			PyObject_SetAttrString(pyMultExp3D, "cpp_ptr", py_multExp3D_ref);

			Py_DECREF(py_multExp3D_ref);
		}

		if(nArgs == 1){
			PyObject* pyMultExp3D = PyTuple_GetItem(args,0);
			MultipoleExpansion3D* cpp_multExp3D = new MultipoleExpansion3D();

			PyObject* py_multExp3D_ref = PyCObject_FromVoidPtr((void *) cpp_multExp3D, NULL);
			PyObject_SetAttrString(pyMultExp3D, "cpp_ptr", py_multExp3D_ref);

			Py_DECREF(py_multExp3D_ref);
		}

    Py_INCREF(Py_None);
    return Py_None;
  }


	//-----------------------------------------------------
  //This method returns the tuple with 3 components of the magnetic field
  //-----------------------------------------------------
  static PyObject* 	MultExp3D_getMagnField(PyObject *self, PyObject *args){

		int nArgs = PyTuple_Size(args);

    if(nArgs != 5 ){
      std::cerr<<"Stop. Error in MultipoleExpansion3D class method getMagneticField."<<std::endl;
    }

		double t = 0.0;
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;

		double bx = 0.0;
		double by = 0.0;
		double bz = 0.0;

		PyObject* pyMultExp3D = NULL;

		PyArg_ParseTuple(	args,"Odddd:pyMultExp3D",&pyMultExp3D,&t,&x,&y,&z);

    PyObject* py_multExp3D_ref = PyObject_GetAttrString( pyMultExp3D ,"cpp_ptr");
		MultipoleExpansion3D* multExp3D = (MultipoleExpansion3D*) PyCObject_AsVoidPtr(py_multExp3D_ref);
		multExp3D->getMagneticField(t,x,y,z,bx,by,bz);

		Py_DECREF(py_multExp3D_ref);

		//create tuple with names
		PyObject* resTuple = PyTuple_New(3);

		PyTuple_SetItem(resTuple,0,Py_BuildValue("d",bx));
		PyTuple_SetItem(resTuple,1,Py_BuildValue("d",by));
		PyTuple_SetItem(resTuple,2,Py_BuildValue("d",bz));

		return resTuple;
	}

  //-----------------------------------------------------
  //destructor for python MultipoleExpansion3D class
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static PyObject* MultExp3D_del(PyObject *self, PyObject *args){
    //std::cerr<<"The SyncParticle __del__ has been called!"<<std::endl;

    PyObject* pyMultExp3D = PyTuple_GetItem(args,0);
    PyObject* py_multExp3D_ref = PyObject_GetAttrString( pyMultExp3D ,"cpp_ptr");
		MultipoleExpansion3D* multExp3D = (MultipoleExpansion3D*) PyCObject_AsVoidPtr(py_multExp3D_ref);

		delete multExp3D;
    Py_DECREF(py_multExp3D_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  static PyMethodDef multExp3DClassMethods[] = {
    //--------------------------------------------------------
    // class MultipoleExpansion3D wrapper                        START
    //--------------------------------------------------------
    { "__init__",              MultExp3D_init           ,METH_VARARGS,"Constructor of MultipoleExpansion3D class"},
    { "getMagneticField",      MultExp3D_getMagnField   ,METH_VARARGS,"Returns the tuple of 3 components of the magnetic field"},
    { "__del__",               MultExp3D_del            ,METH_VARARGS,"Destructor of MultipoleExpansion3D class"},
    {NULL,NULL}
    //--------------------------------------------------------
    // class MultipoleExpansion3D wrapper       STOP
    //--------------------------------------------------------
  };

#ifdef __cplusplus
extern "C" {
#endif

  void initmultexpansion3D(PyObject* module){


    //create new module
    PyObject* moduleDict = PyModule_GetDict(module);

    //create multExp3D class object
    PyObject* multExp3DClassDict = PyDict_New();
    PyObject* multExp3DClassName = PyString_FromString("MultipoleExpansion3D");
    PyObject* multExp3DClass = PyClass_New(NULL,multExp3DClassDict,multExp3DClassName);

    //add Bunch class to the bunch module
    PyDict_SetItemString(moduleDict,"MultipoleExpansion3D",multExp3DClass);

    //clear unnecessary references
    Py_DECREF(multExp3DClassDict);
    Py_DECREF(multExp3DClassName);
    Py_DECREF(multExp3DClass);


    //adds methods to the Bunch class
    PyMethodDef* def;
    for( def = multExp3DClassMethods; def->ml_name != NULL; def++){
      PyObject* func = PyCFunction_New(def,NULL);
      PyObject* method = PyMethod_New(func,NULL,multExp3DClass);
      PyDict_SetItemString(multExp3DClassDict,def->ml_name,method);
      Py_DECREF(func);
      Py_DECREF(method);
    }
  }


#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_multexpansion3D
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
