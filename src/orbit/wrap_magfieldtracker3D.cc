///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "wrap_magfieldtracker3D.hh"
#include "MagneticFieldTracker3D.hh"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <string>

namespace wrap_orbit_magfieldtracker3D{

	//---------------------------------------------------------
	//Python MagneticFieldTracker3D class definition
	//---------------------------------------------------------

  //constructor for python  MagneticFieldTracker3D class
  //this is implementation of the __init__ method
  static PyObject* MagFieldTracker3D_init(PyObject *self, PyObject *args){

    int nArgs = PyTuple_Size(args);

    if(nArgs != 2 && nArgs != 1 ){
      std::cerr<<"Stop. Error in MagneticFieldTracker3D class constructor"<<std::endl;
    }

		if(nArgs == 2){
			const double step_size = 0.0;
			PyObject* pyMagFieldTracker3D = NULL;

			PyArg_ParseTuple(	args,"Od:pyMagFieldTracker3D",&pyMagFieldTracker3D,&step_size);

			MagneticFieldTracker3D* cpp_magFieldTracker3D = new MagneticFieldTracker3D(step_size);

			PyObject* py_magFieldTracker3D_ref = PyCObject_FromVoidPtr((void *) cpp_magFieldTracker3D, NULL);
			PyObject_SetAttrString(pyMagFieldTracker3D, "cpp_ptr", py_magFieldTracker3D_ref);

			Py_DECREF(py_magFieldTracker3D_ref);
		}

		if(nArgs == 1){
			PyObject* pyMagFieldTracker3D = PyTuple_GetItem(args,0);
			MagneticFieldTracker3D* cpp_magFieldTracker3D = new MagneticFieldTracker3D();

			PyObject* py_magFieldTracker3D_ref = PyCObject_FromVoidPtr((void *) cpp_magFieldTracker3D, NULL);
			PyObject_SetAttrString(pyMagFieldTracker3D, "cpp_ptr", py_magFieldTracker3D_ref);

			Py_DECREF(py_magFieldTracker3D_ref);
		}

    Py_INCREF(Py_None);
    return Py_None;
  }


  //-----------------------------------------------------
  //This method returns the stepSize of the tarcker
  //-----------------------------------------------------
static PyObject* 	MagFieldTracker3D_getStepSize(PyObject *self, PyObject *args)
{        
   int nVars = PyTuple_Size(args);
   if(nVars != 1)
   {
      std::cerr << "Stop.  Error in MagneticFieldTracker3D class method getStepSize" << endl;
   }

   double stepSize = 0.0;

   PyObject* pyMagFieldTracker3D = PyTuple_GetItem(args,0);
   PyObject* py_magfieldtracker3D_ref = PyObject_GetAttrString(pyMagFieldTracker3D,"cpp_ptr");
   MagneticFieldTracker3D* cpp_magfieldtracker3D = (MagneticFieldTracker3D*) PyCObject_AsVoidPtr(py_magfieldtracker3D_ref);

   stepSize = cpp_magfieldtracker3D->getStepSize();
   
   Py_DECREF(py_magfieldtracker3D_ref);   
   return Py_BuildValue("d",stepSize);
  }

  
  //-----------------------------------------------------
  //This method returns the stepSize of the tarcker
  //-----------------------------------------------------
static PyObject* 	MagFieldTracker3D_setStepSize(PyObject *self, PyObject *args)
{
   int nVars = PyTuple_Size(args);
   if(nVars != 2)
   {
      std::cerr << "Stop.  Error in MagneticFieldTracker3D class method setStepSize" << endl;
   }
  
   double stepSize = 0.0;
   PyObject* pyMagFieldTracker3D = NULL;
   PyArg_ParseTuple(args,"Od:pyMagFieldTracker3D",&pyMagFieldTracker3D,&stepSize);
   
   PyObject* py_magfieldtracker3D_ref = PyObject_GetAttrString(pyMagFieldTracker3D,"cpp_ptr");
   MagneticFieldTracker3D* cpp_magfieldtracker3D = (MagneticFieldTracker3D*)PyCObject_AsVoidPtr(py_magfieldtracker3D_ref);
   cpp_magfieldtracker3D->setStepSize(stepSize);
   
   Py_DECREF(py_magfieldtracker3D_ref);

   Py_INCREF(Py_None);
   return Py_None;
}



  //-----------------------------------------------------
  //This method tracks the bunch through a given magnetic field
  //-----------------------------------------------------  
    static PyObject* MagFieldTracker3D_track(PyObject *self, PyObject *args){
    int nVars = PyTuple_Size(args);
    if(nVars != 11){
        std::cerr << "Stop.  Error in MagneticFieldTracker3D class method track" << endl;
     }

    double zmin = 0.0;
    double zmax = 0.0;
    double x_0 = 0.0;
    double y_0 = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double gamma = 0.0;
    int method_flag = 0;
    
    PyObject* pyBaseFieldSource = NULL;
    PyObject* pyBunch = NULL;
    PyObject* pyMagFieldTracker3D = NULL;
	 
    PyArg_ParseTuple(args,"OdddddddOOi:pyMagFieldTracker3D",&pyMagFieldTracker3D,&zmin,
                   &zmax,&x_0,&y_0,&alpha,&beta,&gamma,&pyBunch,&pyBaseFieldSource, &method_flag);  
               
    PyObject* py_magfieldtracker3D_ref = PyObject_GetAttrString(pyMagFieldTracker3D,"cpp_ptr");
    MagneticFieldTracker3D* cpp_magFieldTracker3D = (MagneticFieldTracker3D*) PyCObject_AsVoidPtr(py_magfieldtracker3D_ref);
    
    PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
    Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);
    
    PyObject* py_basefieldsource_ref = PyObject_GetAttrString(pyBaseFieldSource,"cpp_ptr");
    BaseFieldSource* cpp_baseFieldSource = (BaseFieldSource*) PyCObject_AsVoidPtr(py_basefieldsource_ref);
    
    cpp_magFieldTracker3D->track(zmin,zmax,x_0,y_0,alpha,beta,gamma,cpp_bunch,cpp_baseFieldSource, method_flag);

    Py_DECREF(py_magfieldtracker3D_ref);
    Py_DECREF(py_bunch_ref);
    Py_DECREF(py_basefieldsource_ref);
    
    //create tuple with names
    PyObject* resTuple = PyTuple_New(3);

    PyTuple_SetItem(resTuple,0,Py_BuildValue("d",alpha));
    PyTuple_SetItem(resTuple,1,Py_BuildValue("d",beta));
    PyTuple_SetItem(resTuple,2,Py_BuildValue("d",gamma));

    return resTuple;
  }
  
  

  //-----------------------------------------------------
  //destructor for python MagneticFieldTracker3D class
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static PyObject* MagFieldTracker3D_del(PyObject *self, PyObject *args){

    PyObject* pyMagFieldTracker3D = PyTuple_GetItem(args,0);
    PyObject* py_magFieldTracker3D_ref = PyObject_GetAttrString( pyMagFieldTracker3D ,"cpp_ptr");
		MagneticFieldTracker3D* magFieldTracker3D = (MagneticFieldTracker3D*) PyCObject_AsVoidPtr(py_magFieldTracker3D_ref);

		delete magFieldTracker3D;
    Py_DECREF(py_magFieldTracker3D_ref);

    Py_INCREF(Py_None);
    return Py_None;
  }

  static PyMethodDef magFieldTracker3DClassMethods[] = {
    //--------------------------------------------------------
    // class MagneticFieldTracker3D wrapper                        START
    //--------------------------------------------------------
    { "__init__",              MagFieldTracker3D_init           ,METH_VARARGS,"Constructor. Creates MagneticFieldTracker3D class"},
    { "track",                 MagFieldTracker3D_track          ,METH_VARARGS,"Tracks bunch through given magnetic field returning tuple of new alpha,beta,and gamma"},
    { "getStepSize",           MagFieldTracker3D_getStepSize    ,METH_VARARGS,"Returns the step size used by the tracker"},
    { "setStepSize",           MagFieldTracker3D_setStepSize    ,METH_VARARGS,"Sets the step size to be used during track"},
    { "__del__",               MagFieldTracker3D_del            ,METH_VARARGS,"Destructor of MagneticFieldTracker3D class"},
    {NULL,NULL}
    //--------------------------------------------------------
    // class MagneticFieldTracker3D wrapper       STOP
    //--------------------------------------------------------
  };

#ifdef __cplusplus
extern "C" {
#endif

  void initmagfieldtracker3D(PyObject* module){


    //create new module
    PyObject* moduleDict = PyModule_GetDict(module);

    //create magFieldTracker3D class object
    PyObject* magFieldTracker3DClassDict = PyDict_New();
    PyObject* magFieldTracker3DClassName = PyString_FromString("MagneticFieldTracker3D");
    PyObject* magFieldTracker3DClass = PyClass_New(NULL,magFieldTracker3DClassDict,magFieldTracker3DClassName);

    //add MagneticFieldTracker3D class to the fieldtracker? module
    PyDict_SetItemString(moduleDict,"MagneticFieldTracker3D",magFieldTracker3DClass);

    //clear unnecessary references
    Py_DECREF(magFieldTracker3DClassDict);
    Py_DECREF(magFieldTracker3DClassName);
    Py_DECREF(magFieldTracker3DClass);


    //adds methods to the MagneticFieldTracker class
    PyMethodDef* def;
    for( def = magFieldTracker3DClassMethods; def->ml_name != NULL; def++){
      PyObject* func = PyCFunction_New(def,NULL);
      PyObject* method = PyMethod_New(func,NULL,magFieldTracker3DClass);
      PyDict_SetItemString(magFieldTracker3DClassDict,def->ml_name,method);
      Py_DECREF(func);
      Py_DECREF(method);
    }
  }


#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_magfieldtracker3D
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
