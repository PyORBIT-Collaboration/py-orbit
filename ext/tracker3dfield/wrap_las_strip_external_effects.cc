#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_las_strip_external_effects.hh"

#include <iostream>
#include <string>

#include "LasStripExternalEffects.hh"

using namespace OrbitUtils;
using namespace Tracker3DField;

namespace wrap_tracker3dfield_lasstrip_external_effects{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppExternalEffects class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppExternalEffects instance
	//It never will be called directly
	static PyObject* LasStripExternalEffects_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  PyExternalEffects class
  //this is implementation of the __init__ method
  static int LasStripExternalEffects_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  
	  double half_Angle;
	  double Power;
	  double Lambda;
	  int states;
	  char* addressEG;
	  
		 if(!PyArg_ParseTuple(	args,"dddsi:",&half_Angle,&Power ,&Lambda,&addressEG,&states)){
			 		          error("LaserExternalEfects(half_Angle,LaserPower,Lambda) - params. are needed");
			 			 		        }  
			else self->cpp_obj =  new  LasStripExternalEffects(half_Angle,Power,Lambda,addressEG,states);
    return 0;
  }
  
 /*
	void	LasStripExternalEffects::setLaserHalfAngle(double a);
	double	LasStripExternalEffects::getLaserHalfAngle();
	void	LasStripExternalEffects::setLaserPower(double a);
	double	LasStripExternalEffects::getLaserPower();
	void	LasStripExternalEffects::setLaser_lambda(double a);
	double	LasStripExternalEffects::getLaser_lambda();
  */
  
static PyObject* LasStripExternalEffects_setLaserHalfAngle(PyObject *self, PyObject *args){
 		LasStripExternalEffects* cpp_LasStripExternalEffects = (LasStripExternalEffects*)((pyORBIT_Object*) self)->cpp_obj;
 		  double half_Angle;	       
 		        if(!PyArg_ParseTuple(	args,"d:",&half_Angle)){
 		          error("LaserExternalEfects - setLaserHalfAngle(half_Angle) - param. is needed");
 		        }
 		        else cpp_LasStripExternalEffects->setLaserHalfAngle(half_Angle);
 		      
 	   
 		    Py_INCREF(Py_None);
 		    return Py_None;	  
  }		

static PyObject* LasStripExternalEffects_setLaserPower(PyObject *self, PyObject *args){
		LasStripExternalEffects* cpp_LasStripExternalEffects = (LasStripExternalEffects*)((pyORBIT_Object*) self)->cpp_obj;
		  double LaserPower;	       
		        if(!PyArg_ParseTuple(	args,"d:",&LaserPower)){
		          error("LaserExternalEfects - setLaserPower(LaserPower) - param. is needed");
		        }
		        else cpp_LasStripExternalEffects->setLaserPower(LaserPower);
		      
	   
		    Py_INCREF(Py_None);
		    return Py_None;	  
}

static PyObject* LasStripExternalEffects_setLaser_lambda(PyObject *self, PyObject *args){
		LasStripExternalEffects* cpp_LasStripExternalEffects = (LasStripExternalEffects*)((pyORBIT_Object*) self)->cpp_obj;
		  double Laser_lambda;	       
		        if(!PyArg_ParseTuple(	args,"d:",&Laser_lambda)){
		          error("LaserExternalEfects - setLaser_lambda(Laser_lambda) - param. is needed");
		        }
		        else cpp_LasStripExternalEffects->setLaserPower(Laser_lambda);
		      
	   
		    Py_INCREF(Py_None);
		    return Py_None;	  
}




static PyObject* LasStripExternalEffects_getLaserHalfAngle(PyObject *self, PyObject *args){
		LasStripExternalEffects* cpp_LasStripExternalEffects = (LasStripExternalEffects*)((pyORBIT_Object*) self)->cpp_obj;
				       
		  return Py_BuildValue("d",cpp_LasStripExternalEffects->getLaserHalfAngle()); 
}		
    

static PyObject* LasStripExternalEffects_getLaserPower(PyObject *self, PyObject *args){
		LasStripExternalEffects* cpp_LasStripExternalEffects = (LasStripExternalEffects*)((pyORBIT_Object*) self)->cpp_obj;
				       
		  return Py_BuildValue("d",cpp_LasStripExternalEffects->getLaserPower()); 
}		    

static PyObject* LasStripExternalEffects_getLaser_lambda(PyObject *self, PyObject *args){
		LasStripExternalEffects* cpp_LasStripExternalEffects = (LasStripExternalEffects*)((pyORBIT_Object*) self)->cpp_obj;
				       
		  return Py_BuildValue("d",cpp_LasStripExternalEffects->getLaser_lambda()); 
}		    

		
	  
 		  


  
  
	// name([name]) - sets or returns the name of the External Effeects class 
  static PyObject* LasStripExternalEffects_name(PyObject *self, PyObject *args){
		LasStripExternalEffects* cpp_LasStripExternalEffects = (LasStripExternalEffects*) ((pyORBIT_Object*) self)->cpp_obj;
    const char* name = NULL;
    if(!PyArg_ParseTuple(	args,"|s:name",&name)){
      error("LasStripExternalEffects - call should be - name([name]).");
    }
		if(name != NULL){
      std::string name_str(name);
      cpp_LasStripExternalEffects->setName(name_str);
		}
		return Py_BuildValue("s",cpp_LasStripExternalEffects->getName().c_str());
  }	
	
  //-----------------------------------------------------
  //destructor for python PyExternalEffects class (__del__ method).
  //-----------------------------------------------------
  static void LasStripExternalEffects_del(pyORBIT_Object* self){
		//std::cerr<<"The LasStripExternalEffects __del__ has been called!"<<std::endl;
		delete ((LasStripExternalEffects*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyExternalEffects wrapper class
	// they will be vailable from python level
  static PyMethodDef LasStripExternalEffectsClassMethods[] = {
		{ "name",         LasStripExternalEffects_name,         METH_VARARGS,"Sets or returns the name of effects."},
		{ "setLaserHalfAngle",         LasStripExternalEffects_setLaserHalfAngle,         METH_VARARGS,"Sets or returns the name of effects."},
		{ "setLaserPower",         LasStripExternalEffects_setLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
		{ "setLaser_lambda",         LasStripExternalEffects_setLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
		{ "getLaserHalfAngle",         LasStripExternalEffects_getLaserHalfAngle,         METH_VARARGS,"Sets or returns the name of effects."},
		{ "getLaserPower",         LasStripExternalEffects_getLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
		{ "getLaser_lambda",         LasStripExternalEffects_getLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
    {NULL}
  };

	// defenition of the memebers of the python PyExternalEffects wrapper class
	// they will be vailable from python level
	static PyMemberDef LasStripExternalEffectsClassMembers [] = {
		{NULL}
	};

	//new python PyExternalEffects wrapper type definition
	static PyTypeObject pyORBIT_LasStripExternalEffects_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"CppExternalEffects", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) LasStripExternalEffects_del , /*tp_dealloc*/
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
		"The LasStripExternalEffects python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		LasStripExternalEffectsClassMethods, /* tp_methods */
		LasStripExternalEffectsClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) LasStripExternalEffects_init, /* tp_init */
		0, /* tp_alloc */
		LasStripExternalEffects_new, /* tp_new */
	};	


		
	//--------------------------------------------------
	//Initialization function of the pyPyExternalEffects class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initLasStripExternalEffects(PyObject* module){
		if (PyType_Ready(&pyORBIT_LasStripExternalEffects_Type) < 0) return;
		Py_INCREF(&pyORBIT_LasStripExternalEffects_Type);
		PyModule_AddObject(module, "LasStripExternalEffects", (PyObject *)&pyORBIT_LasStripExternalEffects_Type);
				
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_tracker3dfield
}
