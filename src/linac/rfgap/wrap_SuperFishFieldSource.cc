#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_SuperFishFieldSource.hh"
#include "wrap_linacmodule.hh"

#include <iostream>

#include "wrap_utils.hh"
#include "wrap_spacecharge.hh"
#include "SuperFishFieldSource.hh"

using namespace OrbitUtils;

namespace wrap_linac{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python SuperFishFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping SuperFishFieldSource instance
	//It never will be called directly
	static PyObject* SuperFishFieldSource_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The SuperFishFieldSource new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  SuperFishFieldSource class
  //this is implementation of the __init__ method
  static int SuperFishFieldSource_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new SuperFishFieldSource();	
		((SuperFishFieldSource*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		return 0;
  }

	//SuperFishFieldSource_getEMField - returns components of the electric and magnetic filds.
  static PyObject* SuperFishFieldSource_getElectricMagneticField(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;
	  double t, x, y, z, E_x, E_y,  E_z, H_x, H_y,  H_z;		
		if(!PyArg_ParseTuple(args,"dddd:getElectricMagneticField",&x,&y,&z,&t)){
			ORBIT_MPI_Finalize("PySuperFishFieldSource - getElectricMagneticField(x,y,z,t) - parameters are needed.");
		}
		cpp_SuperFishFieldSource->getElectricMagneticField(x,y,z,t,E_x,E_y,E_z,H_x,H_y,H_z);
		return Py_BuildValue("(dddddd)",E_x,E_y,E_z,H_x,H_y,H_z);
	}		
	
	//setGrid2D_Fields(grid2D_Ez,grid2D_Er,grid2D_H) - sets the Grid2D instances with Ez, Er, H fields.
  static PyObject* SuperFishFieldSource_setGrid2D_Fields(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;
		PyObject* pyGrid2D_Ez;
		PyObject* pyGrid2D_Er;
		PyObject* pyGrid2D_H;
		if(!PyArg_ParseTuple(args,"OOO:setGrid2D_Fields",&pyGrid2D_Ez,&pyGrid2D_Er,&pyGrid2D_H)){
			ORBIT_MPI_Finalize("PySuperFishFieldSource - setGrid2D_Fields(grid2D_Ez,grid2D_Er,grid2D_H) - parameters are needed.");
		}		
		PyObject* pyORBIT_Grid2D_Type = getSpaceChargeType("Grid2D");
		if(!PyObject_IsInstance(pyGrid2D_Ez,pyORBIT_Grid2D_Type) || 
			 !PyObject_IsInstance(pyGrid2D_Er,pyORBIT_Grid2D_Type) ||
		   !PyObject_IsInstance(pyGrid2D_H,pyORBIT_Grid2D_Type)) {
			ORBIT_MPI_Finalize("PySuperFishFieldSource - setGrid2D_Fields(grid2D_Ez,grid2D_Er,grid2D_H) - the params. should be Grid2D.");
		}		
		Grid2D* grid2d_Ez = (Grid2D*) ((pyORBIT_Object*)pyGrid2D_Ez)->cpp_obj;
		Grid2D* grid2d_Er = (Grid2D*) ((pyORBIT_Object*)pyGrid2D_Er)->cpp_obj;
		Grid2D* grid2d_H  = (Grid2D*) ((pyORBIT_Object*)pyGrid2D_H)->cpp_obj;
		cpp_SuperFishFieldSource->setGrid2D_Fields(grid2d_Ez,grid2d_Er, grid2d_H);
		Py_INCREF(pyGrid2D_Ez);
		Py_INCREF(pyGrid2D_Er);
		Py_INCREF(pyGrid2D_H);
		Py_INCREF(Py_None);
    return Py_None;				
	}	
	
	//getGrid2D_Fields returns the Grid2D instances with Ez, Er, H fields.
  static PyObject* SuperFishFieldSource_getGrid2D_Fields(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;
		Grid2D* grid2d_Ez = cpp_SuperFishFieldSource->getGrid2D_Ez();
		Grid2D* grid2d_Er = cpp_SuperFishFieldSource->getGrid2D_Er();
		Grid2D* grid2d_H  = cpp_SuperFishFieldSource->getGrid2D_H();
		PyObject* pyGrid2D_Ez = (PyObject*) grid2d_Ez->getPyWrapper();
		PyObject* pyGrid2D_Er = (PyObject*) grid2d_Er->getPyWrapper();
		PyObject* pyGrid2D_H  = (PyObject*) grid2d_H->getPyWrapper();	
		if(pyGrid2D_Ez == NULL || pyGrid2D_Er == NULL || pyGrid2D_H == NULL){
			Py_INCREF(Py_None);
			return Py_None;	
		}
		return Py_BuildValue("(OOO)",pyGrid2D_Ez,pyGrid2D_Er,pyGrid2D_H);		
	}		

	//setFrequency(frequency) sets the RF frequency.
  static PyObject* SuperFishFieldSource_setFrequency(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		double frequency;
		if(!PyArg_ParseTuple(args,"d:setFrequency()",&frequency)){
			ORBIT_MPI_Finalize("PySuperFishFieldSource - setFrequency(frequency) - parameter is needed.");
		}		
		cpp_SuperFishFieldSource->setFrequency(frequency);
		Py_INCREF(Py_None);
		return Py_None;			
	}
	
	//getFrequency() returns the RF frequency.
  static PyObject* SuperFishFieldSource_getFrequency(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		return Py_BuildValue("d",cpp_SuperFishFieldSource->getFrequency());
	}
	
	//setAmplitude(amplitude) sets the RF amplitude.
  static PyObject* SuperFishFieldSource_setAmplitude(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		double amplitude;
		if(!PyArg_ParseTuple(args,"d:setAmplitude()",&amplitude)){
			ORBIT_MPI_Finalize("PySuperFishFieldSource - setAmplitude(amplitude) - parameter is needed.");
		}		
		cpp_SuperFishFieldSource->setAmplitude(amplitude);
		Py_INCREF(Py_None);
		return Py_None;			
	}
		
	//getAmplitude() returns the RF amplitude.
  static PyObject* SuperFishFieldSource_getAmplitude(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		return Py_BuildValue("d",cpp_SuperFishFieldSource->getAmplitude());
	}

	//getAvgField() returns the average e_z field.
  static PyObject* SuperFishFieldSource_getAvgField(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		return Py_BuildValue("d",cpp_SuperFishFieldSource->getAvgField());
	}	
	
	//getLength() returns the length of the field.
  static PyObject* SuperFishFieldSource_getLength(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		return Py_BuildValue("d",cpp_SuperFishFieldSource->getLength());
	}	
		
	//setPhase(phase) sets the RF phase.
  static PyObject* SuperFishFieldSource_setPhase(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		double phase;
		if(!PyArg_ParseTuple(args,"d:setPhase()",&phase)){
			ORBIT_MPI_Finalize("PySuperFishFieldSource - setPhase(phase) - parameter is needed.");
		}		
		cpp_SuperFishFieldSource->setPhase(phase);
		Py_INCREF(Py_None);
		return Py_None;			
	}
		
	//getPhase() returns the RF phase.
  static PyObject* SuperFishFieldSource_getPhase(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		return Py_BuildValue("d",cpp_SuperFishFieldSource->getPhase());
	}	
	
	//setFieldCenterPos(pos) sets the RF center pos.
  static PyObject* SuperFishFieldSource_setFieldCenterPos(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		double pos;
		if(!PyArg_ParseTuple(args,"d:setFieldCenterPos()",&pos)){
			ORBIT_MPI_Finalize("PySuperFishFieldSource - setFieldCenterPos(pos) - parameter is needed.");
		}		
		cpp_SuperFishFieldSource->setFieldCenterPos(pos);
		Py_INCREF(Py_None);
		return Py_None;			
	}
		
	//getFieldCenterPos() returns the RF center pos.
  static PyObject* SuperFishFieldSource_getFieldCenterPos(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		return Py_BuildValue("d",cpp_SuperFishFieldSource->getFieldCenterPos());
	}	
	
	//setDirectionZ(directionZ) sets directionZ = +1 by default and -1 in the inversion z case.
  static PyObject* SuperFishFieldSource_setDirectionZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		int directionZ;
		if(!PyArg_ParseTuple(args,"i:setDirectionZ()",&directionZ)){
			ORBIT_MPI_Finalize("PySuperFishFieldSource - setDirectionZ(directionZ) - parameter is needed.");
		}		
		cpp_SuperFishFieldSource->setDirectionZ(directionZ);
		Py_INCREF(Py_None);
		return Py_None;			
	}
		
	//getDirectionZ() returns the directionZ = +1 by default and -1 in the inversion z case.
  static PyObject* SuperFishFieldSource_getDirectionZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		return Py_BuildValue("i",cpp_SuperFishFieldSource->getDirectionZ());
	}		
	
	//setSymmetry(symm) sets the RF symmetry 0-no symmetry 1-it is there.
  static PyObject* SuperFishFieldSource_setSymmetry(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		int symm;
		if(!PyArg_ParseTuple(args,"i:setSymmetry()",&symm)){
			ORBIT_MPI_Finalize("PySuperFishFieldSource - setSymmetry(symm) - parameter is needed.");
		}		
		cpp_SuperFishFieldSource->setSymmetry(symm);
		Py_INCREF(Py_None);
		return Py_None;			
	}
		
	//getSymmetry() returns the RF symmetry 0-no symmetry 1-it is there.
  static PyObject* SuperFishFieldSource_getSymmetry(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		return Py_BuildValue("i",cpp_SuperFishFieldSource->getSymmetry());
	}	
	
	//setTimeInit(time_init) sets the RF initial time.
  static PyObject* SuperFishFieldSource_setTimeInit(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		double time_init;
		if(!PyArg_ParseTuple(args,"d:setTimeInit()",&time_init)){
			ORBIT_MPI_Finalize("PySuperFishFieldSource - setTimeInit(time_init) - parameter is needed.");
		}		
		cpp_SuperFishFieldSource->setTimeInit(time_init);
		Py_INCREF(Py_None);
		return Py_None;			
	}
		
	//getTimeInit() returns the RF initial time.
  static PyObject* SuperFishFieldSource_getTimeInit(PyObject *self, PyObject *args){
    pyORBIT_Object* pySuperFishFieldSource = (pyORBIT_Object*) self;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) pySuperFishFieldSource->cpp_obj;	
		return Py_BuildValue("d",cpp_SuperFishFieldSource->getTimeInit());
	}	
	
  //-----------------------------------------------------
  //destructor for python SuperFishFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void SuperFishFieldSource_del(pyORBIT_Object* self){
		//std::cerr<<"The SuperFishFieldSource __del__ has been called!"<<std::endl;
		SuperFishFieldSource* cpp_SuperFishFieldSource = (SuperFishFieldSource*) self->cpp_obj;	
		delete cpp_SuperFishFieldSource;
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python SuperFishFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef SuperFishFieldSourceClassMethods[] = {
		{ "setGrid2D_Fields",         SuperFishFieldSource_setGrid2D_Fields,           METH_VARARGS,"sets Grid2D instances with Ez, Er, H fields."},
		{ "getGrid2D_Fields",         SuperFishFieldSource_getGrid2D_Fields,           METH_VARARGS,"returns Grid2D instances with Ez, Er, H fields."},
		{ "getElectricMagneticField", SuperFishFieldSource_getElectricMagneticField,   METH_VARARGS,"returns Ex,Ey,Ez, Bx,By,Bz field components."},
		{ "getFrequency",             SuperFishFieldSource_getFrequency,               METH_VARARGS,"get RF frequency."},
		{ "setFrequency",             SuperFishFieldSource_setFrequency,               METH_VARARGS,"set RF frequency."},
		{ "getAmplitude",             SuperFishFieldSource_getAmplitude,               METH_VARARGS,"get the RF amplitude."},
		{ "setAmplitude",             SuperFishFieldSource_setAmplitude,               METH_VARARGS,"set the RF amplitude."},
		{ "getAvgField",              SuperFishFieldSource_getAvgField,                METH_VARARGS,"get the average e_z filed."},
		{ "getLength",                SuperFishFieldSource_getLength,                  METH_VARARGS,"get the length of the filed."},
		{ "getPhase",                 SuperFishFieldSource_getPhase,                   METH_VARARGS,"get the RF phase."},
		{ "setPhase",                 SuperFishFieldSource_setPhase,                   METH_VARARGS,"set the RF phase."},
		{ "getFieldCenterPos",        SuperFishFieldSource_getFieldCenterPos,          METH_VARARGS,"get the RF center position."},
		{ "setFieldCenterPos",        SuperFishFieldSource_setFieldCenterPos,          METH_VARARGS,"set the RF center position."},
		{ "getDirectionZ",            SuperFishFieldSource_getDirectionZ,              METH_VARARGS,"get the RF directionZ."},
		{ "setDirectionZ",            SuperFishFieldSource_setDirectionZ,              METH_VARARGS,"set the RF directionZ."},		
		{ "getSymmetry",              SuperFishFieldSource_getSymmetry,                METH_VARARGS,"get the RF z symmetry."},
		{ "setSymmetry",              SuperFishFieldSource_setSymmetry,                METH_VARARGS,"set the RF z symmetry."},		
		{ "getTimeInit",              SuperFishFieldSource_getTimeInit,               METH_VARARGS,"get the RF initial time."},
		{ "setTimeInit",              SuperFishFieldSource_setTimeInit,               METH_VARARGS,"set the RF initial time."},
    {NULL}
  };

	// defenition of the memebers of the python SuperFishFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef SuperFishFieldSourceClassMembers [] = {
		{NULL}
	};

	//new python SuperFishFieldSource wrapper type definition
	static PyTypeObject pyORBIT_SuperFishFieldSource_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SuperFishFieldSource", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SuperFishFieldSource_del , /*tp_dealloc*/
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
		"The SuperFishFieldSource python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SuperFishFieldSourceClassMethods, /* tp_methods */
		SuperFishFieldSourceClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SuperFishFieldSource_init, /* tp_init */
		0, /* tp_alloc */
		SuperFishFieldSource_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pySuperFishFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initSuperFishFieldSource(PyObject* module){
		if (PyType_Ready(&pyORBIT_SuperFishFieldSource_Type) < 0) return;
		Py_INCREF(&pyORBIT_SuperFishFieldSource_Type);
		PyModule_AddObject(module, "SuperFishFieldSource", (PyObject *)&pyORBIT_SuperFishFieldSource_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_linac
}
