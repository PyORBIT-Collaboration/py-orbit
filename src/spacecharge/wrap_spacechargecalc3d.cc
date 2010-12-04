#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#
#include "wrap_spacechargecalc3d.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "SpaceChargeCalc3D.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python SpaceChargeCalc3D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping SpaceChargeCalc3D instance
	//It never will be called directly

	static PyObject* SpaceChargeCalc3D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The SpaceChargeCalc3D new has been called!"<<std::endl;
		return (PyObject *) self;
	}
	
  //initializator for python SpaceChargeCalc3D class
  //this is implementation of the __init__ method 
  static int SpaceChargeCalc3D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
  	int xSize,ySize,zSize;
		if(!PyArg_ParseTuple(args,"iii:__init__",&xSize,&ySize,&zSize)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc3D - SpaceChargeCalc3D(xSize,ySize,xzSize) - constructor needs parameters.");
		}
		self->cpp_obj = new SpaceChargeCalc3D(xSize,ySize,zSize);
		((SpaceChargeCalc3D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		return 0;
	}
  
  //Grid3D* getRhoGrid() returns the 3D grid with charge density
  static PyObject* SpaceChargeCalc3D_getRhoGrid(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalc3D = (pyORBIT_Object*) self;
		SpaceChargeCalc3D* cpp_SpaceChargeCalc3D = (SpaceChargeCalc3D*) pySpaceChargeCalc3D->cpp_obj;
		Grid3D* cpp_grid3d = cpp_SpaceChargeCalc3D->getRhoGrid();
		if(cpp_grid3d->getPyWrapper() != NULL){
			Py_INCREF(cpp_grid3d->getPyWrapper());
			return cpp_grid3d->getPyWrapper();
		}
		//It will create a pyGrid3D object
		PyObject* mod = PyImport_ImportModule("spacecharge");
		PyObject* pyGrid3D = PyObject_CallMethod(mod,"Grid3D","iii",cpp_grid3d->getSizeX(),cpp_grid3d->getSizeY(),cpp_grid3d->getSizeZ());		
		//delete the c++ reference to the internal Grid3D inside pyGrid3D and assign the new one
		delete ((Grid3D*)((pyORBIT_Object*) pyGrid3D)->cpp_obj);
		((pyORBIT_Object*) pyGrid3D)->cpp_obj = cpp_grid3d;
		cpp_grid3d->setPyWrapper(pyGrid3D);
		Py_INCREF(cpp_grid3d->getPyWrapper());
		Py_DECREF(mod);
		return pyGrid3D;
  }	
	
  //Grid3D* getPhiGrid() returns the 3D grid with potential
  static PyObject* SpaceChargeCalc3D_getPhiGrid(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalc3D = (pyORBIT_Object*) self;
		SpaceChargeCalc3D* cpp_SpaceChargeCalc3D = (SpaceChargeCalc3D*) pySpaceChargeCalc3D->cpp_obj;
		Grid3D* cpp_grid3d = cpp_SpaceChargeCalc3D->getPhiGrid();
		if(cpp_grid3d->getPyWrapper() != NULL){
			Py_INCREF(cpp_grid3d->getPyWrapper());
			return cpp_grid3d->getPyWrapper();
		}
		//It will create a pyGrid3D object
		PyObject* mod = PyImport_ImportModule("spacecharge");
		PyObject* pyGrid3D = PyObject_CallMethod(mod,"Grid3D","iii",cpp_grid3d->getSizeX(),cpp_grid3d->getSizeY(),cpp_grid3d->getSizeZ());		
		//delete the c++ reference to the internal Grid3D inside pyGrid3D and assign the new one
		delete ((Grid3D*)((pyORBIT_Object*) pyGrid3D)->cpp_obj);
		((pyORBIT_Object*) pyGrid3D)->cpp_obj = cpp_grid3d;
		cpp_grid3d->setPyWrapper(pyGrid3D);
		Py_INCREF(cpp_grid3d->getPyWrapper());
		Py_DECREF(mod);
		return pyGrid3D;
  }		
		
  //trackBunch(Bunch* bunch, double length)
  static PyObject* SpaceChargeCalc3D_trackBunch(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalc3D = (pyORBIT_Object*) self;
		SpaceChargeCalc3D* cpp_SpaceChargeCalc3D = (SpaceChargeCalc3D*) pySpaceChargeCalc3D->cpp_obj;
		PyObject* pyBunch;
		double length;
		if(!PyArg_ParseTuple(args,"Od:trackBunch",&pyBunch,&length)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc3D.trackBunch(pyBunch,length) - method needs parameters.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc3D.trackBunch(pyBunch,length) - pyBunch is not Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;		
		cpp_SpaceChargeCalc3D->trackBunch(cpp_bunch,length);
		Py_INCREF(Py_None);
		return Py_None;
  }

	//setRatioLimit(double ratioLimit) sets the ratio change of x to y and x to z to recalculate Green Functions 
  static PyObject* SpaceChargeCalc3D_setRatioLimit(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalc3D = (pyORBIT_Object*) self;
		SpaceChargeCalc3D* cpp_SpaceChargeCalc3D = (SpaceChargeCalc3D*) pySpaceChargeCalc3D->cpp_obj;
		double ratioLimit;
		if(!PyArg_ParseTuple(args,"d:setRatioLimit",&ratioLimit)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc3D.setRatioLimit(ratioLimit) - method needs a parameter.");
		}
		cpp_SpaceChargeCalc3D->setRatioLimit(ratioLimit);
		return Py_BuildValue("d",ratioLimit);;
  }	
	
	//getRatioLimit() returns the ratio change of x to y and x to z to recalculate Green Functions 
  static PyObject* SpaceChargeCalc3D_getRatioLimit(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalc3D = (pyORBIT_Object*) self;
		SpaceChargeCalc3D* cpp_SpaceChargeCalc3D = (SpaceChargeCalc3D*) pySpaceChargeCalc3D->cpp_obj;
		double ratioLimit = cpp_SpaceChargeCalc3D->getRatioLimit();
		return Py_BuildValue("d",ratioLimit);;
  }		
	
  //-----------------------------------------------------
  //destructor for python SpaceChargeCalc3D class (__del__ method).
  //-----------------------------------------------------
  static void SpaceChargeCalc3D_del(pyORBIT_Object* self){
		SpaceChargeCalc3D* cpp_SpaceChargeCalc3D = (SpaceChargeCalc3D*) self->cpp_obj;
		if(cpp_SpaceChargeCalc3D != NULL){
			delete cpp_SpaceChargeCalc3D;
		}
		self->ob_type->tp_free((PyObject*)self);
  }	
  
  // defenition of the methods of the python SpaceChargeCalc3D wrapper class
  // they will be vailable from python level
  static PyMethodDef SpaceChargeCalc3DClassMethods[] = {
		{ "trackBunch",  SpaceChargeCalc3D_trackBunch, METH_VARARGS,"track the bunch - trackBunch(pyBunch,length,pipe_radius)"},
		{ "getRhoGrid",  SpaceChargeCalc3D_getRhoGrid, METH_VARARGS,"returns the Grid3D with a space charge density"},
		{ "getPhiGrid",  SpaceChargeCalc3D_getPhiGrid, METH_VARARGS,"returns the Grid3D with a space charge potential"},
		{ "setRatioLimit",	SpaceChargeCalc3D_setRatioLimit, METH_VARARGS,"sets the ratio change of x to y and x to z to recalculate Green Functions."},
		{ "getRatioLimit",	SpaceChargeCalc3D_getRatioLimit, METH_VARARGS,"returns the ratio change of x to y and x to z to recalculate Green Functions."},
		{NULL}
  };
  
  // defenition of the memebers of the python SpaceChargeCalc3D wrapper class
  // they will be vailable from python level
  static PyMemberDef SpaceChargeCalc3DClassMembers [] = {
		{NULL}
  };

	//new python SpaceChargeCalc3D wrapper type definition
	static PyTypeObject pyORBIT_SpaceChargeCalc3D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SpaceChargeCalc3D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SpaceChargeCalc3D_del , /*tp_dealloc*/
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
		"The SpaceChargeCalc3D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SpaceChargeCalc3DClassMethods, /* tp_methods */
		SpaceChargeCalc3DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SpaceChargeCalc3D_init, /* tp_init */
		0, /* tp_alloc */
		SpaceChargeCalc3D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pySpaceChargeCalc3D class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initSpaceChargeCalc3D(PyObject* module){
		if (PyType_Ready(&pyORBIT_SpaceChargeCalc3D_Type) < 0) return;
		Py_INCREF(&pyORBIT_SpaceChargeCalc3D_Type);
		PyModule_AddObject(module, "SpaceChargeCalc3D", (PyObject *)&pyORBIT_SpaceChargeCalc3D_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
