#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#
#include "wrap_spacechargecalc2p5d_rb.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "SpaceChargeCalc2p5Drb.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python SpaceChargeCalc2p5Drb class definition
	//---------------------------------------------------------

	//constructor for python class wrapping SpaceChargeCalc2p5Drb instance
	//It never will be called directly

	static PyObject* SpaceChargeCalc2p5Drb_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The SpaceChargeCalc2p5Drb new has been called!"<<std::endl;
		return (PyObject *) self;
	}
	
  //initializator for python SpaceChargeCalc2p5Drb class
  //this is implementation of the __init__ method VSpaceChargeCalc2p5Drb(int xSize, int ySize, int zSize [, double xy_ratio])
  static int SpaceChargeCalc2p5Drb_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
  	int xSize,ySize,zSize;
  	double xy_ratio = -1.0;
		if(!PyArg_ParseTuple(args,"iii|d:__init__",&xSize,&ySize,&zSize,&xy_ratio)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc2p5Drb - SpaceChargeCalc2p5Drb(xSize,ySize,xzSize[,xy_ratio = 1.0]) - constructor needs parameters.");
		}
		if(xy_ratio > 0.){
			self->cpp_obj = new SpaceChargeCalc2p5Drb(xSize,ySize,zSize,xy_ratio);
		} else {
			self->cpp_obj = new SpaceChargeCalc2p5Drb(xSize,ySize,zSize);
		}
		((SpaceChargeCalc2p5Drb*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		return 0;
	}
  
  //Grid2D* getRhoGrid() returns the 2D grid with charge density
  static PyObject* SpaceChargeCalc2p5Drb_getRhoGrid(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalc2p5Drb = (pyORBIT_Object*) self;
		SpaceChargeCalc2p5Drb* cpp_SpaceChargeCalc2p5Drb = (SpaceChargeCalc2p5Drb*) pySpaceChargeCalc2p5Drb->cpp_obj;
		Grid2D* cpp_grid2d = cpp_SpaceChargeCalc2p5Drb->getRhoGrid();
		if(cpp_grid2d->getPyWrapper() != NULL){
			Py_INCREF(cpp_grid2d->getPyWrapper());
			return cpp_grid2d->getPyWrapper();
		}
		//It will create a pyGrid2D object
		PyObject* mod = PyImport_ImportModule("spacecharge");
		PyObject* pyGrid2D = PyObject_CallMethod(mod,"Grid2D","ii",cpp_grid2d->getSizeX(),cpp_grid2d->getSizeY());		
		//delete the c++ reference to the internal Grid2D inside pyGrid2D and assign the new one
		delete ((Grid2D*)((pyORBIT_Object*) pyGrid2D)->cpp_obj);
		((pyORBIT_Object*) pyGrid2D)->cpp_obj = cpp_grid2d;
		cpp_grid2d->setPyWrapper(pyGrid2D);
		Py_INCREF(cpp_grid2d->getPyWrapper());
		Py_DECREF(mod);
		return pyGrid2D;
  }	
	
  //Grid2D* getPhiGrid() returns the 2D grid with potential
  static PyObject* SpaceChargeCalc2p5Drb_getPhiGrid(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalc2p5Drb = (pyORBIT_Object*) self;
		SpaceChargeCalc2p5Drb* cpp_SpaceChargeCalc2p5Drb = (SpaceChargeCalc2p5Drb*) pySpaceChargeCalc2p5Drb->cpp_obj;
		Grid2D* cpp_grid2d = cpp_SpaceChargeCalc2p5Drb->getPhiGrid();
		if(cpp_grid2d->getPyWrapper() != NULL){
			Py_INCREF(cpp_grid2d->getPyWrapper());
			return cpp_grid2d->getPyWrapper();
		}
		//It will create a pyGrid2D object
		PyObject* mod = PyImport_ImportModule("spacecharge");
		PyObject* pyGrid2D = PyObject_CallMethod(mod,"Grid2D","ii",cpp_grid2d->getSizeX(),cpp_grid2d->getSizeY());		
		//delete the c++ reference to the internal Grid2D inside pyGrid2D and assign the new one
		delete ((Grid2D*)((pyORBIT_Object*) pyGrid2D)->cpp_obj);
		((pyORBIT_Object*) pyGrid2D)->cpp_obj = cpp_grid2d;
		cpp_grid2d->setPyWrapper(pyGrid2D);
		Py_INCREF(cpp_grid2d->getPyWrapper());
		Py_DECREF(mod);
		return pyGrid2D;
  }		
	
  //Grid1D* getLongGrid() returns the 1D grid with longitudinal density
  static PyObject* SpaceChargeCalc2p5Drb_getLongGrid(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalc2p5Drb = (pyORBIT_Object*) self;
		SpaceChargeCalc2p5Drb* cpp_SpaceChargeCalc2p5Drb = (SpaceChargeCalc2p5Drb*) pySpaceChargeCalc2p5Drb->cpp_obj;
		Grid1D* cpp_grid1d = cpp_SpaceChargeCalc2p5Drb->getLongGrid();
		if(cpp_grid1d->getPyWrapper() != NULL){
			Py_INCREF(cpp_grid1d->getPyWrapper());
			return cpp_grid1d->getPyWrapper();
		}
		//It will create a pyGrid2D object
		PyObject* mod = PyImport_ImportModule("spacecharge");
		PyObject* pyGrid1D = PyObject_CallMethod(mod,"Grid1D","i",cpp_grid1d->getSizeZ());		
		//delete the c++ reference to the internal Grid1D inside pyGrid1D and assign the new one
		delete ((Grid1D*)((pyORBIT_Object*) pyGrid1D)->cpp_obj);
		((pyORBIT_Object*) pyGrid1D)->cpp_obj = cpp_grid1d;
		cpp_grid1d->setPyWrapper(pyGrid1D);
		Py_INCREF(cpp_grid1d->getPyWrapper());
		Py_DECREF(mod);
		return pyGrid1D;
  }			
	
  //trackBunch(Bunch* bunch, double length, double pipe_radius)
  static PyObject* SpaceChargeCalc2p5Drb_trackBunch(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalc2p5Drb = (pyORBIT_Object*) self;
		SpaceChargeCalc2p5Drb* cpp_SpaceChargeCalc2p5Drb = (SpaceChargeCalc2p5Drb*) pySpaceChargeCalc2p5Drb->cpp_obj;
		PyObject* pyBunch;
		double length, pipe_radius;
		if(!PyArg_ParseTuple(args,"Odd:trackBunch",&pyBunch,&length,&pipe_radius)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc2p5Drb.trackBunch(pyBunch,length,pipe_radius) - method needs parameters.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc2p5Drb.trackBunch(pyBunch,length,pipe_radius) - pyBunch is not Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;		
		cpp_SpaceChargeCalc2p5Drb->trackBunch(cpp_bunch,length,pipe_radius);
		Py_INCREF(Py_None);
		return Py_None;
  }

  //-----------------------------------------------------
  //destructor for python SpaceChargeCalc2p5Drb class (__del__ method).
  //-----------------------------------------------------
  static void SpaceChargeCalc2p5Drb_del(pyORBIT_Object* self){
		SpaceChargeCalc2p5Drb* cpp_SpaceChargeCalc2p5Drb = (SpaceChargeCalc2p5Drb*) self->cpp_obj;
		if(cpp_SpaceChargeCalc2p5Drb != NULL){
			delete cpp_SpaceChargeCalc2p5Drb;
		}
		self->ob_type->tp_free((PyObject*)self);
  }	
  
  // defenition of the methods of the python SpaceChargeCalc2p5Drb wrapper class
  // they will be vailable from python level
  static PyMethodDef SpaceChargeCalc2p5DrbClassMethods[] = {
		{ "trackBunch",  SpaceChargeCalc2p5Drb_trackBunch, METH_VARARGS,"track the bunch - trackBunch(pyBunch,length,pipe_radius)"},
		{ "getRhoGrid",  SpaceChargeCalc2p5Drb_getRhoGrid, METH_VARARGS,"returns the Grid2D with a space charge density"},
		{ "getPhiGrid",  SpaceChargeCalc2p5Drb_getPhiGrid, METH_VARARGS,"returns the Grid2D with a space charge potential"},
		{ "getLongGrid", SpaceChargeCalc2p5Drb_getLongGrid, METH_VARARGS,"returns the Grid1D with a longitudinal space charge density"},
		{NULL}
  };
  
  // defenition of the memebers of the python SpaceChargeCalc2p5Drb wrapper class
  // they will be vailable from python level
  static PyMemberDef SpaceChargeCalc2p5DrbClassMembers [] = {
		{NULL}
  };

	//new python SpaceChargeCalc2p5Drb wrapper type definition
	static PyTypeObject pyORBIT_SpaceChargeCalc2p5Drb_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SpaceChargeCalc2p5Drb", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SpaceChargeCalc2p5Drb_del , /*tp_dealloc*/
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
		"The SpaceChargeCalc2p5Drb python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SpaceChargeCalc2p5DrbClassMethods, /* tp_methods */
		SpaceChargeCalc2p5DrbClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SpaceChargeCalc2p5Drb_init, /* tp_init */
		0, /* tp_alloc */
		SpaceChargeCalc2p5Drb_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pySpaceChargeCalc2p5Drb class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initSpaceChargeCalc2p5Drb(PyObject* module){
		if (PyType_Ready(&pyORBIT_SpaceChargeCalc2p5Drb_Type) < 0) return;
		Py_INCREF(&pyORBIT_SpaceChargeCalc2p5Drb_Type);
		PyModule_AddObject(module, "SpaceChargeCalc2p5Drb", (PyObject *)&pyORBIT_SpaceChargeCalc2p5Drb_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
