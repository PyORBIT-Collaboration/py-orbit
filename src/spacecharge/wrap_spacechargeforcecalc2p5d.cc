#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#
#include "wrap_spacechargeforcecalc2p5d.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "SpaceChargeForceCalc2p5D.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python SpaceChargeForceCalc2p5D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping SpaceChargeForceCalc2p5D instance
	//It never will be called directly

	static PyObject* SpaceChargeForceCalc2p5D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The SpaceChargeForceCalc2p5D new has been called!"<<std::endl;
		return (PyObject *) self;
	}
	
  //initializator for python SpaceChargeForceCalc2p5D class
  //this is implementation of the __init__ method VSpaceChargeForceCalc2p5D(int xSize, int ySize, int zSize)
  static int SpaceChargeForceCalc2p5D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
  	int xSize,ySize,zSize;
		if(!PyArg_ParseTuple(args,"iii:__init__",&xSize,&ySize,&zSize)){
			ORBIT_MPI_Finalize("PySpaceChargeForceCalc2p5D - SpaceChargeForceCalc2p5D(xSize,ySize,zSize) - constructor needs parameters.");
		}
		self->cpp_obj = new SpaceChargeForceCalc2p5D(xSize,ySize,zSize);
		
		((SpaceChargeForceCalc2p5D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The SpaceChargeForceCalc2p5D __init__ has been called!"<<std::endl;
		return 0;
	}
  
  //Grid2D* getRhoGrid() returns the 2D grid with charge density
  static PyObject* SpaceChargeForceCalc2p5D_getRhoGrid(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeForceCalc2p5D = (pyORBIT_Object*) self;
		SpaceChargeForceCalc2p5D* cpp_SpaceChargeForceCalc2p5D = (SpaceChargeForceCalc2p5D*) pySpaceChargeForceCalc2p5D->cpp_obj;
		Grid2D* cpp_grid2d = cpp_SpaceChargeForceCalc2p5D->getRhoGrid();
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
	
  //Grid2D* getPhiGrid() returns the 2D grid with horizontal force
  static PyObject* SpaceChargeForceCalc2p5D_getForceGridX(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeForceCalc2p5D = (pyORBIT_Object*) self;
		SpaceChargeForceCalc2p5D* cpp_SpaceChargeForceCalc2p5D = (SpaceChargeForceCalc2p5D*) pySpaceChargeForceCalc2p5D->cpp_obj;
		Grid2D* cpp_grid2d = cpp_SpaceChargeForceCalc2p5D->getForceGridX();
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
	
	
	//Grid2D* getForceGridX() returns the 2D grid with horizontal force
	static PyObject* SpaceChargeForceCalc2p5D_getForceGridY(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeForceCalc2p5D = (pyORBIT_Object*) self;
		SpaceChargeForceCalc2p5D* cpp_SpaceChargeForceCalc2p5D = (SpaceChargeForceCalc2p5D*) pySpaceChargeForceCalc2p5D->cpp_obj;
		Grid2D* cpp_grid2d = cpp_SpaceChargeForceCalc2p5D->getForceGridY();
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
  static PyObject* SpaceChargeForceCalc2p5D_getLongGrid(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeForceCalc2p5D = (pyORBIT_Object*) self;
		SpaceChargeForceCalc2p5D* cpp_SpaceChargeForceCalc2p5D = (SpaceChargeForceCalc2p5D*) pySpaceChargeForceCalc2p5D->cpp_obj;
		Grid1D* cpp_grid1d = cpp_SpaceChargeForceCalc2p5D->getLongGrid();
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
	
    //trackBunch(Bunch* bunch, double length[,BaseBoundary2D* boundary])
  static PyObject* SpaceChargeForceCalc2p5D_trackBunch(PyObject *self, PyObject *args){
		int nVars = PyTuple_Size(args);
		pyORBIT_Object* pySpaceChargeForceCalc2p5D = (pyORBIT_Object*) self;
		SpaceChargeForceCalc2p5D* cpp_SpaceChargeForceCalc2p5D = (SpaceChargeForceCalc2p5D*) pySpaceChargeForceCalc2p5D->cpp_obj;
		PyObject* pyBunch;
		double length;
		
		if(!PyArg_ParseTuple(args,"Od:trackBunch",&pyBunch,&length)){
			ORBIT_MPI_Finalize("PySpaceChargeForceCalc2p5D.trackBunch(pyBunch,length - method needs parameters.");
		    }
		    PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		    if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PySpaceChargeForceCalc2p5D.trackBunch(pyBunch,length) - pyBunch is not Bunch.");
		    }
		    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		    cpp_SpaceChargeForceCalc2p5D->trackBunch(cpp_bunch,length);
		  
	
		Py_INCREF(Py_None);
		return Py_None;  
  }

  //-----------------------------------------------------
  //destructor for python SpaceChargeForceCalc2p5D class (__del__ method).
  //-----------------------------------------------------
  static void SpaceChargeForceCalc2p5D_del(pyORBIT_Object* self){
		SpaceChargeForceCalc2p5D* cpp_SpaceChargeForceCalc2p5D = (SpaceChargeForceCalc2p5D*) self->cpp_obj;
		if(cpp_SpaceChargeForceCalc2p5D != NULL){
			delete cpp_SpaceChargeForceCalc2p5D;
		}
		self->ob_type->tp_free((PyObject*)self);
  }	
  
  // defenition of the methods of the python SpaceChargeForceCalc2p5D wrapper class
  // they will be vailable from python level
  static PyMethodDef SpaceChargeForceCalc2p5DClassMethods[] = {
		{ "trackBunch",  SpaceChargeForceCalc2p5D_trackBunch, METH_VARARGS,"track the bunch - trackBunch(pyBunch,length,boundary)"},
		{ "getRhoGrid",  SpaceChargeForceCalc2p5D_getRhoGrid, METH_VARARGS,"returns the Grid2D with a space charge density"},
		{ "getForceGridX",  SpaceChargeForceCalc2p5D_getForceGridX, METH_VARARGS,"returns the Grid2D with a space charge horizontal force"},
		{ "getForceGridY",  SpaceChargeForceCalc2p5D_getForceGridY, METH_VARARGS,"returns the Grid2D with a space charge horizontal force"},
		{ "getLongGrid", SpaceChargeForceCalc2p5D_getLongGrid, METH_VARARGS,"returns the Grid1D with a longitudinal space charge density"},
		{NULL}
  };
  
  // defenition of the memebers of the python SpaceChargeForceCalc2p5D wrapper class
  // they will be vailable from python level
  static PyMemberDef SpaceChargeForceCalc2p5DClassMembers [] = {
		{NULL}
  };

	//new python SpaceChargeForceCalc2p5D wrapper type definition
	static PyTypeObject pyORBIT_SpaceChargeForceCalc2p5D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SpaceChargeForceCalc2p5D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SpaceChargeForceCalc2p5D_del , /*tp_dealloc*/
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
		"The SpaceChargeForceCalc2p5D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SpaceChargeForceCalc2p5DClassMethods, /* tp_methods */
		SpaceChargeForceCalc2p5DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SpaceChargeForceCalc2p5D_init, /* tp_init */
		0, /* tp_alloc */
		SpaceChargeForceCalc2p5D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pySpaceChargeForceCalc2p5D class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initSpaceChargeForceCalc2p5D(PyObject* module){
		if (PyType_Ready(&pyORBIT_SpaceChargeForceCalc2p5D_Type) < 0) return;
		Py_INCREF(&pyORBIT_SpaceChargeForceCalc2p5D_Type);
		PyModule_AddObject(module, "SpaceChargeForceCalc2p5D", (PyObject *)&pyORBIT_SpaceChargeForceCalc2p5D_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
