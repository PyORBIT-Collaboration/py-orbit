#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#
#include "wrap_spacechargecalc_slicebyslice_2D.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "SpaceChargeCalcSliceBySlice2D.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python SpaceChargeCalcSliceBySlice2D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping SpaceChargeCalcSliceBySlice2D instance
	//It never will be called directly

	static PyObject* SpaceChargeCalcSliceBySlice2D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The SpaceChargeCalcSliceBySlice2D new has been called!"<<std::endl;
		return (PyObject *) self;
	}
	
  //initializator for python SpaceChargeCalcSliceBySlice2D class
  //this is implementation of the __init__ method VSpaceChargeCalcSliceBySlice2D(int xSize, int ySize, int zSize [, double xy_ratio_in])
  static int SpaceChargeCalcSliceBySlice2D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
  	int xSize,ySize,zSize;
  	double xy_ratio = -1.0;
		if(!PyArg_ParseTuple(args,"iii|d:__init__",&xSize,&ySize,&zSize,&xy_ratio)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcSliceBySlice2D - SpaceChargeCalcSliceBySlice2D(xSize,ySize,xzSize[,xy_ratio = 1.0]) - constructor needs parameters.");
		}
		if(xy_ratio > 0.){
			self->cpp_obj = new SpaceChargeCalcSliceBySlice2D(xSize,ySize,zSize,xy_ratio);
		} else {
			self->cpp_obj = new SpaceChargeCalcSliceBySlice2D(xSize,ySize,zSize);
		}
		((SpaceChargeCalcSliceBySlice2D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The SpaceChargeCalcSliceBySlice2D __init__ has been called!"<<std::endl;
		return 0;
	}
  
  //Grid3D* getRhoGrid() returns the 3D grid with charge density
  static PyObject* SpaceChargeCalcSliceBySlice2D_getRhoGrid(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalcSliceBySlice2D = (pyORBIT_Object*) self;
		SpaceChargeCalcSliceBySlice2D* cpp_SpaceChargeCalcSliceBySlice2D = (SpaceChargeCalcSliceBySlice2D*) pySpaceChargeCalcSliceBySlice2D->cpp_obj;
		Grid3D* cpp_grid3d = cpp_SpaceChargeCalcSliceBySlice2D->getRhoGrid();
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
  static PyObject* SpaceChargeCalcSliceBySlice2D_getPhiGrid(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalcSliceBySlice2D = (pyORBIT_Object*) self;
		SpaceChargeCalcSliceBySlice2D* cpp_SpaceChargeCalcSliceBySlice2D = (SpaceChargeCalcSliceBySlice2D*) pySpaceChargeCalcSliceBySlice2D->cpp_obj;
		Grid3D* cpp_grid3d = cpp_SpaceChargeCalcSliceBySlice2D->getPhiGrid();
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
  	
  //trackBunch(Bunch* bunch, double length[,BaseBoundary2D* boundary])
  static PyObject* SpaceChargeCalcSliceBySlice2D_trackBunch(PyObject *self, PyObject *args){
		int nVars = PyTuple_Size(args);
		pyORBIT_Object* pySpaceChargeCalcSliceBySlice2D = (pyORBIT_Object*) self;
		SpaceChargeCalcSliceBySlice2D* cpp_SpaceChargeCalcSliceBySlice2D = (SpaceChargeCalcSliceBySlice2D*) pySpaceChargeCalcSliceBySlice2D->cpp_obj;
		PyObject* pyBunch;
		PyObject* pyBoundary;
		double length;
		
		if(nVars == 2 ||  nVars == 3){
		  if (nVars == 2){
		    if(!PyArg_ParseTuple(args,"Od:trackBunch",&pyBunch,&length)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcSliceBySlice2D.trackBunch(pyBunch,length,NULL) - method needs parameters.");
		    }
		    PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		    if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcSliceBySlice2D.trackBunch(pyBunch,length,NULL) - pyBunch is not Bunch.");
		    }
		    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;		
		    cpp_SpaceChargeCalcSliceBySlice2D->trackBunch(cpp_bunch,length,NULL);
		  }
		  else{		    
		    if(!PyArg_ParseTuple(args,"OdO:trackBunch",&pyBunch,&length,&pyBoundary)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcSliceBySlice2D.trackBunch(pyBunch,length,boundary) - method needs parameters.");
		    }
		    PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		    PyObject* pyORBIT_Boundary_Type = getSpaceChargeType("Boundary2D");
		    if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type) || !PyObject_IsInstance(pyBoundary,pyORBIT_Boundary_Type)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcSliceBySlice2D.trackBunch(pyBunch,length,boundary) - pyBunch is not Bunch.");
		    }
		    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;		
		    BaseBoundary2D* cpp_boundary = (BaseBoundary2D*) ((pyORBIT_Object*)pyBoundary)->cpp_obj;
		    cpp_SpaceChargeCalcSliceBySlice2D->trackBunch(cpp_bunch,length,cpp_boundary);	
		    //std::cerr<<"The boundary has been called!"<<std::endl;
		  }
		}
		else{
		  ORBIT_MPI_Finalize("PyBoundary. You should call trackBunch(pyBunch,length) or trackBunch(pyBunch,length,boundary)");  
		}		
		Py_INCREF(Py_None);
		return Py_None;  
  }

  //-----------------------------------------------------
  //destructor for python SpaceChargeCalcSliceBySlice2D class (__del__ method).
  //-----------------------------------------------------
  static void SpaceChargeCalcSliceBySlice2D_del(pyORBIT_Object* self){
		SpaceChargeCalcSliceBySlice2D* cpp_SpaceChargeCalcSliceBySlice2D = (SpaceChargeCalcSliceBySlice2D*) self->cpp_obj;
		if(cpp_SpaceChargeCalcSliceBySlice2D != NULL){
			delete cpp_SpaceChargeCalcSliceBySlice2D;
		}
		self->ob_type->tp_free((PyObject*)self);
  }	
  
  // defenition of the methods of the python SpaceChargeCalcSliceBySlice2D wrapper class
  // they will be vailable from python level
  static PyMethodDef SpaceChargeCalcSliceBySlice2DClassMethods[] = {
		{ "trackBunch",  SpaceChargeCalcSliceBySlice2D_trackBunch, METH_VARARGS,"track the bunch - trackBunch(pyBunch,length,boundary)"},
		{ "getRhoGrid",  SpaceChargeCalcSliceBySlice2D_getRhoGrid, METH_VARARGS,"returns the Grid3D with a space charge density"},
		{ "getPhiGrid",  SpaceChargeCalcSliceBySlice2D_getPhiGrid, METH_VARARGS,"returns the Grid3D with a space charge potential"},
		{NULL}
  };
  
  // defenition of the memebers of the python SpaceChargeCalcSliceBySlice2D wrapper class
  // they will be vailable from python level
  static PyMemberDef SpaceChargeCalcSliceBySlice2DClassMembers [] = {
		{NULL}
  };

	//new python SpaceChargeCalcSliceBySlice2D wrapper type definition
	static PyTypeObject pyORBIT_SpaceChargeCalcSliceBySlice2D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SpaceChargeCalcSliceBySlice2D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SpaceChargeCalcSliceBySlice2D_del , /*tp_dealloc*/
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
		"The SpaceChargeCalcSliceBySlice2D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SpaceChargeCalcSliceBySlice2DClassMethods, /* tp_methods */
		SpaceChargeCalcSliceBySlice2DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SpaceChargeCalcSliceBySlice2D_init, /* tp_init */
		0, /* tp_alloc */
		SpaceChargeCalcSliceBySlice2D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pySpaceChargeCalcSliceBySlice2D class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initSpaceChargeCalcSliceBySlice2D(PyObject* module){
		if (PyType_Ready(&pyORBIT_SpaceChargeCalcSliceBySlice2D_Type) < 0) return;
		Py_INCREF(&pyORBIT_SpaceChargeCalcSliceBySlice2D_Type);
		PyModule_AddObject(module, "SpaceChargeCalcSliceBySlice2D", (PyObject *)&pyORBIT_SpaceChargeCalcSliceBySlice2D_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
