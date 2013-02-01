#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "ForceSolverFFT2D.hh"
#include "Grid2D.hh"

#include "wrap_forcesolverfft2d.hh"
#include "wrap_spacecharge.hh"

#include <iostream>

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python forceSolverFFT2D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping forceSolverFFT2D instance
	//It never will be called directly
	static PyObject* ForceSolverFFT2D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The ForceSolverFFT2D new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  ForceSolverFFT2D class
  //this is implementation of the __init__ method
  static int ForceSolverFFT2D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		int xSize, ySize;
		if(!PyArg_ParseTuple(args,"ii:__init__",&xSize,&ySize)){
			ORBIT_MPI_Finalize("PyForceSolverFFT2D - ForceSolverFFT2D(nX,nY) - constructor needs parameters.");
		}	
		self->cpp_obj = new ForceSolverFFT2D(xSize,ySize);
		((ForceSolverFFT2D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The ForceSolverFFT2D __init__ has been called!"<<std::endl;
		return 0;
  }
		
	//get grid size in X 
	static PyObject* ForceSolverFFT2D_getSizeX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyForceSolverFFT2D = (pyORBIT_Object*) self;
		ForceSolverFFT2D* cpp_ForceSolverFFT2D = (ForceSolverFFT2D*) pyForceSolverFFT2D->cpp_obj;
		return Py_BuildValue("i",cpp_ForceSolverFFT2D->getSizeX());
	}
	
	//get grid size in Y 
	static PyObject* ForceSolverFFT2D_getSizeY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyForceSolverFFT2D = (pyORBIT_Object*) self;
		ForceSolverFFT2D* cpp_ForceSolverFFT2D = (ForceSolverFFT2D*) pyForceSolverFFT2D->cpp_obj;
		return Py_BuildValue("i",cpp_ForceSolverFFT2D->getSizeY());
	}
	
	// getMaxX()
	static PyObject* ForceSolverFFT2D_getMaxX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyForceSolverFFT2D = (pyORBIT_Object*) self;
		ForceSolverFFT2D* cpp_ForceSolverFFT2D = (ForceSolverFFT2D*) pyForceSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_ForceSolverFFT2D->getMaxX());
	}
	
	// getMaxY()
	static PyObject* ForceSolverFFT2D_getMaxY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyForceSolverFFT2D = (pyORBIT_Object*) self;
		ForceSolverFFT2D* cpp_ForceSolverFFT2D = (ForceSolverFFT2D*) pyForceSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_ForceSolverFFT2D->getMaxY());
	}
	
	// getMinX()
	static PyObject* ForceSolverFFT2D_getMinX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyForceSolverFFT2D = (pyORBIT_Object*) self;
		ForceSolverFFT2D* cpp_ForceSolverFFT2D = (ForceSolverFFT2D*) pyForceSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_ForceSolverFFT2D->getMinX());
	}
	
	// getMinY()
	static PyObject* ForceSolverFFT2D_getMinY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyForceSolverFFT2D = (pyORBIT_Object*) self;
		ForceSolverFFT2D* cpp_ForceSolverFFT2D = (ForceSolverFFT2D*) pyForceSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_ForceSolverFFT2D->getMinY());
	}	
	
	// getStepX()
	static PyObject* ForceSolverFFT2D_getStepX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyForceSolverFFT2D = (pyORBIT_Object*) self;
		ForceSolverFFT2D* cpp_ForceSolverFFT2D = (ForceSolverFFT2D*) pyForceSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_ForceSolverFFT2D->getStepX());
	}	
	
	// getStepY()
	static PyObject* ForceSolverFFT2D_getStepY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyForceSolverFFT2D = (pyORBIT_Object*) self;
		ForceSolverFFT2D* cpp_ForceSolverFFT2D = (ForceSolverFFT2D*) pyForceSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_ForceSolverFFT2D->getStepY());
	}			
		
	//findForce(Grid2D* rhoGrid2D,Grid2D* forceGridX, Grid2D* forceGridY)
  static PyObject* ForceSolverFFT2D_findForce(PyObject *self, PyObject *args){
    pyORBIT_Object* pyForceSolverFFT2D = (pyORBIT_Object*) self;
		ForceSolverFFT2D* cpp_ForceSolverFFT2D = (ForceSolverFFT2D*) pyForceSolverFFT2D->cpp_obj;
		PyObject* pyRhoG;
		PyObject* pyForceGX;
		PyObject* pyForceGY;
		if(!PyArg_ParseTuple(args,"OOO:__init__",&pyRhoG,&pyForceGX,&pyForceGY)){
			ORBIT_MPI_Finalize("PyForceSolverFFT2D.findForce(Grid2D rhoGrid2D,Grid2D forceXGrid2D, forceYGrid2D) - method needs parameters.");
		}	
		PyObject* pyORBIT_Grid2D_Type = getSpaceChargeType("Grid2D");
		if(!PyObject_IsInstance(pyRhoG,pyORBIT_Grid2D_Type) || !PyObject_IsInstance(pyForceGX,pyORBIT_Grid2D_Type) || !PyObject_IsInstance(pyForceGY,pyORBIT_Grid2D_Type)){
			ORBIT_MPI_Finalize("PyForceSolverFFT2D.findForce(Grid2D rhoGrid2D,Grid2D forceXGrid2D, forceYGrid2D) - method needs parameters.");
		}			
		Grid2D* grid2D_rho = (Grid2D*)(((pyORBIT_Object*) pyRhoG)->cpp_obj);
		Grid2D* grid2D_xforce = (Grid2D*)(((pyORBIT_Object*) pyForceGX)->cpp_obj);
		Grid2D* grid2D_yforce = (Grid2D*)(((pyORBIT_Object*) pyForceGY)->cpp_obj);

		cpp_ForceSolverFFT2D->findForce(grid2D_rho,grid2D_xforce,grid2D_yforce);
		Py_INCREF(Py_None);
    return Py_None;
	}		
	
  //-----------------------------------------------------
  //destructor for python ForceSolverFFT2D class (__del__ method).
  //-----------------------------------------------------
  static void ForceSolverFFT2D_del(pyORBIT_Object* self){
	  //std::cerr<<"The ForceSolverFFT2D __del__ has been called! 0"<<std::endl;
		ForceSolverFFT2D* cpp_ForceSolverFFT2D = (ForceSolverFFT2D*) self->cpp_obj;
		if(cpp_ForceSolverFFT2D != NULL){
			delete cpp_ForceSolverFFT2D;
		}
		self->ob_type->tp_free((PyObject*)self);
		//std::cerr<<"The ForceSolverFFT2D __del__ has been called! 1"<<std::endl;
  }
	
	// defenition of the methods of the python ForceSolverFFT2D wrapper class
	// they will be vailable from python level
  static PyMethodDef ForceSolverFFT2DClassMethods[] = {
		{ "getSizeX",            ForceSolverFFT2D_getSizeX,            METH_VARARGS,"returns the number of grid points in x-direction"},
		{ "getSizeY",            ForceSolverFFT2D_getSizeY,            METH_VARARGS,"returns the number of grid points in y-direction"},
		{ "getMaxX",             ForceSolverFFT2D_getMaxX,             METH_VARARGS,"returns max grid value in x-direction"},
		{ "getMaxY",             ForceSolverFFT2D_getMaxY,             METH_VARARGS,"returns max grid value in y-direction"},
		{ "getMinX",             ForceSolverFFT2D_getMinX,             METH_VARARGS,"returns min grid value in x-direction"},
		{ "getMinY",             ForceSolverFFT2D_getMinY,             METH_VARARGS,"returns min grid value in y-direction"},
		{ "getStepX",            ForceSolverFFT2D_getStepX,            METH_VARARGS,"returns grid step in x-direction"},
		{ "getStepY",            ForceSolverFFT2D_getStepY,            METH_VARARGS,"returns grid step in y-direction"},
		{ "findForce",       ForceSolverFFT2D_findForce,       METH_VARARGS,"findForce(Grid2D rhoGrid2D,Grid2D forceXGrid2D, Grid2D forceYGrid2D)"},
    {NULL}
  };

	// defenition of the memebers of the python ForceSolverFFT2D wrapper class
	// they will be vailable from python level
	static PyMemberDef ForceSolverFFT2DClassMembers [] = {
		{NULL}
	};

	//new python ForceSolverFFT2D wrapper type definition
	static PyTypeObject pyORBIT_ForceSolverFFT2D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"ForceSolverFFT2D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) ForceSolverFFT2D_del , /*tp_dealloc*/
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
		"The ForceSolverFFT2D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		ForceSolverFFT2DClassMethods, /* tp_methods */
		ForceSolverFFT2DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) ForceSolverFFT2D_init, /* tp_init */
		0, /* tp_alloc */
		ForceSolverFFT2D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyForceSolverFFT2D class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initForceSolverFFT2D(PyObject* module){
		if (PyType_Ready(&pyORBIT_ForceSolverFFT2D_Type) < 0) return;
		Py_INCREF(&pyORBIT_ForceSolverFFT2D_Type);
		PyModule_AddObject(module, "ForceSolverFFT2D", (PyObject *)&pyORBIT_ForceSolverFFT2D_Type);
		//std::cout<<"debug ForceSolverFFT2D added! "<<std::endl;
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
