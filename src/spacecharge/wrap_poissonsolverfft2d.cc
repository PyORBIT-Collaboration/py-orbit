#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "PoissonSolverFFT2D.hh"
#include "Grid2D.hh"

#include "wrap_poissonsolverfft2d.hh"
#include "wrap_spacecharge.hh"

#include <iostream>

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python PoissonSolverFFT2D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping PoissonSolverFFT2D instance
	//It never will be called directly
	static PyObject* PoissonSolverFFT2D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The PoissonSolverFFT2D new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  PoissonSolverFFT2D class
  //this is implementation of the __init__ method
  static int PoissonSolverFFT2D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		int xSize, ySize;
		double xMin = -1.0, xMax = +1.0;
		double yMin = -1.0, yMax = +1.0;
		if(!PyArg_ParseTuple(args,"ii|dddd:__init__",&xSize,&ySize,&xMin,&xMax,&yMin,&yMax)){
			ORBIT_MPI_Finalize("PyPoissonSolverFFT2D - PoissonSolverFFT2D(nX,nY[,xMin,xMax,yMin,yMax]) - constructor needs parameters.");
		}	
		self->cpp_obj = new PoissonSolverFFT2D(xSize,ySize,xMin,xMax,yMin,yMax);
		((PoissonSolverFFT2D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The PoissonSolverFFT2D __init__ has been called!"<<std::endl;
		return 0;
  }
		
	//get grid size in X 
	static PyObject* PoissonSolverFFT2D_getSizeX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT2D = (pyORBIT_Object*) self;
		PoissonSolverFFT2D* cpp_PoissonSolverFFT2D = (PoissonSolverFFT2D*) pyPoissonSolverFFT2D->cpp_obj;
		return Py_BuildValue("i",cpp_PoissonSolverFFT2D->getSizeX());
	}
	
	//get grid size in Y 
	static PyObject* PoissonSolverFFT2D_getSizeY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT2D = (pyORBIT_Object*) self;
		PoissonSolverFFT2D* cpp_PoissonSolverFFT2D = (PoissonSolverFFT2D*) pyPoissonSolverFFT2D->cpp_obj;
		return Py_BuildValue("i",cpp_PoissonSolverFFT2D->getSizeY());
	}
	
	// getMaxX()
	static PyObject* PoissonSolverFFT2D_getMaxX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT2D = (pyORBIT_Object*) self;
		PoissonSolverFFT2D* cpp_PoissonSolverFFT2D = (PoissonSolverFFT2D*) pyPoissonSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT2D->getMaxX());
	}
	
	// getMaxY()
	static PyObject* PoissonSolverFFT2D_getMaxY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT2D = (pyORBIT_Object*) self;
		PoissonSolverFFT2D* cpp_PoissonSolverFFT2D = (PoissonSolverFFT2D*) pyPoissonSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT2D->getMaxY());
	}
	
	// getMinX()
	static PyObject* PoissonSolverFFT2D_getMinX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT2D = (pyORBIT_Object*) self;
		PoissonSolverFFT2D* cpp_PoissonSolverFFT2D = (PoissonSolverFFT2D*) pyPoissonSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT2D->getMinX());
	}
	
	// getMinY()
	static PyObject* PoissonSolverFFT2D_getMinY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT2D = (pyORBIT_Object*) self;
		PoissonSolverFFT2D* cpp_PoissonSolverFFT2D = (PoissonSolverFFT2D*) pyPoissonSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT2D->getMinY());
	}	
	
	// getStepX()
	static PyObject* PoissonSolverFFT2D_getStepX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT2D = (pyORBIT_Object*) self;
		PoissonSolverFFT2D* cpp_PoissonSolverFFT2D = (PoissonSolverFFT2D*) pyPoissonSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT2D->getStepX());
	}	
	
	// getStepY()
	static PyObject* PoissonSolverFFT2D_getStepY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT2D = (pyORBIT_Object*) self;
		PoissonSolverFFT2D* cpp_PoissonSolverFFT2D = (PoissonSolverFFT2D*) pyPoissonSolverFFT2D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT2D->getStepY());
	}			
		
	//findPotential(Grid2D* rhoGrid2D,Grid2D* phiGrid2D)
  static PyObject* PoissonSolverFFT2D_findPotential(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPoissonSolverFFT2D = (pyORBIT_Object*) self;
		PoissonSolverFFT2D* cpp_PoissonSolverFFT2D = (PoissonSolverFFT2D*) pyPoissonSolverFFT2D->cpp_obj;
		PyObject* pyRhoG;
		PyObject* pyPhiG;
		if(!PyArg_ParseTuple(args,"OO:__init__",&pyRhoG,&pyPhiG)){		
			ORBIT_MPI_Finalize("PyPoissonSolverFFT2D.findPotential(Grid2D rhoGrid2D,Grid2D phiGrid2D) - method needs parameters.");
		}	
		PyObject* pyORBIT_Grid2D_Type = getSpaceChargeType("Grid2D");
		if(!PyObject_IsInstance(pyRhoG,pyORBIT_Grid2D_Type) || !PyObject_IsInstance(pyPhiG,pyORBIT_Grid2D_Type)){
			ORBIT_MPI_Finalize("PyPoissonSolverFFT2D.findPotential(Grid2D rhoGrid2D,Grid2D phiGrid2D) - method needs parameters.");
		}			
		Grid2D* grid2D_rho = (Grid2D*)(((pyORBIT_Object*) pyRhoG)->cpp_obj);
		Grid2D* grid2D_phi = (Grid2D*)(((pyORBIT_Object*) pyPhiG)->cpp_obj);
		cpp_PoissonSolverFFT2D->findPotential(grid2D_rho,grid2D_phi);
		Py_INCREF(Py_None);
    return Py_None;
	}		
	
  //-----------------------------------------------------
  //destructor for python PoissonSolverFFT2D class (__del__ method).
  //-----------------------------------------------------
  static void PoissonSolverFFT2D_del(pyORBIT_Object* self){
	  //std::cerr<<"The PoissonSolverFFT2D __del__ has been called! 0"<<std::endl;
		PoissonSolverFFT2D* cpp_PoissonSolverFFT2D = (PoissonSolverFFT2D*) self->cpp_obj;
		if(cpp_PoissonSolverFFT2D != NULL){
			delete cpp_PoissonSolverFFT2D;
		}
		self->ob_type->tp_free((PyObject*)self);
		//std::cerr<<"The PoissonSolverFFT2D __del__ has been called! 1"<<std::endl;
  }
	
	// defenition of the methods of the python PoissonSolverFFT2D wrapper class
	// they will be vailable from python level
  static PyMethodDef PoissonSolverFFT2DClassMethods[] = {
		{ "getSizeX",            PoissonSolverFFT2D_getSizeX,            METH_VARARGS,"returns the number of grid points in x-direction"},
		{ "getSizeY",            PoissonSolverFFT2D_getSizeY,            METH_VARARGS,"returns the number of grid points in y-direction"},
		{ "getMaxX",             PoissonSolverFFT2D_getMaxX,             METH_VARARGS,"returns max grid value in x-direction"},
		{ "getMaxY",             PoissonSolverFFT2D_getMaxY,             METH_VARARGS,"returns max grid value in y-direction"},
		{ "getMinX",             PoissonSolverFFT2D_getMinX,             METH_VARARGS,"returns min grid value in x-direction"},
		{ "getMinY",             PoissonSolverFFT2D_getMinY,             METH_VARARGS,"returns min grid value in y-direction"},
		{ "getStepX",            PoissonSolverFFT2D_getStepX,            METH_VARARGS,"returns grid step in x-direction"},
		{ "getStepY",            PoissonSolverFFT2D_getStepY,            METH_VARARGS,"returns grid step in y-direction"},
		{ "findPotential",       PoissonSolverFFT2D_findPotential,       METH_VARARGS,"findPotential(Grid2D rhoGrid2D,Grid2D phiGrid2D)"},
    {NULL}
  };

	// defenition of the memebers of the python PoissonSolverFFT2D wrapper class
	// they will be vailable from python level
	static PyMemberDef PoissonSolverFFT2DClassMembers [] = {
		{NULL}
	};

	//new python PoissonSolverFFT2D wrapper type definition
	static PyTypeObject pyORBIT_PoissonSolverFFT2D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"PoissonSolverFFT2D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) PoissonSolverFFT2D_del , /*tp_dealloc*/
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
		"The PoissonSolverFFT2D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		PoissonSolverFFT2DClassMethods, /* tp_methods */
		PoissonSolverFFT2DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) PoissonSolverFFT2D_init, /* tp_init */
		0, /* tp_alloc */
		PoissonSolverFFT2D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPoissonSolverFFT2D class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initPoissonSolverFFT2D(PyObject* module){
		if (PyType_Ready(&pyORBIT_PoissonSolverFFT2D_Type) < 0) return;
		Py_INCREF(&pyORBIT_PoissonSolverFFT2D_Type);
		PyModule_AddObject(module, "PoissonSolverFFT2D", (PyObject *)&pyORBIT_PoissonSolverFFT2D_Type);
		//std::cout<<"debug PoissonSolverFFT2D added! "<<std::endl;
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
