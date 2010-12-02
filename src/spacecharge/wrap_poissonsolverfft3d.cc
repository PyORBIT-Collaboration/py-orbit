#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "PoissonSolverFFT3D.hh"
#include "Grid3D.hh"

#include "wrap_poissonsolverfft2d.hh"
#include "wrap_spacecharge.hh"

#include <iostream>

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python PoissonSolverFFT3D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping PoissonSolverFFT3D instance
	//It never will be called directly
	static PyObject* PoissonSolverFFT3D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The PoissonSolverFFT3D new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  PoissonSolverFFT3D class
  //this is implementation of the __init__ method
  static int PoissonSolverFFT3D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		int xSize, ySize, zSize;
		double xMin = -1.0, xMax = +1.0;
		double yMin = -1.0, yMax = +1.0;
		double zMin = -1.0, zMax = +1.0;
		if(!PyArg_ParseTuple(args,"iii|dddddd:__init__",&xSize,&ySize,&zSize,&xMin,&xMax,&yMin,&yMax,&zMin,&zMax)){
			ORBIT_MPI_Finalize("PyPoissonSolverFFT3D - PoissonSolverFFT3D(nX,nY,nZ,[,xMin,xMax,yMin,yMax,zMin,zMax]) - constructor needs parameters.");
		}	
		self->cpp_obj = new PoissonSolverFFT3D(xSize,ySize,zSize,xMin,xMax,yMin,yMax,zMin,zMax);
		((PoissonSolverFFT3D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The PoissonSolverFFT3D __init__ has been called!"<<std::endl;
		return 0;
  }
		
	//get grid size in X 
	static PyObject* PoissonSolverFFT3D_getSizeX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("i",cpp_PoissonSolverFFT3D->getSizeX());
	}
	
	//get grid size in Y 
	static PyObject* PoissonSolverFFT3D_getSizeY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("i",cpp_PoissonSolverFFT3D->getSizeY());
	}

	//get grid size in Z
	static PyObject* PoissonSolverFFT3D_getSizeZ(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("i",cpp_PoissonSolverFFT3D->getSizeZ());
	}	

	// getMaxX()
	static PyObject* PoissonSolverFFT3D_getMaxX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT3D->getMaxX());
	}
	
	// getMaxY()
	static PyObject* PoissonSolverFFT3D_getMaxY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT3D->getMaxY());
	}
	
	// getMaxZ()
	static PyObject* PoissonSolverFFT3D_getMaxZ(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT3D->getMaxZ());
	}
	
	// getMinX()
	static PyObject* PoissonSolverFFT3D_getMinX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT3D->getMinX());
	}
	
	// getMinY()
	static PyObject* PoissonSolverFFT3D_getMinY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT3D->getMinY());
	}	
	
	// getMinZ()
	static PyObject* PoissonSolverFFT3D_getMinZ(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT3D->getMinZ());
	}	

	// getStepX()
	static PyObject* PoissonSolverFFT3D_getStepX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT3D->getStepX());
	}	
	
	// getStepY()
	static PyObject* PoissonSolverFFT3D_getStepY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT3D->getStepY());
	}		
	
	// getStepZ()
	static PyObject* PoissonSolverFFT3D_getStepZ(PyObject *self, PyObject *args){
		pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		return Py_BuildValue("d",cpp_PoissonSolverFFT3D->getStepZ());
	}			
		
	//findPotential(Grid3D* rhoGrid3D,Grid3D* phiGrid3D)
  static PyObject* PoissonSolverFFT3D_findPotential(PyObject *self, PyObject *args){
    pyORBIT_Object* pyPoissonSolverFFT3D = (pyORBIT_Object*) self;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) pyPoissonSolverFFT3D->cpp_obj;
		PyObject* pyRhoG;
		PyObject* pyPhiG;
		if(!PyArg_ParseTuple(args,"OO:__init__",&pyRhoG,&pyPhiG)){		
			ORBIT_MPI_Finalize("PyPoissonSolverFFT3D.findPotential(Grid3D rhoGrid3D,Grid3D phiGrid3D) - method needs parameters.");
		}	
		PyObject* pyORBIT_Grid3D_Type = getSpaceChargeType("Grid3D");
		if(!PyObject_IsInstance(pyRhoG,pyORBIT_Grid3D_Type) || !PyObject_IsInstance(pyPhiG,pyORBIT_Grid3D_Type)){
			ORBIT_MPI_Finalize("PyPoissonSolverFFT3D.findPotential(Grid3D rhoGrid3D,Grid3D phiGrid3D) - method needs parameters.");
		}			
		Grid3D* grid3D_rho = (Grid3D*)(((pyORBIT_Object*) pyRhoG)->cpp_obj);
		Grid3D* grid3D_phi = (Grid3D*)(((pyORBIT_Object*) pyPhiG)->cpp_obj);
		cpp_PoissonSolverFFT3D->findPotential(grid3D_rho,grid3D_phi);
		Py_INCREF(Py_None);
    return Py_None;
	}		
	
  //-----------------------------------------------------
  //destructor for python PoissonSolverFFT3D class (__del__ method).
  //-----------------------------------------------------
  static void PoissonSolverFFT3D_del(pyORBIT_Object* self){
	  //std::cerr<<"The PoissonSolverFFT3D __del__ has been called! 0"<<std::endl;
		PoissonSolverFFT3D* cpp_PoissonSolverFFT3D = (PoissonSolverFFT3D*) self->cpp_obj;
		if(cpp_PoissonSolverFFT3D != NULL){
			delete cpp_PoissonSolverFFT3D;
		}
		self->ob_type->tp_free((PyObject*)self);
		//std::cerr<<"The PoissonSolverFFT3D __del__ has been called! 1"<<std::endl;
  }
	
	// defenition of the methods of the python PoissonSolverFFT3D wrapper class
	// they will be vailable from python level
  static PyMethodDef PoissonSolverFFT3DClassMethods[] = {
		{ "getSizeX",            PoissonSolverFFT3D_getSizeX,            METH_VARARGS,"returns the number of grid points in x-direction"},
		{ "getSizeY",            PoissonSolverFFT3D_getSizeY,            METH_VARARGS,"returns the number of grid points in y-direction"},
		{ "getSizeZ",            PoissonSolverFFT3D_getSizeZ,            METH_VARARGS,"returns the number of grid points in z-direction"},
		{ "getMaxX",             PoissonSolverFFT3D_getMaxX,             METH_VARARGS,"returns max grid value in x-direction"},
		{ "getMaxY",             PoissonSolverFFT3D_getMaxY,             METH_VARARGS,"returns max grid value in y-direction"},
		{ "getMaxZ",             PoissonSolverFFT3D_getMaxZ,             METH_VARARGS,"returns max grid value in z-direction"},
		{ "getMinX",             PoissonSolverFFT3D_getMinX,             METH_VARARGS,"returns min grid value in x-direction"},
		{ "getMinY",             PoissonSolverFFT3D_getMinY,             METH_VARARGS,"returns min grid value in y-direction"},
		{ "getMinZ",             PoissonSolverFFT3D_getMinZ,             METH_VARARGS,"returns min grid value in z-direction"},
		{ "getStepX",            PoissonSolverFFT3D_getStepX,            METH_VARARGS,"returns grid step in x-direction"},
		{ "getStepY",            PoissonSolverFFT3D_getStepY,            METH_VARARGS,"returns grid step in y-direction"},
		{ "getStepZ",            PoissonSolverFFT3D_getStepZ,            METH_VARARGS,"returns grid step in z-direction"},
		{ "findPotential",       PoissonSolverFFT3D_findPotential,       METH_VARARGS,"findPotential(Grid3D rhoGrid3D,Grid3D phiGrid3D)"},
    {NULL}
  };

	// defenition of the memebers of the python PoissonSolverFFT3D wrapper class
	// they will be vailable from python level
	static PyMemberDef PoissonSolverFFT3DClassMembers [] = {
		{NULL}
	};

	//new python PoissonSolverFFT3D wrapper type definition
	static PyTypeObject pyORBIT_PoissonSolverFFT3D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"PoissonSolverFFT3D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) PoissonSolverFFT3D_del , /*tp_dealloc*/
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
		"The PoissonSolverFFT3D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		PoissonSolverFFT3DClassMethods, /* tp_methods */
		PoissonSolverFFT3DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) PoissonSolverFFT3D_init, /* tp_init */
		0, /* tp_alloc */
		PoissonSolverFFT3D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPoissonSolverFFT3D class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initPoissonSolverFFT3D(PyObject* module){
		if (PyType_Ready(&pyORBIT_PoissonSolverFFT3D_Type) < 0) return;
		Py_INCREF(&pyORBIT_PoissonSolverFFT3D_Type);
		PyModule_AddObject(module, "PoissonSolverFFT3D", (PyObject *)&pyORBIT_PoissonSolverFFT3D_Type);
		//std::cout<<"debug PoissonSolverFFT3D added! "<<std::endl;
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
