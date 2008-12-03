#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_boundary2D.hh"
#include "wrap_spacecharge.hh"

#include <iostream>

#include "Boundary2D.hh"
#include "Grid2D.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python Boundary2D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping Boundary2D instance
	//It never will be called directly
	static PyObject* Boundary2D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The Boundary2D new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  Boundary2D class
  //this is implementation of the __init__ method
  static int Boundary2D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		int xBins, yBins;
		double xSize, ySize;
		int BPoints;
		int BModes;
		const char* BShape_name = NULL;
		if(!PyArg_ParseTuple(args,"iiddisi:__init__",&xBins,&yBins,&xSize,&ySize,&BPoints,&BShape_name,&BModes)){
			ORBIT_MPI_Finalize("PyBoundary2D - Boundary2D(nBinX,nBinY,xSize,ySize, nBpoints, BShape, BModes) - constructor needs parameters.");
		}	
		string BShape(BShape_name);
		self->cpp_obj = new Boundary2D(xBins,yBins,xSize,ySize, BPoints, BShape, BModes);
		((Boundary2D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The Boundary2D __init__ has been called!"<<std::endl;
		return 0;
  }
	
	//calculatePhi(double x, double y)
  static PyObject* Boundary2D_calculatePhi(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		double x,y;
		if(!PyArg_ParseTuple(args,"dd:calculatePhi",&x,&y)){
			ORBIT_MPI_Finalize("PyBoundary2D - calculatePhi(x,y) - parameters are needed.");
		}
		return Py_BuildValue("d",cpp_Boundary2D->calculatePhi(x,y));
	}
	
	//get bins in X 
	static PyObject* Boundary2D_getBinsX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("i",cpp_Boundary2D->getBinsX());
	}
	
	//get bins in Y 
	static PyObject* Boundary2D_getBinsY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("i",cpp_Boundary2D->getBinsY());
	}
	
	// getGridMaxX()
	static PyObject* Boundary2D_getGridMaxX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getGridMaxX());
	}
	
	// getGridMaxY()
	static PyObject* Boundary2D_getGridMaxY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getGridMaxY());
	}
	
	// getGridMinX()
	static PyObject* Boundary2D_getGridMinX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getGridMinX());
	}
	
	// getGridMinY()
	static PyObject* Boundary2D_getGridMinY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getGridMinY());
	}	
	
	// getStepX()
	static PyObject* Boundary2D_getStepX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getStepX());
	}	
	
	// getStepY()
	static PyObject* Boundary2D_getStepY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getStepY());
	}			
	
	// getBoundaryShape()
	static PyObject* Boundary2D_getBoundaryShape(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		if(cpp_Boundary2D->getBoundaryShape() == 1){
			return Py_BuildValue("s","Circle");
		}
		if(cpp_Boundary2D->getBoundaryShape() == 2){
			return Py_BuildValue("s","Ellipse");
		}
		if(cpp_Boundary2D->getBoundaryShape() == 3){
			return Py_BuildValue("s","Rectangle");
		}
		Py_INCREF(Py_None);
    return Py_None;	
	}		
	
	// getBoundarySizeX()
	static PyObject* Boundary2D_getBoundarySizeX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getBoundarySizeX());
	}	
	
	// getBoundarySizeY()
	static PyObject* Boundary2D_getBoundarySizeY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getBoundarySizeY());
	}		
	
	// getNumbBPoints()
	static PyObject* Boundary2D_getNumbBPoints(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("i",cpp_Boundary2D->getNumbBPoints());
	}		
	
	//getBoundPointX(int i)
	static PyObject* Boundary2D_getBoundPointX(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		int i;
		if(!PyArg_ParseTuple(args,"i:getBoundPointX",&i)){
			ORBIT_MPI_Finalize("PyBoundary2D - getBoundPointX(int i) - parameter is needed.");
		}		
		return Py_BuildValue("d",cpp_Boundary2D->getBoundPointX(i));
	}		
	
	//getBoundPointY(int i)
	static PyObject* Boundary2D_getBoundPointY(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		int i;
		if(!PyArg_ParseTuple(args,"i:getBoundPointY",&i)){
			ORBIT_MPI_Finalize("PyBoundary2D - getBoundPointY(int i) - parameter is needed.");
		}		
		return Py_BuildValue("d",cpp_Boundary2D->getBoundPointY(i));
	}		
		
	//findPotential(Grid2D* rhoGrid2D,Grid2D* phiGrid2D)
  static PyObject* Boundary2D_findPotential(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		PyObject* pyRhoG;
		PyObject* pyPhiG;
		if(!PyArg_ParseTuple(args,"OO:__init__",&pyRhoG,&pyPhiG)){		
			ORBIT_MPI_Finalize("PyBoundary2D.findPotential(Grid2D rhoGrid2D,Grid2D phiGrid2D) - method needs parameters.");
		}	
		PyObject* pyORBIT_Grid2D_Type = getSpaceChargeType("Grid2D");
		if(!PyObject_IsInstance(pyRhoG,pyORBIT_Grid2D_Type) || !PyObject_IsInstance(pyPhiG,pyORBIT_Grid2D_Type)){
			ORBIT_MPI_Finalize("PyBoundary2D.findPotential(Grid2D rhoGrid2D,Grid2D phiGrid2D) - method needs parameters.");
		}			
		Grid2D* grid2D_rho = (Grid2D*)(((pyORBIT_Object*) pyRhoG)->cpp_obj);
		Grid2D* grid2D_phi = (Grid2D*)(((pyORBIT_Object*) pyPhiG)->cpp_obj);
		if(cpp_Boundary2D->getBinsX() != grid2D_rho->getBinsX() || 
			 cpp_Boundary2D->getBinsY() != grid2D_rho->getBinsY() ||
			 cpp_Boundary2D->getBinsX() != grid2D_phi->getBinsX() || 
			 cpp_Boundary2D->getBinsY() != grid2D_phi->getBinsY()){
		    ORBIT_MPI_Finalize("PyBoundary2D.findPotential(Grid2D rhoGrid2D,Grid2D phiGrid2D) - grid sizes are wrong.");
			 }
		cpp_Boundary2D->findPotential(grid2D_rho->getArr(),grid2D_phi->getArr());
		Py_INCREF(Py_None);
    return Py_None;	
	}		
	
	//addBoundaryPotential(double** phisc)
  static PyObject* Boundary2D_addBoundaryPotential(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) pyBoundary2D->cpp_obj;
		PyObject* pyPhiG;
		if(!PyArg_ParseTuple(args,"O:__init__",&pyPhiG)){		
			ORBIT_MPI_Finalize("PyBoundary2D.addBoundaryPotential(Grid2D phiGrid2D) - method needs parameter.");
		}	
		PyObject* pyORBIT_Grid2D_Type = getSpaceChargeType("Grid2D");
		if(!PyObject_IsInstance(pyPhiG,pyORBIT_Grid2D_Type)){
			ORBIT_MPI_Finalize("PyBoundary2D.addBoundaryPotential(Grid2D phiGrid2D) - method needs parameter.");
		}			
		Grid2D* grid2D_phi = (Grid2D*)(((pyORBIT_Object*) pyPhiG)->cpp_obj);
		if(cpp_Boundary2D->getBinsX() != grid2D_phi->getBinsX() || 
			 cpp_Boundary2D->getBinsY() != grid2D_phi->getBinsY()){
		    ORBIT_MPI_Finalize("PyBoundary2D.addBoundaryPotential(Grid2D phiGrid2D) - grid sizes are wrong.");
			 }
		cpp_Boundary2D->addBoundaryPotential(grid2D_phi->getArr());
		Py_INCREF(Py_None);
    return Py_None;	
	}	
	
  //-----------------------------------------------------
  //destructor for python Boundary2D class (__del__ method).
  //-----------------------------------------------------
  static void Boundary2D_del(pyORBIT_Object* self){
	  //std::cerr<<"The Boundary2D __del__ has been called! 0"<<std::endl;
		Boundary2D* cpp_Boundary2D = (Boundary2D*) self->cpp_obj;
		delete (cpp_Boundary2D);
		self->ob_type->tp_free((PyObject*)self);
		//std::cerr<<"The Boundary2D __del__ has been called! 1"<<std::endl;
  }
	
	// defenition of the methods of the python Boundary2D wrapper class
	// they will be vailable from python level
  static PyMethodDef Boundary2DClassMethods[] = {
		{ "calculatePhi",        Boundary2D_calculatePhi,        METH_VARARGS,"returns potential from boundary for (x,y) point "},
		{ "getBinsX",            Boundary2D_getBinsX,            METH_VARARGS,"returns the number of grid points in x-direction"},
		{ "getBinsY",            Boundary2D_getBinsY,            METH_VARARGS,"returns the number of grid points in y-direction"},
		{ "getGridMaxX",         Boundary2D_getGridMaxX,         METH_VARARGS,"returns max grid value in x-direction"},
		{ "getGridMaxY",         Boundary2D_getGridMaxY,         METH_VARARGS,"returns max grid value in y-direction"},
		{ "getGridMinX",         Boundary2D_getGridMinX,         METH_VARARGS,"returns min grid value in x-direction"},
		{ "getGridMinY",         Boundary2D_getGridMinY,         METH_VARARGS,"returns min grid value in y-direction"},
		{ "getStepX",            Boundary2D_getStepX,            METH_VARARGS,"returns grid step in x-direction"},
		{ "getStepY",            Boundary2D_getStepY,            METH_VARARGS,"returns grid step in y-direction"},
		{ "getBoundaryShape",    Boundary2D_getBoundaryShape,    METH_VARARGS,"returns boundary shape Circle,Ellipse, or Rectangle"},
		{ "getBoundarySizeX",    Boundary2D_getBoundarySizeX,    METH_VARARGS,"returns size (from min to max) in x-direction"},
		{ "getBoundarySizeY",    Boundary2D_getBoundarySizeY,    METH_VARARGS,"returns size (from min to max) in y-direction"},
		{ "getNumbBPoints",      Boundary2D_getNumbBPoints,      METH_VARARGS,"returns number of points on boundary"},
		{ "getBoundPointX",      Boundary2D_getBoundPointX,      METH_VARARGS,"returns x-coord. of the i-th boundary point"},
		{ "getBoundPointY",      Boundary2D_getBoundPointY,      METH_VARARGS,"returns y-coord. of the i-th boundary point"},
		{ "findPotential",       Boundary2D_findPotential,       METH_VARARGS,"findPotential(Grid2D rhoGrid2D,Grid2D phiGrid2D)"},
		{ "addBoundaryPotential",Boundary2D_addBoundaryPotential,METH_VARARGS,"addBoundaryPotential(Grid2D phiGrid2D)"},

    {NULL}
  };

	// defenition of the memebers of the python Boundary2D wrapper class
	// they will be vailable from python level
	static PyMemberDef Boundary2DClassMembers [] = {
		{NULL}
	};

	//new python Boundary2D wrapper type definition
	static PyTypeObject pyORBIT_Boundary2D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Boundary2D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Boundary2D_del , /*tp_dealloc*/
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
		"The Boundary2D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		Boundary2DClassMethods, /* tp_methods */
		Boundary2DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Boundary2D_init, /* tp_init */
		0, /* tp_alloc */
		Boundary2D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyBoundary2D class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initBoundary2D(PyObject* module){
		if (PyType_Ready(&pyORBIT_Boundary2D_Type) < 0) return;
		Py_INCREF(&pyORBIT_Boundary2D_Type);
		PyModule_AddObject(module, "Boundary2D", (PyObject *)&pyORBIT_Boundary2D_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
