#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_boundary2d.hh"
#include "wrap_spacecharge.hh"

#include <iostream>

#include "BaseBoundary2D.hh"
#include "ShapedBoundary2D.hh"
#include "Grid2D.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python BaseBoundary2D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping BaseBoundary2D instance
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
		int nPoints, nModes;
    //if nVars == 2 this is BaseBoundary2D(nPoints, nModes)
    //if nVars == 4 or 5 this is ShapedBoundary2D(nPoints,nModes,shape,xDim,yDim)
    int nVars = PyTuple_Size(args);		
		if(nVars == 2){
			if(!PyArg_ParseTuple(args,"ii:__init__",&nPoints,&nModes)){
				ORBIT_MPI_Finalize("PyBoundary2D - Boundary2D(nPoints,nModes) - constructor needs parameters.");
			}	
			self->cpp_obj = new BaseBoundary2D(nPoints, nModes);
			((BaseBoundary2D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
			//std::cerr<<"The Boundary2D __init__ has been called!"<<std::endl;
			return 0;
		}
		if(nVars == 4 || nVars == 5){
			const char* shape_name = NULL;
			double xDim,yDim;
			if(!PyArg_ParseTuple(args,"iisd|d:__init__",&nPoints,&nModes,&shape_name,&xDim,&yDim)){
				ORBIT_MPI_Finalize("PyBoundary2D - Boundary2D(nPoints,nModes,shape,xDim,yDim) - constructor needs parameters.");
			}			
			string shape(shape_name);	
			if(nVars == 4){ yDim = xDim;}
			self->cpp_obj = new ShapedBoundary2D(nPoints, nModes,shape,xDim,yDim);
			((BaseBoundary2D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
			//std::cerr<<"The Boundary2D __init__ has been called!"<<std::endl;
			return 0;			
		}
		ORBIT_MPI_Finalize("PyBoundary2D - Boundary2D(nPoints,nModes|shape,xDim,yDim) - constructor needs parameters.");
		return 0;		
  }
	
	//getMaxX()
  static PyObject* Boundary2D_getMaxX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getMaxX());
	}
	
	//getMinX()
  static PyObject* Boundary2D_getMinX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getMinX());
	}
	//getMaxY()
  static PyObject* Boundary2D_getMaxY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getMaxY());
	}
	
	//getMinY()
  static PyObject* Boundary2D_getMinY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_Boundary2D->getMinY());
	}
	
	//getNumberOfPoints()
	static PyObject* Boundary2D_getNumberOfPoints(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("i",cpp_Boundary2D->getNumberOfPoints());
	}
	
	//getNumberOfModes()
	static PyObject* Boundary2D_getNumberOfModes(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("i",cpp_Boundary2D->getNumberOfModes());
	}
	
	//getShapeName()
	static PyObject* Boundary2D_getShapeName(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("s",cpp_Boundary2D->getShapeName().c_str());
	}
	
	
	//isInitialized()
	static PyObject* Boundary2D_isInitialized(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		return Py_BuildValue("i",cpp_Boundary2D->isInitialized());
	}
	
	//getBoundaryPoint(index) - returns tuple (x,y)
	static PyObject* Boundary2D_getBoundaryPoint(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		int index = -1;
		if(!PyArg_ParseTuple(args,"i:getBoundaryPoint",&index)){
			ORBIT_MPI_Finalize("PyBoundary2D - getBoundaryPoint(index) - constructor needs a parameter.");
		}			
		double x = cpp_Boundary2D->getBoundaryPointX(index);
		double y = cpp_Boundary2D->getBoundaryPointY(index);
		return Py_BuildValue("(dd)",x,y);
	}
	
	//setBoundaryPoint(int index, double x, double y)
	static PyObject* Boundary2D_setBoundaryPoint(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		double x,y;
		int index = -1;
		if(!PyArg_ParseTuple(args,"idd:setBoundaryPoint",&index,&x,&y)){
			ORBIT_MPI_Finalize("PyBoundary2D - setBoundaryPoint(index,x,y) - method needs parameters.");
		}			
		cpp_Boundary2D->setBoundaryPoint(index,x,y);
		Py_INCREF(Py_None);
    return Py_None;
	}
	
	//isInside(double x, double y) - returns -1 or +1 for "no" and "yes"
	static PyObject* Boundary2D_isInside(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		double x,y;
		if(!PyArg_ParseTuple(args,"dd:isInside",&x,&y)){
			ORBIT_MPI_Finalize("PyBoundary2D - isInside(x,y) - method needs parameters.");
		}			
		return Py_BuildValue("i",cpp_Boundary2D->isInside(x,y));
	}	
	
	// initialize() - initializes the all arrays
	static PyObject* Boundary2D_initialize(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		cpp_Boundary2D->initializeBPs();
		Py_INCREF(Py_None);
    return Py_None;
	}	
	
	// addBoundaryPotential(Grid2D* rhoGrid,Grid2D*  phiGrid);
  static PyObject* Boundary2D_addBoundaryPotential(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) pyBoundary2D->cpp_obj;
		PyObject* pyRhoG;
		PyObject* pyPhiG;
		if(!PyArg_ParseTuple(args,"OO:__init__",&pyRhoG,&pyPhiG)){		
			ORBIT_MPI_Finalize("PyBoundary2D.addBoundaryPotential(rhoGrid2D,phiGrid2D) - method needs parameters.");
		}	
		PyObject* pyORBIT_Grid2D_Type = getSpaceChargeType("Grid2D");
		if(!PyObject_IsInstance(pyRhoG,pyORBIT_Grid2D_Type) || !PyObject_IsInstance(pyPhiG,pyORBIT_Grid2D_Type)){
			ORBIT_MPI_Finalize("PyBoundary2D.addBoundaryPotential(rhoGrid2D,phiGrid2D) - method needs parameters.");
		}			
		Grid2D* grid2D_rho = (Grid2D*)(((pyORBIT_Object*) pyRhoG)->cpp_obj);
		Grid2D* grid2D_phi = (Grid2D*)(((pyORBIT_Object*) pyPhiG)->cpp_obj);
		cpp_Boundary2D->addBoundaryPotential(grid2D_rho,grid2D_phi);
		Py_INCREF(Py_None);
    return Py_None;	
	}		
		
  //-----------------------------------------------------
  //destructor for python Boundary2D class (__del__ method).
  //-----------------------------------------------------
  static void Boundary2D_del(pyORBIT_Object* self){
	  //std::cerr<<"The Boundary2D __del__ has been called! 0"<<std::endl;
		BaseBoundary2D* cpp_Boundary2D = (BaseBoundary2D*) self->cpp_obj;
		if(cpp_Boundary2D != NULL){
			delete cpp_Boundary2D;
		}
		self->ob_type->tp_free((PyObject*)self);
		//std::cerr<<"The Boundary2D __del__ has been called! 1"<<std::endl;
  }
	
	// defenition of the methods of the python Boundary2D wrapper class
	// they will be vailable from python level
  static PyMethodDef Boundary2DClassMethods[] = {
		{ "getMaxX",             Boundary2D_getMaxX,             METH_VARARGS,"returns max value in x-direction"},
		{ "getMinX",             Boundary2D_getMinX,             METH_VARARGS,"returns min value in x-direction"},
		{ "getMaxY",             Boundary2D_getMaxY,             METH_VARARGS,"returns max value in y-direction"},
		{ "getMinY",             Boundary2D_getMinY,             METH_VARARGS,"returns min value in y-direction"},
		{ "getNumberOfModes",    Boundary2D_getNumberOfModes,    METH_VARARGS,"returns the number of free-space modes"},
		{ "getNumberOfPoints",   Boundary2D_getNumberOfPoints,   METH_VARARGS,"returns the number of boundary points"},
		{ "getShapeName",        Boundary2D_getShapeName,        METH_VARARGS,"returns the shape name"},
		{ "isInside",            Boundary2D_isInside,            METH_VARARGS,"returns -1 or +1 if (x,y) point is inside the boundary"},
		{ "isInitialized",       Boundary2D_isInitialized,       METH_VARARGS,"returns 1 if initialized and 0 if not"},
		{ "initialize",          Boundary2D_initialize,          METH_VARARGS,"initializes all boundary arrays"},
		{ "getBoundaryPoint",    Boundary2D_getBoundaryPoint,    METH_VARARGS,"returns tuple (x,y) - a boundary point with index"},
		{ "setBoundaryPoint",    Boundary2D_setBoundaryPoint,    METH_VARARGS,"sets a boundary point (ind,x,y)"},
		{ "addBoundaryPotential",Boundary2D_addBoundaryPotential,METH_VARARGS,"addBoundaryPotential(Grid2D rhoGrid2D,Grid2D phiGrid2D)"},
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
	//It will be called from SpaceCharge wrapper initialization
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
