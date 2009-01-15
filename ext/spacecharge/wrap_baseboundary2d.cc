#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_baseboundary2d.hh"
#include "wrap_spacecharge.hh"

#include <iostream>

#include "BaseBoundary2D.hh"
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
	static PyObject* BaseBoundary2D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The BaseBoundary2D new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  BaseBoundary2D class
  //this is implementation of the __init__ method
  static int BaseBoundary2D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		int nPoints, nModes;
		if(!PyArg_ParseTuple(args,"ii:__init__",&nPoints,&nModes)){
			ORBIT_MPI_Finalize("PyBaseBoundary2D - BaseBoundary2D(nPoints,nModes) - constructor needs parameters.");
		}	
		self->cpp_obj = new BaseBoundary2D(nPoints, nModes);
		((BaseBoundary2D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The BaseBoundary2D __init__ has been called!"<<std::endl;
		return 0;
  }
	
	//getMaxX()
  static PyObject* BaseBoundary2D_getMaxX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_BaseBoundary2D->getMaxX());
	}
	
	//getMinX()
  static PyObject* BaseBoundary2D_getMinX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_BaseBoundary2D->getMinX());
	}
	//getMaxY()
  static PyObject* BaseBoundary2D_getMaxY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_BaseBoundary2D->getMaxY());
	}
	
	//getMinY()
  static PyObject* BaseBoundary2D_getMinY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		return Py_BuildValue("d",cpp_BaseBoundary2D->getMinY());
	}
	
	//getNumberOfPoints()
	static PyObject* BaseBoundary2D_getNumberOfPoints(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		return Py_BuildValue("i",cpp_BaseBoundary2D->getNumberOfPoints());
	}
	
	//getNumberOfModes()
	static PyObject* BaseBoundary2D_getNumberOfModes(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		return Py_BuildValue("i",cpp_BaseBoundary2D->getNumberOfModes());
	}
	
	//isInitialized()
	static PyObject* BaseBoundary2D_isInitialized(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		return Py_BuildValue("i",cpp_BaseBoundary2D->isInitialized());
	}
	
	//getBoundaryPoint(index) - returns tuple (x,y)
	static PyObject* BaseBoundary2D_getBoundaryPoint(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		int index = -1;
		if(!PyArg_ParseTuple(args,"i:getBoundaryPoint",&index)){
			ORBIT_MPI_Finalize("PyBaseBoundary2D - getBoundaryPoint(index) - constructor needs a parameter.");
		}			
		double x = cpp_BaseBoundary2D->getBoundaryPointX(index);
		double y = cpp_BaseBoundary2D->getBoundaryPointY(index);
		return Py_BuildValue("(dd)",x,y);
	}
	
	//setBoundaryPoint(int index, double x, double y)
	static PyObject* BaseBoundary2D_setBoundaryPoint(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		double x,y;
		int index = -1;
		if(!PyArg_ParseTuple(args,"idd:setBoundaryPoint",&index,&x,&y)){
			ORBIT_MPI_Finalize("PyBaseBoundary2D - setBoundaryPoint(index,x,y) - method needs parameters.");
		}			
		cpp_BaseBoundary2D->setBoundaryPoint(index,x,y);
		Py_INCREF(Py_None);
    return Py_None;
	}
	
	// initialize() - initializes the all arrays
	static PyObject* BaseBoundary2D_initialize(PyObject *self, PyObject *args){
		pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		cpp_BaseBoundary2D->initializeBPs();
		Py_INCREF(Py_None);
    return Py_None;
	}	
	
	// addBoundaryPotential(Grid2D* rhoGrid,Grid2D*  phiGrid);
  static PyObject* BaseBoundary2D_addBoundaryPotential(PyObject *self, PyObject *args){
    pyORBIT_Object* pyBaseBoundary2D = (pyORBIT_Object*) self;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) pyBaseBoundary2D->cpp_obj;
		PyObject* pyRhoG;
		PyObject* pyPhiG;
		if(!PyArg_ParseTuple(args,"OO:__init__",&pyRhoG,&pyPhiG)){		
			ORBIT_MPI_Finalize("PyBaseBoundary2D.addBoundaryPotential(Grid2D rhoGrid2D,Grid2D phiGrid2D) - method needs parameters.");
		}	
		PyObject* pyORBIT_Grid2D_Type = getSpaceChargeType("Grid2D");
		if(!PyObject_IsInstance(pyRhoG,pyORBIT_Grid2D_Type) || !PyObject_IsInstance(pyPhiG,pyORBIT_Grid2D_Type)){
			ORBIT_MPI_Finalize("PyBaseBoundary2D.addBoundaryPotential(Grid2D rhoGrid2D,Grid2D phiGrid2D) - method needs parameters.");
		}			
		Grid2D* grid2D_rho = (Grid2D*)(((pyORBIT_Object*) pyRhoG)->cpp_obj);
		Grid2D* grid2D_phi = (Grid2D*)(((pyORBIT_Object*) pyPhiG)->cpp_obj);
		cpp_BaseBoundary2D->addBoundaryPotential(grid2D_rho,grid2D_phi);
		Py_INCREF(Py_None);
    return Py_None;	
	}		
		
  //-----------------------------------------------------
  //destructor for python BaseBoundary2D class (__del__ method).
  //-----------------------------------------------------
  static void BaseBoundary2D_del(pyORBIT_Object* self){
	  //std::cerr<<"The BaseBoundary2D __del__ has been called! 0"<<std::endl;
		BaseBoundary2D* cpp_BaseBoundary2D = (BaseBoundary2D*) self->cpp_obj;
		if(cpp_BaseBoundary2D != NULL){
			delete cpp_BaseBoundary2D;
		}
		self->ob_type->tp_free((PyObject*)self);
		//std::cerr<<"The BaseBoundary2D __del__ has been called! 1"<<std::endl;
  }
	
	// defenition of the methods of the python BaseBoundary2D wrapper class
	// they will be vailable from python level
  static PyMethodDef BaseBoundary2DClassMethods[] = {
		{ "addBoundaryPotential",BaseBoundary2D_addBoundaryPotential,METH_VARARGS,"adds boundary potential to the phiGrid"},
		{ "getMaxX",             BaseBoundary2D_getMaxX,             METH_VARARGS,"returns max value in x-direction"},
		{ "getMinX",             BaseBoundary2D_getMinX,             METH_VARARGS,"returns min value in x-direction"},
		{ "getMaxY",             BaseBoundary2D_getMaxY,             METH_VARARGS,"returns max value in y-direction"},
		{ "getMinY",             BaseBoundary2D_getMinY,             METH_VARARGS,"returns min value in y-direction"},
		{ "getNumberOfModes",    BaseBoundary2D_getNumberOfModes,    METH_VARARGS,"returns the number of free-space modes"},
		{ "getNumberOfPoints",   BaseBoundary2D_getNumberOfPoints,   METH_VARARGS,"returns the number of boundary points"},
		{ "isInitialized",       BaseBoundary2D_isInitialized,       METH_VARARGS,"returns 1 if initialized and 0 if not"},
		{ "initialize",          BaseBoundary2D_initialize,          METH_VARARGS,"initializes all boundary arrays"},
		{ "getBoundaryPoint",    BaseBoundary2D_getBoundaryPoint,    METH_VARARGS,"returns tuple (x,y) - a boundary point with index"},
		{ "setBoundaryPoint",    BaseBoundary2D_setBoundaryPoint,    METH_VARARGS,"sets a boindary point (ind,x,y)"},
		{ "addBoundaryPotential",BaseBoundary2D_addBoundaryPotential,METH_VARARGS,"addBoundaryPotential(Grid2D phiGrid2D)"},
    {NULL}
  };

	// defenition of the memebers of the python BaseBoundary2D wrapper class
	// they will be vailable from python level
	static PyMemberDef BaseBoundary2DClassMembers [] = {
		{NULL}
	};

	//new python BaseBoundary2D wrapper type definition
	static PyTypeObject pyORBIT_BaseBoundary2D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"BaseBoundary2D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) BaseBoundary2D_del , /*tp_dealloc*/
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
		"The BaseBoundary2D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		BaseBoundary2DClassMethods, /* tp_methods */
		BaseBoundary2DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) BaseBoundary2D_init, /* tp_init */
		0, /* tp_alloc */
		BaseBoundary2D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyBaseBoundary2D class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initBaseBoundary2D(PyObject* module){
		if (PyType_Ready(&pyORBIT_BaseBoundary2D_Type) < 0) return;
		Py_INCREF(&pyORBIT_BaseBoundary2D_Type);
		PyModule_AddObject(module, "BaseBoundary2D", (PyObject *)&pyORBIT_BaseBoundary2D_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
