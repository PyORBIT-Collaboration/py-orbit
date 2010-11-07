#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_grid3D.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "Grid3D.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python Grid3D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping Grid3D instance
	//It never will be called directly
	static PyObject* Grid3D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  Grid3D class
  //this is implementation of the __init__ method
  static int Grid3D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
   int binX, binY, binZ;
	 if(!PyArg_ParseTuple(args,"iii:__init__",&binX,&binY,&binZ)){
				ORBIT_MPI_Finalize("PyGrid3D - Grid3D(nX,nY,nZ) - constructor needs parameters.");
		}		
		self->cpp_obj = new Grid3D(binX,binY,binZ);	
		((Grid3D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		return 0;
  }
		
  //setZero()
  static PyObject* Grid3D_setZero(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		cpp_Grid3D->setZero();
		Py_INCREF(Py_None);
		return Py_None;	
	}
	
	//getValue(double x, double y, double z)
  static PyObject* Grid3D_getValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		double x,y,z;
		if(!PyArg_ParseTuple(args,"ddd:getValue",&x,&y,&z)){
			ORBIT_MPI_Finalize("PyGrid3D - getValue(x,y,z) - parameters are needed.");
		}
		return Py_BuildValue("d",cpp_Grid3D->getValue(x,y,z));
	}
	
	//setValue(double value, int ix, int iy)
  static PyObject* Grid3D_setValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		double val;
		int ix,iy,iz;
		if(!PyArg_ParseTuple(args,"diii:setValue",&val,&ix,&iy,&iz)){
			ORBIT_MPI_Finalize("PyGrid3D - setValue(val,ix,iy,iz) - parameters are needed.");
		}
		cpp_Grid3D->setValue(val,ix,iy,iz);
		Py_INCREF(Py_None);
		return Py_None;	
	}

	//getValueOnGrid(int ix, int iy, int iz)
  static PyObject* Grid3D_getValueOnGrid(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		int ix,iy,iz;
		if(!PyArg_ParseTuple(args,"iii:getValueOnGrid",&ix,&iy,&iz)){
			ORBIT_MPI_Finalize("PyGrid3D - getValueOnGrid(ix,iy,iz) - parameters are needed.");
		}
		return Py_BuildValue("d",cpp_Grid3D->getValueOnGrid(ix,iy,iz));
	}
	
	//setGridX(double min, double max)
  static PyObject* Grid3D_setGridX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		double min,max;
		if(!PyArg_ParseTuple(args,"dd:setGridX",&min,&max)){
			ORBIT_MPI_Finalize("PyGrid3D - setGridX(min,max) - parameters are needed.");
		}
		cpp_Grid3D->setGridX(min,max);
		Py_INCREF(Py_None);
		return Py_None;	
	}
	
	//setGridY(double min, double max, int n)
  static PyObject* Grid3D_setGridY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		double min,max;
		if(!PyArg_ParseTuple(args,"dd:setGridY",&min,&max)){
			ORBIT_MPI_Finalize("PyGrid3D - setGridY(min,max) - parameters are needed.");
		}
		cpp_Grid3D->setGridY(min,max);
		Py_INCREF(Py_None);
		return Py_None;	
	}	
	
	//setGridZ(double min, double max)
  static PyObject* Grid3D_setGridZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		double min,max;
		if(!PyArg_ParseTuple(args,"dd:setGridZ",&min,&max)){
			ORBIT_MPI_Finalize("PyGrid3D - setGridZ(min,max) - parameters are needed.");
		}
		cpp_Grid3D->setGridZ(min,max);
		Py_INCREF(Py_None);
		return Py_None;	
	}	
	
	//getGridX(ix)
  static PyObject* Grid3D_getGridX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(args,"i:getGridX",&ind) || ind < 0 || ind >= cpp_Grid3D->getSizeX()){
			ORBIT_MPI_Finalize("PyGrid3D - getGridX(ix) - parameter is needed. [0 - sizeX[");
		}		
		return Py_BuildValue("d",cpp_Grid3D->getGridX(ind));
	}	
	
	//getGridY(iy)
  static PyObject* Grid3D_getGridY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(args,"i:getGridY",&ind) || ind < 0 || ind >= cpp_Grid3D->getSizeY()){
			ORBIT_MPI_Finalize("PyGrid3D - getGridY(iy) - parameter is needed. [0 - sizeY[");
		}		
		return Py_BuildValue("d",cpp_Grid3D->getGridY(ind));
	}	

	//getGridZ(iy)
  static PyObject* Grid3D_getGridZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(args,"i:getGridZ",&ind) || ind < 0 || ind >= cpp_Grid3D->getSizeZ()){
			ORBIT_MPI_Finalize("PyGrid3D - getGridZ(iy) - parameter is needed. [0 - sizeZ[");
		}		
		return Py_BuildValue("d",cpp_Grid3D->getGridZ(ind));
	}	
	
	//getSizeX()
  static PyObject* Grid3D_getSizeX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;	
		return Py_BuildValue("i",cpp_Grid3D->getSizeX());
	}	
	
	//getSizeY()
  static PyObject* Grid3D_getSizeY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;	
		return Py_BuildValue("i",cpp_Grid3D->getSizeY());
	}	
	
	//getSizeZ()
  static PyObject* Grid3D_getSizeZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;	
		return Py_BuildValue("i",cpp_Grid3D->getSizeZ());
	}	

  //It will synchronize through the MPI communicator
  static PyObject* Grid3D_synchronizeMPI(PyObject *self, PyObject *args){
     pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		int nVars = PyTuple_Size(args);
		if(nVars == 0){
			cpp_Grid3D->synchronizeMPI(NULL);
		}
		else {
			PyObject* py_mpi_comm_type = wrap_orbit_mpi_comm::getMPI_CommType("MPI_Comm");
			PyObject* pyMPIComm = PyTuple_GetItem(args,0);			
			if((!PyObject_IsInstance(pyMPIComm,py_mpi_comm_type))){
				ORBIT_MPI_Finalize("Grid3D.synchronizeMPI(MPI_Comm) - input parameter is not MPI_Comm");
			}					
			cpp_Grid3D->synchronizeMPI((pyORBIT_MPI_Comm*) pyMPIComm);
		}
	 	Py_INCREF(Py_None);
		return Py_None; 
  }	
	
	//getMinX()
  static PyObject* Grid3D_getMinX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;	
		return Py_BuildValue("d",cpp_Grid3D->getMinX());
	}	
	
	//getMaxX()
  static PyObject* Grid3D_getMaxX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;	
		return Py_BuildValue("d",cpp_Grid3D->getMaxX());
	}	
	
	//getMinY()
  static PyObject* Grid3D_getMinY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;	
		return Py_BuildValue("d",cpp_Grid3D->getMinY());
	}	
	
	//getMaxY()
  static PyObject* Grid3D_getMaxY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;	
		return Py_BuildValue("d",cpp_Grid3D->getMaxY());
	}		
	
	//getMinZ()
  static PyObject* Grid3D_getMinZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;	
		return Py_BuildValue("d",cpp_Grid3D->getMinZ());
	}	
	
	//getMaxZ()
  static PyObject* Grid3D_getMaxZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;	
		return Py_BuildValue("d",cpp_Grid3D->getMaxZ());
	}		
	
	
	//binBunch(Bunch* bunch)
  static PyObject* Grid3D_binBunch(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		PyObject* pyBunch;
		if(!PyArg_ParseTuple(args,"O:binBunch",&pyBunch)){
			ORBIT_MPI_Finalize("PyGrid3D - binBunch(Bunch* bunch) - parameter are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PyGrid3D - binBunch(Bunch* bunch) - method needs a Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		cpp_Grid3D->binBunch(cpp_bunch);
		Py_INCREF(Py_None);
    return Py_None;	
	}		
	
	//binValue(double value, double x, double y, double z)
  static PyObject* Grid3D_binValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		double val,x,y,z;
		if(!PyArg_ParseTuple(args,"dddd:binValue",&val,&x,&y,&z)){
			ORBIT_MPI_Finalize("PyGrid3D - binValue(val,x,y,z) - parameters are needed.");
		}
		cpp_Grid3D->binValue(val,x,y,z);
		Py_INCREF(Py_None);
    return Py_None;	
	}	
	
	//calcGradient(double x, double y, double y)
  static PyObject* Grid3D_calcGradient(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid3D = (pyORBIT_Object*) self;
		Grid3D* cpp_Grid3D = (Grid3D*) pyGrid3D->cpp_obj;
		double x,y,z;
		double ex, ey, ez;
		if(!PyArg_ParseTuple(args,"ddd:calcGradient",&x,&y,&z)){
			ORBIT_MPI_Finalize("PyGrid3D - calcGradient(x,y,z) - parameters are needed.");
		}
		cpp_Grid3D->calcGradient(x,ex,y,ey,z,ez);
		return Py_BuildValue("(ddd)",ex,ey,ez);
	}
	
  //-----------------------------------------------------
  //destructor for python Grid3D class (__del__ method).
  //-----------------------------------------------------
  static void Grid3D_del(pyORBIT_Object* self){
		//std::cerr<<"The Grid3D __del__ has been called!"<<std::endl;
		Grid3D* cpp_Grid3D = (Grid3D*) self->cpp_obj;
		delete cpp_Grid3D;
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python Grid3D wrapper class
	// they will be vailable from python level
  static PyMethodDef Grid3DClassMethods[] = {
		{ "setZero",       Grid3D_setZero,       METH_VARARGS,"sets all points on the grid to zero"},
		{ "getValue",      Grid3D_getValue,      METH_VARARGS,"returns value for (x,y,z) point"},
		{ "getValueOnGrid",Grid3D_getValueOnGrid,METH_VARARGS,"returns value for indeces (ix,iy,iz) "},
		{ "setValue",      Grid3D_setValue,      METH_VARARGS,"sets value for (ix,iy,iz) point - (val,ix,iy,iz)"},
		{ "setGridX",      Grid3D_setGridX,      METH_VARARGS,"sets the X grid with min,max"},
		{ "setGridY",      Grid3D_setGridY,      METH_VARARGS,"sets the Y grid with min,max"},
		{ "setGridZ",      Grid3D_setGridZ,      METH_VARARGS,"sets the Z grid with min,max"},
		{ "getGridX",      Grid3D_getGridX,      METH_VARARGS,"returns the x-grid point with index ind"},
		{ "getGridY",      Grid3D_getGridY,      METH_VARARGS,"returns the x-grid point with index ind"},
		{ "getGridZ",      Grid3D_getGridZ,      METH_VARARGS,"sets the Z grid with min,max"},
		{ "getSizeX",      Grid3D_getSizeX,      METH_VARARGS,"returns the size of grid in X dir."},
		{ "getSizeY",      Grid3D_getSizeY,      METH_VARARGS,"returns the size of grid in Y dir."},
		{ "getSizeZ",      Grid3D_getSizeZ,      METH_VARARGS,"returns the size of grid in Z dir."},
		{ "getMinX",       Grid3D_getMinX,       METH_VARARGS,"returns the min grid point in X dir."},
		{ "getMaxX",       Grid3D_getMaxX,       METH_VARARGS,"returns the max grid point in X dir."},
		{ "getMinY",       Grid3D_getMinY,       METH_VARARGS,"returns the min grid point in Y dir."},
		{ "getMaxY",       Grid3D_getMaxY,       METH_VARARGS,"returns the max grid point in Y dir."},
		{ "getMinZ",       Grid3D_getMinZ,       METH_VARARGS,"returns the min grid point in Z dir."},
		{ "getMaxZ",       Grid3D_getMaxZ,       METH_VARARGS,"returns the max grid point in Z dir."},
		{ "binValue",      Grid3D_binValue,      METH_VARARGS,"bins the value into the 3D mesh"},
		{ "binBunch",      Grid3D_binBunch,      METH_VARARGS,"bins the Bunch instance into the 3D mesh"},
		{ "calcGradient",  Grid3D_calcGradient,  METH_VARARGS,"returns gradient as (gx,gy,gz) for point (x,y,z)"},
		{ "synchronizeMPI",Grid3D_synchronizeMPI,METH_VARARGS,"synchronize through the MPI communicator"},		
    {NULL}
  };

	// defenition of the memebers of the python Grid3D wrapper class
	// they will be vailable from python level
	static PyMemberDef Grid3DClassMembers [] = {
		{NULL}
	};

	//new python Grid3D wrapper type definition
	static PyTypeObject pyORBIT_Grid3D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Grid3D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Grid3D_del , /*tp_dealloc*/
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
		"The Grid3D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		Grid3DClassMethods, /* tp_methods */
		Grid3DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Grid3D_init, /* tp_init */
		0, /* tp_alloc */
		Grid3D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyGrid3D class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initGrid3D(PyObject* module){
		if (PyType_Ready(&pyORBIT_Grid3D_Type) < 0) return;
		Py_INCREF(&pyORBIT_Grid3D_Type);
		PyModule_AddObject(module, "Grid3D", (PyObject *)&pyORBIT_Grid3D_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
