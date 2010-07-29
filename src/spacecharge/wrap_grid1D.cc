#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_grid1D.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"
#include "wrap_mpi_comm.hh"

#include <iostream>

#include "Grid1D.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python Grid1D class definition
	//---------------------------------------------------------
	
	//constructor for python class wrapping Grid1D instance
	//It never will be called directly
	static PyObject* Grid1D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The Grid1D new has been called!"<<std::endl;
		return (PyObject *) self;
	}	
	
  //initializator for python  Grid1D class
  //this is implementation of the __init__ method
  static int Grid1D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
  	  int binZ;
  	  double zMin = -1.0, zMax = +1.0;
  	  if(!PyArg_ParseTuple(args,"i|dd:__init__",&binZ,&zMin,&zMax)){
				ORBIT_MPI_Finalize("PyGrid1D - Grid1D(0Z[,zMin,zMax]) - constructor needs parameters.");
		}		
		self->cpp_obj = new Grid1D(binZ,zMin,zMax);	
		((Grid1D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The Grid1D __init__ has been called!"<<std::endl;
		return 0;
  }
  //setZero()
  static PyObject* Grid1D_setZero(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;
		cpp_Grid1D->setZero();
		Py_INCREF(Py_None);
    return Py_None;	
	}
	
  //getValue(double z)
  static PyObject* Grid1D_getValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;
		double z;
		if(!PyArg_ParseTuple(args,"d:getValue",&z)){
			ORBIT_MPI_Finalize("PyGrid1D - getValue(z) - parameters are needed.");
		}
		return Py_BuildValue("d",cpp_Grid1D->getValue(z));
	}
	
  //setValue(double value, int iz)
  static PyObject* Grid1D_setValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;
		double val;
		int iz;
		if(!PyArg_ParseTuple(args,"di:setValue",&val,&iz)){
			ORBIT_MPI_Finalize("PyGrid1D - setValue(val,iz) - parameters are needed.");
		}
		cpp_Grid1D->setValue(val,iz);
		Py_INCREF(Py_None);
		return Py_None;	
	}
  //binBunch(Bunch* bunch)
  static PyObject* Grid1D_binBunch(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;
		PyObject* pyBunch;
		if(!PyArg_ParseTuple(args,"O:binBunch",&pyBunch)){
			ORBIT_MPI_Finalize("PyGrid1D - binBunch(Bunch* bunch) - parameter are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PyGrid1D - binBunch(Bunch* bunch) - constructor needs a Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		cpp_Grid1D->binBunch(cpp_bunch);
		Py_INCREF(Py_None);
		return Py_None;	
	}	

  //binValue(double value, double z)
  static PyObject* Grid1D_binValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;
		double val,z;
		if(!PyArg_ParseTuple(args,"dd:binValue",&val,&z)){
			ORBIT_MPI_Finalize("PyGrid1D - binValue(val,z) - parameters are needed.");
		}
		cpp_Grid1D->binValue(val,z);
		Py_INCREF(Py_None);
		return Py_None;	
	}	
	
  //calcGradient(double z)
  static PyObject* Grid1D_calcGradient(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;
		double z;
		double ez;
		if(!PyArg_ParseTuple(args,"d:calcGradient",&z)){
			ORBIT_MPI_Finalize("PyGrid1D - calcGradient(z) - parameters are needed.");
		}
		cpp_Grid1D->calcGradient(z,ez);
		return Py_BuildValue("d",ez);
	}	
		
  //getMinZ()
  static PyObject* Grid1D_getMinZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;	
		return Py_BuildValue("d",cpp_Grid1D->getMinZ());
	}	
	
  //getMaxZ()
  static PyObject* Grid1D_getMaxZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;	
		return Py_BuildValue("d",cpp_Grid1D->getMaxZ());
	}
	
  //getStepZ()
  static PyObject* Grid1D_getStepZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;	
		return Py_BuildValue("d",cpp_Grid1D->getStepZ());
	}	

  //setGridZ()	
  static PyObject* Grid1D_setGridZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;
		double min,max;
		if(!PyArg_ParseTuple(args,"dd:setGridZ",&min,&max)){
			ORBIT_MPI_Finalize("PyGrid1D - setGridZ(min,max) - parameters are needed.");
		}
		cpp_Grid1D->setGridZ(min,max);
		Py_INCREF(Py_None);
		return Py_None;	
	}
	
  //getGridZ(iz)
  static PyObject* Grid1D_getGridZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(args,"i:getGridZ",&ind) || ind < 0 || ind >= cpp_Grid1D->getSizeZ()){
			ORBIT_MPI_Finalize("PyGrid1D - getGridZ(iz) - parameter is needed. [0 - sizeZ[");
		}		
		return Py_BuildValue("d",cpp_Grid1D->getGridZ(ind));
	}	
	
  //getSizeZ()
  static PyObject* Grid1D_getSizeZ(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;	
		return Py_BuildValue("i",cpp_Grid1D->getSizeZ());
	}
	
  //It will synchronize through the MPI communicator
  static PyObject* Grid1D_synchronizeMPI(PyObject *self, PyObject *args){
     pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;
		int nVars = PyTuple_Size(args);
		if(nVars == 0){
			cpp_Grid1D->synchronizeMPI(NULL);
		}
		else {
			PyObject* py_mpi_comm_type = wrap_orbit_mpi_comm::getMPI_CommType("MPI_Comm");
			PyObject* pyMPIComm = PyTuple_GetItem(args,0);			
			if((!PyObject_IsInstance(pyMPIComm,py_mpi_comm_type))){
				ORBIT_MPI_Finalize("Grid1D.synchronizeMPI(MPI_Comm) - input parameter is not MPI_Comm");
			}					
			cpp_Grid1D->synchronizeMPI((pyORBIT_MPI_Comm*) pyMPIComm);
		}
	 	Py_INCREF(Py_None);
		return Py_None; 
  }
  
  //getValueOnGrid(int index)
  static PyObject* Grid1D_getValueOnGrid(PyObject *self, PyObject *args){
     pyORBIT_Object* pyGrid1D = (pyORBIT_Object*) self;
		Grid1D* cpp_Grid1D = (Grid1D*) pyGrid1D->cpp_obj;
		int index;
		if(!PyArg_ParseTuple(args,"i:getValueOnGrid",&index)){
			ORBIT_MPI_Finalize("PyGrid1D - getValueOnGrid(index) - parameters are needed.");
		}
		return Py_BuildValue("d",cpp_Grid1D->getValueOnGrid(index));
  }
  //-----------------------------------------------------
  //destructor for python Grid1D class (__del__ method).
  //-----------------------------------------------------
  static void Grid1D_del(pyORBIT_Object* self){
		//std::cerr<<"The Grid1D __del__ has been called!"<<std::endl;
		Grid1D* cpp_Grid1D = (Grid1D*) self->cpp_obj;
		delete cpp_Grid1D;
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python Grid1D wrapper class
	// they will be vailable from python level
  static PyMethodDef Grid1DClassMethods[] = {
		{ "setZero",      Grid1D_setZero,     METH_VARARGS,"sets all points on the grid to zero"},
		{ "getValue",     Grid1D_getValue,    METH_VARARGS,"returns value for (z) point"},
		{ "setValue",     Grid1D_setValue,    METH_VARARGS,"sets value for (iz) point"},
		{ "setGridZ",     Grid1D_setGridZ,    METH_VARARGS,"sets the Z grid with min,max"},
		{ "getGridZ",     Grid1D_getGridZ,    METH_VARARGS,"returns the z-grid point with index ind"},
		{ "getSizeZ",     Grid1D_getSizeZ,    METH_VARARGS,"returns the size of grid in Z dir."},
		{ "getMinZ",      Grid1D_getMinZ,     METH_VARARGS,"returns the min grid point in Z dir."},
		{ "getMaxZ",      Grid1D_getMaxZ,     METH_VARARGS,"returns the max grid point in Z dir."},
		{ "getStepZ",     Grid1D_getStepZ,     METH_VARARGS,"returns the step in Z dir."},
		{ "binValue",     Grid1D_binValue,    METH_VARARGS,"bins the value into the 1D mesh"},
		{ "binBunch",     Grid1D_binBunch,    METH_VARARGS,"bins the Bunch instance into the 1D mesh"},
		{ "calcGradient", Grid1D_calcGradient,METH_VARARGS,"returns gradient as (gz) for point (z)"},
		{ "synchronizeMPI",Grid1D_synchronizeMPI,METH_VARARGS,"synchronize through the MPI communicator"},
    {NULL}
  };

	// defenition of the memebers of the python Grid1D wrapper class
	// they will be vailable from python level
	static PyMemberDef Grid1DClassMembers [] = {
		{NULL}
	};

	//new python Grid1D wrapper type definition
	static PyTypeObject pyORBIT_Grid1D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Grid1D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Grid1D_del , /*tp_dealloc*/
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
		"The Grid1D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		Grid1DClassMethods, /* tp_methods */
		Grid1DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Grid1D_init, /* tp_init */
		0, /* tp_alloc */
		Grid1D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyGrid1D class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initGrid1D(PyObject* module){
		if (PyType_Ready(&pyORBIT_Grid1D_Type) < 0) return;
		Py_INCREF(&pyORBIT_Grid1D_Type);
		PyModule_AddObject(module, "Grid1D", (PyObject *)&pyORBIT_Grid1D_Type);
		//std::cout<<"debug Grid1D added! "<<std::endl;
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
