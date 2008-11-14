#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_grid2D.hh"
#include "wrap_spacecharge.hh"

#include <iostream>

#include "Grid2D.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python Grid2D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping Grid2D instance
	//It never will be called directly
	static PyObject* Grid2D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The Grid2D new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  Grid2D class
  //this is implementation of the __init__ method
  static int Grid2D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		int nVars = PyTuple_Size(args);
		if(nVars == 2){
			int binX, binY;
			if(!PyArg_ParseTuple(args,"ii:__init__",&binX,&binY)){
				ORBIT_MPI_Finalize("PyGrid2D - Grid2D(nBinX,nBinY) - constructor needs parameters.");
			}
			self->cpp_obj = new Grid2D(binX, binY);
		}
		if(nVars == 6){
			int binX, binY;
			double xMin,yMin,xMax,yMax;
			if(!PyArg_ParseTuple(args,"ddiddi:__init__",&xMin,&xMax,&binX,&yMin,&yMax,&binY)){
				ORBIT_MPI_Finalize("PyGrid2D - Grid2D(xMin,xMax,nX,yMin,yMax,nY) - constructor needs parameters.");
			}		
			self->cpp_obj = new Grid2D(xMin,xMax,binX,yMin,yMax,binY);	
		}
		if(nVars == 1){
			PyObject* pyB;
			if(!PyArg_ParseTuple(args,"O:__init__",&pyB)){		
				ORBIT_MPI_Finalize("PyGrid2D - Grid2D(Boundary2D) - constructor needs a parameter.");
			}	
			PyObject* pyORBIT_Boindary2D_Type = getSpaceChargeType("Boundary2D");
			if(!PyObject_IsInstance(pyB,pyORBIT_Boindary2D_Type)){
				ORBIT_MPI_Finalize("PyPhaseVector - Grid2D(Boundary2D) - constructor needs a parameter.");
			}			
			Boundary2D* boundary2D = (Boundary2D*)(((pyORBIT_Object*) pyB)->cpp_obj);
			self->cpp_obj = new Grid2D(boundary2D);	
			Py_INCREF(boundary2D->getPyWrapper());
		}
		if(self->cpp_obj == NULL){
			ORBIT_MPI_Finalize("PyGrid2D - Grid2D(...) - constructor needs parameters.");
		}
		((Grid2D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The Grid2D __init__ has been called!"<<std::endl;
		return 0;
  }
	
	// setBoundary2D
  static PyObject* Grid2D_setBoundary2D(PyObject *self, PyObject *args){
		pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;	
		PyObject* pyB;
		if(!PyArg_ParseTuple(args,"O:setBoundary2D",&pyB)){		
			ORBIT_MPI_Finalize("PyGrid2D - setBoundary2D(Boundary2D) - the parameter is needed");
		}	
		PyObject* pyORBIT_Boindary2D_Type = getSpaceChargeType("Boundary2D");
		if(!PyObject_IsInstance(pyB,pyORBIT_Boindary2D_Type)){
			ORBIT_MPI_Finalize("PyPhaseVector - setBoundary2D(Boundary2D) - the parameter is needed");
		}	
		Boundary2D* boundary2D_init = cpp_Grid2D->getBoundary2D();
		Boundary2D* boundary2D = (Boundary2D*)(((pyORBIT_Object*) pyB)->cpp_obj);		
		if(boundary2D_init != NULL && boundary2D_init != boundary2D){
			Py_DECREF(boundary2D_init->getPyWrapper());
		}
		cpp_Grid2D->setBoundary2D(boundary2D);
		if(boundary2D_init != boundary2D) {Py_INCREF(boundary2D->getPyWrapper());}
		Py_INCREF(Py_None);
    return Py_None;	
	}

	//getBoundary2D
  static PyObject* Grid2D_getBoundary2D(PyObject *self, PyObject *args){
		pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		Boundary2D* boundary2D = cpp_Grid2D->getBoundary2D();
    return boundary2D->getPyWrapper();
	}	
	
	//setZero()
	static PyObject* Grid2D_setZero(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		cpp_Grid2D->setZero();
		Py_INCREF(Py_None);
    return Py_None;	
	}
	
	//getValue(double x, double y)
  static PyObject* Grid2D_getValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double x,y;
		if(!PyArg_ParseTuple(args,"dd:getValue",&x,&y)){
			ORBIT_MPI_Finalize("PyGrid2D - getValue(x,y) - parameters are needed.");
		}
		return Py_BuildValue("d",cpp_Grid2D->getValue(x,y));
	}
	
	//setValue(double value, int ix, int y)
  static PyObject* Grid2D_setValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double val;
		int ix,iy;
		if(!PyArg_ParseTuple(args,"dii:setValue",&val,&ix,&iy)){
			ORBIT_MPI_Finalize("PyGrid2D - setValue(val,ix,iy) - parameters are needed.");
		}
		return Py_BuildValue("d",cpp_Grid2D->getValue(ix,iy));
	}
	
	//setGridX(double min, double max, int n)
  static PyObject* Grid2D_setGridX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double min,max;
		int n;
		if(!PyArg_ParseTuple(args,"ddi:setGridX",&min,&max,&n)){
			ORBIT_MPI_Finalize("PyGrid2D - setGridX(min,max,n) - parameters are needed.");
		}
		cpp_Grid2D->setGridX(min,max,n);
		Py_INCREF(Py_None);
    return Py_None;	
	}
	
	//setGridY(double min, double max, int n)
  static PyObject* Grid2D_setGridY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double min,max;
		int n;
		if(!PyArg_ParseTuple(args,"ddi:setGridY",&min,&max,&n)){
			ORBIT_MPI_Finalize("PyGrid2D - setGridY(min,max,n) - parameters are needed.");
		}
		cpp_Grid2D->setGridY(min,max,n);
		Py_INCREF(Py_None);
    return Py_None;	
	}	
	
	//getGridX()
  static PyObject* Grid2D_getGridX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		//create tuple with names
		PyObject* resTuple = PyTuple_New(cpp_Grid2D->getBinsX());
		for(int i = 0, n = cpp_Grid2D->getBinsX(); i < n; i++){
			PyObject* py_p = Py_BuildValue("d",cpp_Grid2D->getGridX(i));
			if(PyTuple_SetItem(resTuple,i,py_p)){
				ORBIT_MPI_Finalize("PyBunch - getGridX() - cannot create tuple with grid X points");
			}
		}
		return resTuple;	
	}	
	
	//getGridY()
  static PyObject* Grid2D_getGridY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		//create tuple with names
		PyObject* resTuple = PyTuple_New(cpp_Grid2D->getBinsY());
		for(int i = 0, n = cpp_Grid2D->getBinsY(); i < n; i++){
			PyObject* py_p = Py_BuildValue("d",cpp_Grid2D->getGridY(i));
			if(PyTuple_SetItem(resTuple,i,py_p)){
				ORBIT_MPI_Finalize("PyBunch - getGridY() - cannot create tuple with grid Y points");
			}
		}
		return resTuple;	
	}		
	
	//binValue(double value, double x, double y)
  static PyObject* Grid2D_binValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double val,x,y;
		if(!PyArg_ParseTuple(args,"ddd:binValue",&val,&x,&y)){
			ORBIT_MPI_Finalize("PyGrid2D - binValue(val,x,y) - parameters are needed.");
		}
		cpp_Grid2D->binValue(val,x,y);
		Py_INCREF(Py_None);
    return Py_None;	
	}	

	//calcGradient(double x, double y)
  static PyObject* Grid2D_calcGradient(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double x,y;
		double ex, ey;
		if(!PyArg_ParseTuple(args,"dd:calcGradient",&x,&y)){
			ORBIT_MPI_Finalize("PyGrid2D - calcGradient(x,y) - parameters are needed.");
		}
		cpp_Grid2D->calcGradient(x,y,ex,ey);
		return Py_BuildValue("(dd)",ex,ey);
	}	
	
  //-----------------------------------------------------
  //destructor for python Grid2D class (__del__ method).
  //-----------------------------------------------------
  static void Grid2D_del(pyORBIT_Object* self){
		//std::cerr<<"The Grid2D __del__ has been called!"<<std::endl;
		Grid2D* cpp_Grid2D = (Grid2D*) self->cpp_obj;
		if(cpp_Grid2D->getBoundary2D() != NULL){
			Py_DECREF(cpp_Grid2D->getBoundary2D()->getPyWrapper());
		}
		delete (cpp_Grid2D);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python Grid2D wrapper class
	// they will be vailable from python level
  static PyMethodDef Grid2DClassMethods[] = {
		{ "getBoundary2D",   Grid2D_getBoundary2D, METH_VARARGS,"returns the internal Boundary2D instance"},
		{ "setBoundary2D",   Grid2D_setBoundary2D, METH_VARARGS,"sets the internal Boundary2D instance"},
		{ "setZero",         Grid2D_setZero,       METH_VARARGS,"sets all points on the grid to zero"},
		{ "getValue",        Grid2D_getValue,      METH_VARARGS,"returns value for (x,y) point"},
		{ "setValue",        Grid2D_setValue,      METH_VARARGS,"sets value for (ix,iy) point"},
		{ "setGridX",        Grid2D_setGridX,      METH_VARARGS,"sets the X grid with min,max, number of points"},
		{ "setGridY",        Grid2D_setGridY,      METH_VARARGS,"sets the Y grid with min,max, number of points"},
		{ "getGridX",        Grid2D_getGridX,      METH_VARARGS,"returns the X as a tuple"},
		{ "getGridY",        Grid2D_getGridY,      METH_VARARGS,"returns the Y as a tuple"},
		{ "binValue",        Grid2D_binValue,      METH_VARARGS,"bins the value into the 2D mesh"},
		{ "calcGradient",    Grid2D_calcGradient,  METH_VARARGS,"returns gradient as (gx,gy) for point (x,y)"},
    {NULL}
  };

	// defenition of the memebers of the python Grid2D wrapper class
	// they will be vailable from python level
	static PyMemberDef Grid2DClassMembers [] = {
		{NULL}
	};

	//new python Grid2D wrapper type definition
	static PyTypeObject pyORBIT_Grid2D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Grid2D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Grid2D_del , /*tp_dealloc*/
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
		"The Grid2D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		Grid2DClassMethods, /* tp_methods */
		Grid2DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Grid2D_init, /* tp_init */
		0, /* tp_alloc */
		Grid2D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyGrid2D class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initGrid2D(PyObject* module){
		if (PyType_Ready(&pyORBIT_Grid2D_Type) < 0) return;
		Py_INCREF(&pyORBIT_Grid2D_Type);
		PyModule_AddObject(module, "Grid2D", (PyObject *)&pyORBIT_Grid2D_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
