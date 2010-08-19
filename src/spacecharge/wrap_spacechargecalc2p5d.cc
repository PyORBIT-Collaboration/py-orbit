#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#
#include "wrap_spacechargecalc2p5d.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "SpaceChargeCalc2p5D.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python SpaceChargeCalc2p5D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping SpaceChargeCalc2p5D instance
	//It never will be called directly

	static PyObject* SpaceChargeCalc2p5D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The SpaceChargeCalc2p5D new has been called!"<<std::endl;
		return (PyObject *) self;
	}
	
  //initializator for python SpaceChargeCalc2p5D class
  //this is implementation of the __init__ method
  static int SpaceChargeCalc2p5D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
  	int xSize,ySize,zSize;
  	double xy_ratio;
  	double xMin = -1.0, yMin = -1.0, xMax = +1.0, yMax = +1.0;
	//if nVars == 4 this is SpaceChargeCalc2p5D(int xSize, int ySize, int zSize, double xy_ratio)
	//if nVars == 7 this is SpaceChargeCalc2p5D(int xSize, int ySize, int zSize, double xMin, double xMax, double yMin, double yMax)
	int nVars = PyTuple_Size(args);
	if(nVars == 4){
		if(!PyArg_ParseTuple(args,"iiid:__init__",&xSize,&ySize,&zSize,&xy_ratio)){
		ORBIT_MPI_Finalize("PySpaceChargeCalc2p5D - SpaceChargeCalc2p5D(xSize,ySize,xzSize,xy_ratio]) - constructor needs parameters.");
		}
		self->cpp_obj = new SpaceChargeCalc2p5D(xSize,ySize,zSize,xy_ratio);	
	}
	if(nVars == 7){
		if(!PyArg_ParseTuple(args,"iiidddd:__init__",&xSize,&ySize,&zSize,&xMin,&xMax,&yMin,&yMax)){
		ORBIT_MPI_Finalize("PySpaceChargeCalc2p5D - SpaceChargeCalc2p5D(xSize,ySize,xzSize,xMin,xMax,yMin,yMax) - constructor needs parameters.");
		}
		self->cpp_obj = new SpaceChargeCalc2p5D(xSize,ySize,zSize,xMin,xMax,yMin,yMax);
	}
	((SpaceChargeCalc2p5D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
	//std::cerr<<"The SpaceChargeCalc2p5D __init__ has been called!"<<std::endl;
	return 0;
  }
  
  //trackBunch(Bunch* bunch, double length, BaseBoundary2D* boundary)
  static PyObject* SpaceChargeCalc2p5D_trackBunch(PyObject *self, PyObject *args){
  	  pyORBIT_Object* pySpaceChargeCalc2p5D = (pyORBIT_Object*) self;
  	  SpaceChargeCalc2p5D* cpp_SpaceChargeCalc2p5D = (SpaceChargeCalc2p5D*) pySpaceChargeCalc2p5D->cpp_obj;
  	  PyObject* pyBunch;
  	  PyObject* pyBoundary;
  	  double length;
  	  //std::cerr<<"trackBunch begin in wrapper.";
  	  int nVars = PyTuple_Size(args);
  	  if(nVars == 2 || nVars == 3){
  	    if(nVars == 2){
		if(!PyArg_ParseTuple(args,"Od:trackBunch",&pyBunch,&length)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc2p5D.trackBunch(pyBunch,length) - constructor needs parameters.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc2p5D.trackBunch(pyBunch,length) - method needs parameters.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;		
		cpp_SpaceChargeCalc2p5D->trackBunch(cpp_bunch,length);
	    } else{
		if(!PyArg_ParseTuple(args,"OdO:trackBunch",&pyBunch,&length,&pyBoundary)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc2p5D.trackBunch(pyBunch,length,pyBoundary) - constructor needs parameters.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		PyObject* pyORBIT_Boundary_Type = getSpaceChargeType("Boundary2D");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type) || !PyObject_IsInstance(pyBoundary,pyORBIT_Boundary_Type)){
			ORBIT_MPI_Finalize("PySpaceChargeCalc2p5D.trackBunch(pyBunch,pyBoundary,length) - method needs parameters.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		BaseBoundary2D* cpp_boundary = (BaseBoundary2D*) ((pyORBIT_Object*)pyBoundary)->cpp_obj;
		cpp_SpaceChargeCalc2p5D->trackBunch(cpp_bunch,length,cpp_boundary);
	    }  	
	} else{
	    ORBIT_MPI_Finalize("PyBoundary. You should call trackBunch(pyBunch,length) or trackBunch(pyBunch,pyBoundary,length)");    
	  }
	Py_INCREF(Py_None);
	return Py_None;
  }

  //-----------------------------------------------------
  //destructor for python SpaceChargeCalc2p5D class (__del__ method).
  //-----------------------------------------------------
  static void SpaceChargeCalc2p5D_del(pyORBIT_Object* self){
		//std::cerr<<"The SpaceChargeCalc2p5D __del__ has been called!"<<std::endl;
		SpaceChargeCalc2p5D* cpp_SpaceChargeCalc2p5D = (SpaceChargeCalc2p5D*) self->cpp_obj;
		if(cpp_SpaceChargeCalc2p5D != NULL){
			delete cpp_SpaceChargeCalc2p5D;
		}
		self->ob_type->tp_free((PyObject*)self);
  }	
  
  // defenition of the methods of the python SpaceChargeCalc2p5D wrapper class
  // they will be vailable from python level
  static PyMethodDef SpaceChargeCalc2p5DClassMethods[] = {
		{ "trackBunch", SpaceChargeCalc2p5D_trackBunch, METH_VARARGS,"track the Bunch"},
		{NULL}
  };
  
  // defenition of the memebers of the python SpaceChargeCalc2p5D wrapper class
  // they will be vailable from python level
  static PyMemberDef SpaceChargeCalc2p5DClassMembers [] = {
		{NULL}
  };

	//new python SpaceChargeCalc2p5D wrapper type definition
	static PyTypeObject pyORBIT_SpaceChargeCalc2p5D_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SpaceChargeCalc2p5D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SpaceChargeCalc2p5D_del , /*tp_dealloc*/
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
		"The SpaceChargeCalc2p5D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SpaceChargeCalc2p5DClassMethods, /* tp_methods */
		SpaceChargeCalc2p5DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SpaceChargeCalc2p5D_init, /* tp_init */
		0, /* tp_alloc */
		SpaceChargeCalc2p5D_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pySpaceChargeCalc2p5D class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initSpaceChargeCalc2p5D(PyObject* module){
		if (PyType_Ready(&pyORBIT_SpaceChargeCalc2p5D_Type) < 0) return;
		Py_INCREF(&pyORBIT_SpaceChargeCalc2p5D_Type);
		PyModule_AddObject(module, "SpaceChargeCalc2p5D", (PyObject *)&pyORBIT_SpaceChargeCalc2p5D_Type);
		//std::cout<<"debug SpaceChargeCalc2p5D added! "<<std::endl;
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
