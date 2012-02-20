#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_foil.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "Foil.hh"

namespace wrap_foil{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	Constructor for python class wrapping c++ Foil instance.
      It never will be called directly.
	*/
	static PyObject* Foil_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int Foil_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){

	  double xmin = 0.; double xmax = 0.; double ymin = 0.; double ymax = 0.; double thick = 0.;
	  
	  //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
	  if(!PyArg_ParseTuple(	args,"ddddd:arguments",&xmin,&xmax,&ymin,&ymax,&thick)){
		  error("PyBunch - addParticle - cannot parse arguments! It should be (xmin, xmax, ymin, ymax, thick)");	
	  }
		self->cpp_obj =  new Foil(xmin, xmax, ymin, ymax, thick);
	  ((Foil*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
  /** Performs the full model foil scattering of the bunch */
  static PyObject* Foil_traverseFoilFullScatter(PyObject *self, PyObject *args){
	  Foil* cpp_Foil = (Foil*)((pyORBIT_Object*) self)->cpp_obj;
		PyObject* pyBunch;
		PyObject* pyLostBunch;
		if(!PyArg_ParseTuple(args,"OO:traverseFoilFullScatter",&pyBunch, &pyLostBunch)){
			ORBIT_MPI_Finalize("Foil - traverseFoilFullScatter(Bunch* bunch, Bunch* bunch) - parameters are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type) || !PyObject_IsInstance(pyLostBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("Foil - traverseFoilFullScatter(Bunch* bunch, Bunch* bunch) - method needs a Bunch.");
		}
	  
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		Bunch* cpp_lostbunch = (Bunch*) ((pyORBIT_Object*)pyLostBunch)->cpp_obj;
		cpp_Foil->traverseFoilFullScatter(cpp_bunch, cpp_lostbunch);
		Py_INCREF(Py_None);
		return Py_None;
  }


	/** Performs the simplified model foil scattering of the bunch */
	static PyObject* Foil_traverseFoilSimpleScatter(PyObject *self, PyObject *args){
		Foil* cpp_Foil = (Foil*)((pyORBIT_Object*) self)->cpp_obj;
		PyObject* pyBunch;
		if(!PyArg_ParseTuple(args,"O:traverseFoilFullScatter",&pyBunch)){
			ORBIT_MPI_Finalize("Foil - traverseFoilSimpleScatter(Bunch* bunch) - parameters are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("Foil - traverseFoilSimpleScatter(Bunch* bunch) - method needs a Bunch.");
		}
		
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		cpp_Foil->traverseFoilSimpleScatter(cpp_bunch);
		Py_INCREF(Py_None);
		return Py_None;
	}
	
		
	
  //-----------------------------------------------------
  //destructor for python Foil class (__del__ method).
  //-----------------------------------------------------
  static void Foil_del(pyORBIT_Object* self){
		//std::cerr<<"The Foil __del__ has been called!"<<std::endl;
		delete ((Foil*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// definition of the methods of the python Foil wrapper class
	// they will be vailable from python level
	static PyMethodDef FoilClassMethods[] = {
		{ "traverseFoilFullScatter", Foil_traverseFoilFullScatter, METH_VARARGS,"Performs the foil scatter of the bunch."},
		{ "traverseFoilSimpleScatter", Foil_traverseFoilSimpleScatter, METH_VARARGS,"Performs the foil scatter of the bunch."},
   {NULL}
  	};

	static PyMemberDef FoilClassMembers [] = {
		{NULL}
	};
	
	
	//new python Foil wrapper type definition
	static PyTypeObject pyORBIT_Foil_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Foil", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Foil_del , /*tp_dealloc*/
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
		"The Foil python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		FoilClassMethods, /* tp_methods */
		FoilClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Foil_init, /* tp_init */
		0, /* tp_alloc */
		Foil_new, /* tp_new */
	};	
	
	
	
	//--------------------------------------------------
	//Initialization Foil of the pyBunchFoil class
	//--------------------------------------------------

	void initfoil(){
		//check that the Foil wrapper is ready
		if (PyType_Ready(&pyORBIT_Foil_Type) < 0) return;
		Py_INCREF(&pyORBIT_Foil_Type);
		//create new module
		PyObject* module = Py_InitModule("foil",FoilClassMethods);
		PyModule_AddObject(module, "Foil", (PyObject *)&pyORBIT_Foil_Type);			
	}

#ifdef __cplusplus
}
#endif


}
