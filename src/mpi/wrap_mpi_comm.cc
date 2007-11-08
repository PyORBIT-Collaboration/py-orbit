///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"

//this header is from Python package
#include "structmember.h"

//c++ header for cerr and cout
#include <iostream>

#include "wrap_mpi_comm.hh"

namespace wrap_orbit_mpi_comm{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python MPI_Comm class definition
	//---------------------------------------------------------

	//constructor for python class wrapping MPI_Comm instance
	//It never will be called directly
	static PyObject* MPI_Comm_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_MPI_Comm* self;
		self = (pyORBIT_MPI_Comm *) type->tp_alloc(type, 0);
		self->comm = MPI_COMM_WORLD;
		return (PyObject *) self;
	}

  //initializator for python MPI_Comm  class
  //this is implementation of the __init__ method
  static int MPI_Comm_init(pyORBIT_MPI_Comm *self, PyObject *args, PyObject *kwds){
    pyORBIT_MPI_Comm* pyMPI_Comm = (pyORBIT_MPI_Comm*) self;
	  int nArgs = PyTuple_Size(args);
		
    if(nArgs != 0 && nArgs != 1){
      error("MPI_Comm constructor needs nothing or MPI_Comm as input parameter.");
    }
		
		if(nArgs == 1){
			pyORBIT_MPI_Comm* pyMPI_Comm_in = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
			pyMPI_Comm->comm = pyMPI_Comm_in->comm;
		}
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python MPI_Comm class.
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static void MPI_Comm_del(pyORBIT_MPI_Comm* self){
		//std::cerr<<"The MPI_Comm __del__ has been called!"<<std::endl;
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python MPI_Comm wrapper class
	// they will be vailable from python level
  static PyMethodDef MPI_CommClassMethods[] = {
    //{ "test",       MPI_Comm_test      ,METH_VARARGS,"document string"},
    {NULL}
  };

	// defenition of the memebers of the python MPI_Comm wrapper class
	// they will be vailable from python level
	static PyMemberDef MPI_CommClassMembers[] = {
		{NULL}
	};
	
	//new python SyncPart wrapper type definition
	static PyTypeObject pyORBIT_MPI_Comm_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"MPI_Comm", /*tp_name*/
		sizeof(pyORBIT_MPI_Comm), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) MPI_Comm_del , /*tp_dealloc*/
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
		"The MPI_Comm python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		MPI_CommClassMethods, /* tp_methods */
		MPI_CommClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) MPI_Comm_init, /* tp_init */
		0, /* tp_alloc */
		MPI_Comm_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the MPI_Comm class
	//It will be called from orbit_mpi wrapper initialization
	//--------------------------------------------------
  void init_orbit_mpi_comm(PyObject* module){
		if (PyType_Ready(&pyORBIT_MPI_Comm_Type) < 0) return;
		Py_INCREF(&pyORBIT_MPI_Comm_Type);
		
		PyObject * comm_module = PyModule_New("mpi_comm");
		PyModule_AddObject(comm_module, "MPI_Comm", (PyObject *)&pyORBIT_MPI_Comm_Type);
		Py_INCREF(comm_module);

		pyORBIT_MPI_Comm* pyMPI_Comm_WORLD = PyObject_New(pyORBIT_MPI_Comm,&pyORBIT_MPI_Comm_Type);
		pyMPI_Comm_WORLD->comm = MPI_COMM_WORLD;
		Py_INCREF((PyObject *) pyMPI_Comm_WORLD);
		
		pyORBIT_MPI_Comm* pyMPI_Comm_SELF = PyObject_New(pyORBIT_MPI_Comm,&pyORBIT_MPI_Comm_Type);
		pyMPI_Comm_SELF->comm = MPI_COMM_SELF;
		Py_INCREF((PyObject *) pyMPI_Comm_SELF); 

    PyModule_AddObject(comm_module, "MPI_COMM_WORLD", (PyObject *) pyMPI_Comm_WORLD);
    PyModule_AddObject(comm_module, "MPI_COMM_SELF", (PyObject *) pyMPI_Comm_SELF);
		
		PyModule_AddObject(module, "mpi_comm", comm_module);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_mpi_comm
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
