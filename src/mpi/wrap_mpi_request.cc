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

#include "wrap_mpi_request.hh"

namespace wrap_orbit_mpi_request{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python MPI_Request class definition
	//---------------------------------------------------------

	//constructor for python class wrapping MPI_Request instance
	//It never will be called directly
	static PyObject* MPI_Request_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_MPI_Request* self;
		self = (pyORBIT_MPI_Request *) type->tp_alloc(type, 0);
		return (PyObject *) self;
	}

  //initializator for python MPI_Request  class
  //this is implementation of the __init__ method
  static int MPI_Request_init(pyORBIT_MPI_Request *self, PyObject *args, PyObject *kwds){
    //pyORBIT_MPI_Request* pyMPI_Request = (pyORBIT_MPI_Request*) self;		
    if(PyTuple_Size(args) != 0){
      error("MPI_Request constructor needs nothing.");
    }
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python MPI_Request class.
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static void MPI_Request_del(pyORBIT_MPI_Request* self){
		//std::cerr<<"The MPI_Request __del__ has been called!"<<std::endl;	
		self->ob_type->tp_free((PyObject*)self);
  }
		
	// defenition of the methods of the python MPI_Request wrapper class
	// they will be vailable from python level
  static PyMethodDef MPI_RequestClassMethods[] = {
    //{ "test",       MPI_Request_test      ,METH_VARARGS,"document string"},
    {NULL}
  };

	// defenition of the memebers of the python MPI_Request wrapper class
	// they will be vailable from python level
	static PyMemberDef MPI_RequestClassMembers[] = {
		{NULL}
	};
	
	//new python SyncPart wrapper type definition
	static PyTypeObject pyORBIT_MPI_Request_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"MPI_Request", /*tp_name*/
		sizeof(pyORBIT_MPI_Request), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) MPI_Request_del , /*tp_dealloc*/
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
		"The MPI_Request python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		MPI_RequestClassMethods, /* tp_methods */
		MPI_RequestClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) MPI_Request_init, /* tp_init */
		0, /* tp_alloc */
		MPI_Request_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the MPI_Request class
	//It will be called from orbit_mpi wrapper initialization
	//--------------------------------------------------
  void init_orbit_mpi_request(PyObject* module){
		if (PyType_Ready(&pyORBIT_MPI_Request_Type) < 0) return;
		Py_INCREF(&pyORBIT_MPI_Request_Type);
		
		PyObject * request_module = PyModule_New("mpi_request");
		PyModule_AddObject(request_module, "MPI_Request", (PyObject *)&pyORBIT_MPI_Request_Type);
		Py_INCREF(request_module);
		
		PyModule_AddObject(module, "mpi_request", request_module);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_mpi_request
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
