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

#include "wrap_mpi_status.hh"

namespace wrap_orbit_mpi_status{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python MPI_Status class definition
	//---------------------------------------------------------

	//constructor for python class wrapping MPI_Status instance
	//It never will be called directly
	static PyObject* MPI_Status_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_MPI_Status* self;
		self = (pyORBIT_MPI_Status *) type->tp_alloc(type, 0);
		return (PyObject *) self;
	}

  //initializator for python MPI_Status  class
  //this is implementation of the __init__ method
  static int MPI_Status_init(pyORBIT_MPI_Status *self, PyObject *args, PyObject *kwds){
    //pyORBIT_MPI_Status* pyMPI_Status = (pyORBIT_MPI_Status*) self;		
    if(PyTuple_Size(args) != 0){
      error("MPI_Status constructor needs nothing.");
    }
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python MPI_Status class.
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static void MPI_Status_del(pyORBIT_MPI_Status* self){
		//std::cerr<<"The MPI_Status __del__ has been called!"<<std::endl;	
		self->ob_type->tp_free((PyObject*)self);
  }
		
	// defenition of the methods of the python MPI_Status wrapper class
	// they will be vailable from python level
  static PyMethodDef MPI_StatusClassMethods[] = {
    //{ "test",       MPI_Status_test      ,METH_VARARGS,"document string"},
    {NULL}
  };

	// defenition of the memebers of the python MPI_Status wrapper class
	// they will be vailable from python level
	static PyMemberDef MPI_StatusClassMembers[] = {
		{NULL}
	};
	
	//new python SyncPart wrapper type definition
	static PyTypeObject pyORBIT_MPI_Status_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"MPI_Status", /*tp_name*/
		sizeof(pyORBIT_MPI_Status), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) MPI_Status_del , /*tp_dealloc*/
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
		"The MPI_Status python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		MPI_StatusClassMethods, /* tp_methods */
		MPI_StatusClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) MPI_Status_init, /* tp_init */
		0, /* tp_alloc */
		MPI_Status_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the MPI_Status class
	//It will be called from orbit_mpi wrapper initialization
	//--------------------------------------------------
  void init_orbit_mpi_status(PyObject* module){
		if (PyType_Ready(&pyORBIT_MPI_Status_Type) < 0) return;
		Py_INCREF(&pyORBIT_MPI_Status_Type);
		
		PyObject * status_module = PyModule_New("mpi_status");
		PyModule_AddObject(status_module, "MPI_Status", (PyObject *)&pyORBIT_MPI_Status_Type);
		Py_INCREF(status_module);
		
		PyModule_AddObject(module, "mpi_status", status_module);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_mpi_status
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
