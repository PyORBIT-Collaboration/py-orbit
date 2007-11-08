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

#include "wrap_mpi_op.hh"

namespace wrap_orbit_mpi_op{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python MPI_Op class definition
	//---------------------------------------------------------

	//constructor for python class wrapping MPI_Op instance
	//It never will be called directly
	static PyObject* MPI_Op_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_MPI_Op* self;
		self = (pyORBIT_MPI_Op *) type->tp_alloc(type, 0);
		return (PyObject *) self;
	}

  //initializator for python MPI_Op  class
  //this is implementation of the __init__ method
  static int MPI_Op_init(pyORBIT_MPI_Op *self, PyObject *args, PyObject *kwds){
    //pyORBIT_MPI_Op* pyMPI_Op = (pyORBIT_MPI_Op*) self;		
    if(PyTuple_Size(args) != 0){
      error("MPI_Op constructor needs nothing.");
    }
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python MPI_Op class.
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static void MPI_Op_del(pyORBIT_MPI_Op* self){
		//std::cerr<<"The MPI_Op __del__ has been called!"<<std::endl;	
		self->ob_type->tp_free((PyObject*)self);
  }
		
	// defenition of the methods of the python MPI_Op wrapper class
	// they will be vailable from python level
  static PyMethodDef MPI_OpClassMethods[] = {
    //{ "test",       MPI_Op_test      ,METH_VARARGS,"document string"},
    {NULL}
  };

	// defenition of the memebers of the python MPI_Op wrapper class
	// they will be vailable from python level
	static PyMemberDef MPI_OpClassMembers[] = {
		{NULL}
	};
	
	//new python SyncPart wrapper type definition
	static PyTypeObject pyORBIT_MPI_Op_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"MPI_Op", /*tp_name*/
		sizeof(pyORBIT_MPI_Op), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) MPI_Op_del , /*tp_dealloc*/
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
		"The MPI_Op python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		MPI_OpClassMethods, /* tp_methods */
		MPI_OpClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) MPI_Op_init, /* tp_init */
		0, /* tp_alloc */
		MPI_Op_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the MPI_Op class
	//It will be called from orbit_mpi wrapper initialization
	//--------------------------------------------------
  void init_orbit_mpi_op(PyObject* module){
		if (PyType_Ready(&pyORBIT_MPI_Op_Type) < 0) return;
		Py_INCREF(&pyORBIT_MPI_Op_Type);
		
		PyObject * op_module = PyModule_New("mpi_op");
		PyModule_AddObject(op_module, "MPI_Op", (PyObject *)&pyORBIT_MPI_Op_Type);
		Py_INCREF(op_module);
		
		pyORBIT_MPI_Op* pyMPI_Op_MAX = PyObject_New(pyORBIT_MPI_Op,&pyORBIT_MPI_Op_Type);
		pyMPI_Op_MAX->op = MPI_MAX;
		Py_INCREF((PyObject *) pyMPI_Op_MAX); 
    PyModule_AddObject(op_module, "MPI_MAX", (PyObject *) pyMPI_Op_MAX);		
		
		pyORBIT_MPI_Op* pyMPI_Op_MIN = PyObject_New(pyORBIT_MPI_Op,&pyORBIT_MPI_Op_Type);
		pyMPI_Op_MIN->op = MPI_MIN;
		Py_INCREF((PyObject *) pyMPI_Op_MIN); 
    PyModule_AddObject(op_module, "MPI_MIN", (PyObject *) pyMPI_Op_MIN);	
		
		pyORBIT_MPI_Op* pyMPI_Op_SUM = PyObject_New(pyORBIT_MPI_Op,&pyORBIT_MPI_Op_Type);
		pyMPI_Op_SUM->op = MPI_SUM;
		Py_INCREF((PyObject *) pyMPI_Op_SUM); 
    PyModule_AddObject(op_module, "MPI_SUM", (PyObject *) pyMPI_Op_SUM);		
		
		pyORBIT_MPI_Op* pyMPI_Op_PROD = PyObject_New(pyORBIT_MPI_Op,&pyORBIT_MPI_Op_Type);
		pyMPI_Op_PROD->op = MPI_PROD;
		Py_INCREF((PyObject *) pyMPI_Op_PROD); 
    PyModule_AddObject(op_module, "MPI_PROD", (PyObject *) pyMPI_Op_PROD);			
		
		PyModule_AddObject(module, "mpi_op", op_module);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_mpi_op
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
