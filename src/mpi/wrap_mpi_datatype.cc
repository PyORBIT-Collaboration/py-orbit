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

#include "wrap_mpi_datatype.hh"

namespace wrap_orbit_mpi_datatype{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python MPI_Datatype class definition
	//---------------------------------------------------------

	//constructor for python class wrapping MPI_Datatype instance
	//It never will be called directly
	static PyObject* MPI_Datatype_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_MPI_Datatype* self;
		self = (pyORBIT_MPI_Datatype *) type->tp_alloc(type, 0);
		return (PyObject *) self;
	}

  //initializator for python MPI_Datatype  class
  //this is implementation of the __init__ method
  static int MPI_Datatype_init(pyORBIT_MPI_Datatype *self, PyObject *args, PyObject *kwds){
    //pyORBIT_MPI_Datatype* pyMPI_Datatype = (pyORBIT_MPI_Datatype*) self;		
    if(PyTuple_Size(args) != 0){
      error("MPI_Datatype constructor needs nothing.");
    }
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python MPI_Datatype class.
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static void MPI_Datatype_del(pyORBIT_MPI_Datatype* self){
		//std::cerr<<"The MPI_Datatype __del__ has been called!"<<std::endl;	
		self->ob_type->tp_free((PyObject*)self);
  }
		
	// defenition of the methods of the python MPI_Datatype wrapper class
	// they will be vailable from python level
  static PyMethodDef MPI_DatatypeClassMethods[] = {
    //{ "test",       MPI_Datatype_test      ,METH_VARARGS,"document string"},
    {NULL}
  };

	// defenition of the memebers of the python MPI_Datatype wrapper class
	// they will be vailable from python level
	static PyMemberDef MPI_DatatypeClassMembers[] = {
		{NULL}
	};
	
	//new python SyncPart wrapper type definition
	static PyTypeObject pyORBIT_MPI_Datatype_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"MPI_Datatype", /*tp_name*/
		sizeof(pyORBIT_MPI_Datatype), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) MPI_Datatype_del , /*tp_dealloc*/
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
		"The MPI_Datatype python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		MPI_DatatypeClassMethods, /* tp_methods */
		MPI_DatatypeClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) MPI_Datatype_init, /* tp_init */
		0, /* tp_alloc */
		MPI_Datatype_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the MPI_Datatype class
	//It will be called from orbit_mpi wrapper initialization
	//--------------------------------------------------
  void init_orbit_mpi_datatype(PyObject* module){
		if (PyType_Ready(&pyORBIT_MPI_Datatype_Type) < 0) return;
		Py_INCREF(&pyORBIT_MPI_Datatype_Type);
		
		PyObject * datatype_module = PyModule_New("mpi_datatype");
		PyModule_AddObject(datatype_module, "MPI_Datatype", (PyObject *)&pyORBIT_MPI_Datatype_Type);
		Py_INCREF(datatype_module);
		
		pyORBIT_MPI_Datatype* pyMPI_Datatype_CHAR = PyObject_New(pyORBIT_MPI_Datatype,&pyORBIT_MPI_Datatype_Type);
		pyMPI_Datatype_CHAR->datatype = MPI_CHAR;
		Py_INCREF((PyObject *) pyMPI_Datatype_CHAR); 
		
		pyORBIT_MPI_Datatype* pyMPI_Datatype_INT = PyObject_New(pyORBIT_MPI_Datatype,&pyORBIT_MPI_Datatype_Type);
		pyMPI_Datatype_INT->datatype = MPI_INT;
		Py_INCREF((PyObject *) pyMPI_Datatype_INT); 

		pyORBIT_MPI_Datatype* pyMPI_Datatype_DOUBLE = PyObject_New(pyORBIT_MPI_Datatype,&pyORBIT_MPI_Datatype_Type);
		pyMPI_Datatype_DOUBLE->datatype = MPI_DOUBLE;
		Py_INCREF((PyObject *) pyMPI_Datatype_DOUBLE); 

    PyModule_AddObject(datatype_module, "MPI_CHAR", (PyObject *) pyMPI_Datatype_CHAR);
    PyModule_AddObject(datatype_module, "MPI_INT", (PyObject *) pyMPI_Datatype_INT);
    PyModule_AddObject(datatype_module, "MPI_DOUBLE", (PyObject *) pyMPI_Datatype_DOUBLE);		
		
		PyModule_AddObject(module, "mpi_datatype", datatype_module);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_mpi_datatype
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
