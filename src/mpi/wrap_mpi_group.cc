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

#include "wrap_mpi_group.hh"

namespace wrap_orbit_mpi_group{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif
	//---------------------------------------------------------
	//Python MPI_Group class definition
	//---------------------------------------------------------

	//constructor for python class wrapping MPI_Group instance
	//It never will be called directly
	static PyObject* MPI_Group_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_MPI_Group* self;
		self = (pyORBIT_MPI_Group *) type->tp_alloc(type, 0);
		self->group = MPI_GROUP_EMPTY;
		return (PyObject *) self;
	}

  //initializator for python MPI_Group  class
  //this is implementation of the __init__ method
  static int MPI_Group_init(pyORBIT_MPI_Group *self, PyObject *args, PyObject *kwds){
    pyORBIT_MPI_Group* pyMPI_Group = (pyORBIT_MPI_Group*) self;
	  int nArgs = PyTuple_Size(args);
		
    if(nArgs != 0 && nArgs != 1){
      error("MPI_Group constructor needs nothing or MPI_Group as input parameter.");
    }
		
		if(nArgs == 1){
			pyORBIT_MPI_Group* pyMPI_Group_in = (pyORBIT_MPI_Group*) PyTuple_GetItem(args,0);
			pyMPI_Group->group = pyMPI_Group_in->group;
		}
    return 0;
  }

  //-----------------------------------------------------
  //destructor for python MPI_Group class.
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static void MPI_Group_del(pyORBIT_MPI_Group* self){
		//std::cerr<<"The MPI_Group __del__ has been called!"<<std::endl;
		MPI_Comm group = self->group;
		if(group != MPI_GROUP_NULL && group != MPI_GROUP_EMPTY){
			ORBIT_MPI_Group_free(&group);
		}		
		self->ob_type->tp_free((PyObject*)self);
  }
		
	// defenition of the methods of the python MPI_Group wrapper class
	// they will be vailable from python level
  static PyMethodDef MPI_GroupClassMethods[] = {
    //{ "test",       MPI_Group_test      ,METH_VARARGS,"document string"},
    {NULL}
  };

	// defenition of the memebers of the python MPI_Group wrapper class
	// they will be vailable from python level
	static PyMemberDef MPI_GroupClassMembers[] = {
		{NULL}
	};
	
	//new python SyncPart wrapper type definition
	static PyTypeObject pyORBIT_MPI_Group_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"MPI_Group", /*tp_name*/
		sizeof(pyORBIT_MPI_Group), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) MPI_Group_del , /*tp_dealloc*/
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
		"The MPI_Group python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		MPI_GroupClassMethods, /* tp_methods */
		MPI_GroupClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) MPI_Group_init, /* tp_init */
		0, /* tp_alloc */
		MPI_Group_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the MPI_Group class
	//It will be called from orbit_mpi wrapper initialization
	//--------------------------------------------------
  void init_orbit_mpi_group(PyObject* module){
		if (PyType_Ready(&pyORBIT_MPI_Group_Type) < 0) return;
		Py_INCREF(&pyORBIT_MPI_Group_Type);
		
		PyObject * group_module = PyModule_New("mpi_group");
		PyModule_AddObject(group_module, "MPI_Group", (PyObject *)&pyORBIT_MPI_Group_Type);
		Py_INCREF(group_module);

		pyORBIT_MPI_Group* pyMPI_Group_NULL = PyObject_New(pyORBIT_MPI_Group,&pyORBIT_MPI_Group_Type);
		pyMPI_Group_NULL->group = MPI_GROUP_NULL;
		Py_INCREF((PyObject *) pyMPI_Group_NULL);
		
		pyORBIT_MPI_Group* pyMPI_Group_EMPTY = PyObject_New(pyORBIT_MPI_Group,&pyORBIT_MPI_Group_Type);
		pyMPI_Group_EMPTY->group = MPI_GROUP_EMPTY;
		Py_INCREF((PyObject *) pyMPI_Group_EMPTY); 

    PyModule_AddObject(group_module, "MPI_GROUP_NULL", (PyObject *) pyMPI_Group_NULL);
    PyModule_AddObject(group_module, "MPI_GROUP_EMPTY", (PyObject *) pyMPI_Group_EMPTY);
		
		PyModule_AddObject(module, "mpi_group", group_module);
	}
	
	//-----------------------------------------------------------
	//The function that will be exposed as C/C++ API for MPI_Group
	//-----------------------------------------------------------
	pyORBIT_MPI_Group* newMPI_Group(){
		pyORBIT_MPI_Group* pyMPI_Group = PyObject_New(pyORBIT_MPI_Group,&pyORBIT_MPI_Group_Type);
		pyMPI_Group->group = MPI_GROUP_EMPTY;
		Py_INCREF((PyObject *) pyMPI_Group);
    return pyMPI_Group;
	}
	
	void freeMPI_Group(pyORBIT_MPI_Group* pyMPI_Group){
		Py_DECREF(pyMPI_Group);
	}	

#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_mpi_group
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
