#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#
#include "wrap_envsolver_kv.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "EnvSolverKV.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

    //---------------------------------------------------------
    //Python EnvSolverKV class definition
    //---------------------------------------------------------

    //constructor for python class wrapping EnvSolverKV instance
    //It never will be called directly

    static PyObject* EnvSolverKV_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
        pyORBIT_Object* self;
        self = (pyORBIT_Object *) type->tp_alloc(type, 0);
        self->cpp_obj = NULL;
        //std::cerr<<"The EnvSolverKV new has been called!"<<std::endl;
        return (PyObject *) self;
    }
    
  //initializator for python EnvSolverKV class
  //this is implementation of the __init__ method EnvSolverKV(double ex, double ey, double perveance)
  static int EnvSolverKV_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
  {
      double ex, ey, Q;
        if (!PyArg_ParseTuple(args, "ddd:__init__", &ex, &ey, &Q)) {
            ORBIT_MPI_Finalize("PyEnvSolverKV - EnvSolverKV(ex, ey, Q) - constructor needs parameters.");
        }
      self->cpp_obj = new EnvSolverKV(ex, ey, Q);
      ((EnvSolverKV*) self->cpp_obj)->setPyWrapper((PyObject*) self);
      //std::cerr<<"The EnvSolverKV __init__ has been called!"<<std::endl;
      return 0;
    }
    
    //trackBunch(Bunch* bunch, double length)
    static PyObject* EnvSolverKV_trackBunch(PyObject *self, PyObject *args) {
        int nVars = PyTuple_Size(args);
        pyORBIT_Object* pyEnvSolverKV = (pyORBIT_Object*) self;
        EnvSolverKV* cpp_EnvSolverKV = (EnvSolverKV*) pyEnvSolverKV->cpp_obj;
        PyObject* pyBunch;
        double length;
        
        if (!PyArg_ParseTuple(args, "Od:trackBunch", &pyBunch, &length)){
            ORBIT_MPI_Finalize("PyEnvSolverKV.trackBunch(pyBunch, length) - method needs parameters.");
        }
        PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
        if (!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
            ORBIT_MPI_Finalize("PyEnvSolverKV.trackBunch(pyBunch,length) - pyBunch is not Bunch.");
        }
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
        cpp_EnvSolverKV->trackBunch(cpp_bunch, length);
        Py_INCREF(Py_None);
        return Py_None;
  }

  //-----------------------------------------------------
  //destructor for python EnvSolverKV class (__del__ method).
  //-----------------------------------------------------
  static void EnvSolverKV_del(pyORBIT_Object* self){
        EnvSolverKV* cpp_EnvSolverKV = (EnvSolverKV*) self->cpp_obj;
        if(cpp_EnvSolverKV != NULL){
            delete cpp_EnvSolverKV;
        }
        self->ob_type->tp_free((PyObject*)self);
  }
  
  // defenition of the methods of the python EnvSolverKV wrapper class
  // they will be vailable from python level
  static PyMethodDef EnvSolverKVClassMethods[] = {
        { "trackBunch",  EnvSolverKV_trackBunch, METH_VARARGS,"track the bunch - trackBunch(pyBunch, length)"},
        {NULL}
  };
  
  // defenition of the memebers of the python EnvSolverKV wrapper class
  // they will be vailable from python level
  static PyMemberDef EnvSolverKVClassMembers [] = {
        {NULL}
  };

    //new python EnvSolverKV wrapper type definition
    static PyTypeObject pyORBIT_EnvSolverKV_Type = {
        PyObject_HEAD_INIT(NULL)
        0, /*ob_size*/
        "EnvSolverKV", /*tp_name*/
        sizeof(pyORBIT_Object), /*tp_basicsize*/
        0, /*tp_itemsize*/
        (destructor) EnvSolverKV_del , /*tp_dealloc*/
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
        "The EnvSolverKV python wrapper", /* tp_doc */
        0, /* tp_traverse */
        0, /* tp_clear */
        0, /* tp_richcompare */
        0, /* tp_weaklistoffset */
        0, /* tp_iter */
        0, /* tp_iternext */
        EnvSolverKVClassMethods, /* tp_methods */
        EnvSolverKVClassMembers, /* tp_members */
        0, /* tp_getset */
        0, /* tp_base */
        0, /* tp_dict */
        0, /* tp_descr_get */
        0, /* tp_descr_set */
        0, /* tp_dictoffset */
        (initproc) EnvSolverKV_init, /* tp_init */
        0, /* tp_alloc */
        EnvSolverKV_new, /* tp_new */
    };

    //--------------------------------------------------
    //Initialization function of the pyEnvSolverKV class
    //It will be called from SpaceCharge wrapper initialization
    //--------------------------------------------------
  void initEnvSolverKV(PyObject* module){
        if (PyType_Ready(&pyORBIT_EnvSolverKV_Type) < 0) return;
        Py_INCREF(&pyORBIT_EnvSolverKV_Type);
        PyModule_AddObject(module, "EnvSolverKV", (PyObject *)&pyORBIT_EnvSolverKV_Type);
    }

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
