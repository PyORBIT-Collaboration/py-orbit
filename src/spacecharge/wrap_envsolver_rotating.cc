#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#
#include "wrap_envsolver_kv.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "EnvSolverRotating.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

    //---------------------------------------------------------
    //Python EnvSolverRotating class definition
    //---------------------------------------------------------

    //constructor for python class wrapping EnvSolverRotating instance
    //It never will be called directly

    static PyObject* EnvSolverRotating_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
        pyORBIT_Object* self;
        self = (pyORBIT_Object *) type->tp_alloc(type, 0);
        self->cpp_obj = NULL;
        //std::cerr<<"The EnvSolverRotating new has been called!"<<std::endl;
        return (PyObject *) self;
    }
    
  //initializator for python EnvSolverRotating class
  //this is implementation of the __init__ method EnvSolverRotating(double perveance)
  static int EnvSolverRotating_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
  {
      double ex, ey, Q;
        if (!PyArg_ParseTuple(args, "d:__init__", &Q)) {
            ORBIT_MPI_Finalize("PyEnvSolverRotating - EnvSolverRotating(Q) - constructor needs parameters.");
        }
      self->cpp_obj = new EnvSolverRotating(Q);
      ((EnvSolverRotating*) self->cpp_obj)->setPyWrapper((PyObject*) self);
      //std::cerr<<"The EnvSolverRotating __init__ has been called!"<<std::endl;
      return 0;
    }
    
    //trackBunch(Bunch* bunch, double length)
    static PyObject* EnvSolverRotating_trackBunch(PyObject *self, PyObject *args) {
        int nVars = PyTuple_Size(args);
        pyORBIT_Object* pyEnvSolverRotating = (pyORBIT_Object*) self;
        EnvSolverRotating* cpp_EnvSolverRotating = (EnvSolverRotating*) pyEnvSolverRotating->cpp_obj;
        PyObject* pyBunch;
        double length;
        
        if (!PyArg_ParseTuple(args, "Od:trackBunch", &pyBunch, &length)){
            ORBIT_MPI_Finalize("PyEnvSolverRotating.trackBunch(pyBunch, length) - method needs parameters.");
        }
        PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
        if (!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
            ORBIT_MPI_Finalize("PyEnvSolverRotating.trackBunch(pyBunch,length) - pyBunch is not Bunch.");
        }
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
        cpp_EnvSolverRotating->trackBunch(cpp_bunch, length);
        Py_INCREF(Py_None);
        return Py_None;
  }

  //-----------------------------------------------------
  //destructor for python EnvSolverRotating class (__del__ method).
  //-----------------------------------------------------
  static void EnvSolverRotating_del(pyORBIT_Object* self){
        EnvSolverRotating* cpp_EnvSolverRotating = (EnvSolverRotating*) self->cpp_obj;
        if(cpp_EnvSolverRotating != NULL){
            delete cpp_EnvSolverRotating;
        }
        self->ob_type->tp_free((PyObject*)self);
  }
  
  // defenition of the methods of the python EnvSolverRotating wrapper class
  // they will be vailable from python level
  static PyMethodDef EnvSolverRotatingClassMethods[] = {
        { "trackBunch",  EnvSolverRotating_trackBunch, METH_VARARGS,"track the bunch - trackBunch(pyBunch, length)"},
        {NULL}
  };
  
  // defenition of the memebers of the python EnvSolverRotating wrapper class
  // they will be vailable from python level
  static PyMemberDef EnvSolverRotatingClassMembers [] = {
        {NULL}
  };

    //new python EnvSolverRotating wrapper type definition
    static PyTypeObject pyORBIT_EnvSolverRotating_Type = {
        PyObject_HEAD_INIT(NULL)
        0, /*ob_size*/
        "EnvSolverRotating", /*tp_name*/
        sizeof(pyORBIT_Object), /*tp_basicsize*/
        0, /*tp_itemsize*/
        (destructor) EnvSolverRotating_del , /*tp_dealloc*/
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
        "The EnvSolverRotating python wrapper", /* tp_doc */
        0, /* tp_traverse */
        0, /* tp_clear */
        0, /* tp_richcompare */
        0, /* tp_weaklistoffset */
        0, /* tp_iter */
        0, /* tp_iternext */
        EnvSolverRotatingClassMethods, /* tp_methods */
        EnvSolverRotatingClassMembers, /* tp_members */
        0, /* tp_getset */
        0, /* tp_base */
        0, /* tp_dict */
        0, /* tp_descr_get */
        0, /* tp_descr_set */
        0, /* tp_dictoffset */
        (initproc) EnvSolverRotating_init, /* tp_init */
        0, /* tp_alloc */
        EnvSolverRotating_new, /* tp_new */
    };

    //--------------------------------------------------
    //Initialization function of the pyEnvSolverRotating class
    //It will be called from SpaceCharge wrapper initialization
    //--------------------------------------------------
  void initEnvSolverRotating(PyObject* module){
        if (PyType_Ready(&pyORBIT_EnvSolverRotating_Type) < 0) return;
        Py_INCREF(&pyORBIT_EnvSolverRotating_Type);
        PyModule_AddObject(module, "EnvSolverRotating", (PyObject *)&pyORBIT_EnvSolverRotating_Type);
    }

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
