#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_LImpedance.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "LImpedance.hh"

using namespace OrbitUtils;

namespace wrap_impedances
{

#ifdef __cplusplus
extern "C"
{
#endif

  //---------------------------------------------------------
  // Python LImpedance class definition
  //---------------------------------------------------------

  //-----------------------------------------------------
  // Constructor for python class wrapping LImpedance instance
  // It never will be called directly
  //-----------------------------------------------------

  static PyObject* LImpedance_new(PyTypeObject *type,
                                  PyObject *args,
                                  PyObject *kwds)
  {
    pyORBIT_Object* self;
    self = (pyORBIT_Object *) type->tp_alloc(type, 0);
    self->cpp_obj = NULL;
    // std::cerr << "The LImpedance new has been called!" << std::endl;
    return (PyObject *) self;
  }

  //-----------------------------------------------------
  // Initialization of python LImpedance class
  // This is implementation of the
  // __init__ method LImpedance(double length, int nMacrosMin, int nBins)
  //-----------------------------------------------------

  static int LImpedance_init(pyORBIT_Object *self,
                             PyObject *args,
                             PyObject *kwds)
  {
    pyORBIT_Object* pyLImpedance = (pyORBIT_Object*) self;
    LImpedance* cpp_LImpedance = (LImpedance*) pyLImpedance->cpp_obj;

    double length = 1.0;
    int nMacrosMin;
    int nBins;

    if(!PyArg_ParseTuple(args, "dii:arguments",
                         &length, &nMacrosMin, &nBins))
    {
      ORBIT_MPI_Finalize("PyLImpedance - LImpedance(length, nMacrosMin, nBins) - constructor needs parameters.");
    }

    self->cpp_obj = new LImpedance(length, nMacrosMin, nBins);

    ((LImpedance*) self->cpp_obj)->setPyWrapper((PyObject*) self);

    return 0;
  }

  //-----------------------------------------------------
  // assignImpedance: Routine to import a python complex tuple
  // and convert to c++ impedance array
  //-----------------------------------------------------

  static PyObject* LImpedance_assignImpedance(PyObject *self, PyObject *args)
  {
    pyORBIT_Object* pyLImpedance = (pyORBIT_Object*) self;
    LImpedance* cpp_LImpedance = (LImpedance*) pyLImpedance->cpp_obj;
    PyObject* py_cmplx_arr;

    if(!PyArg_ParseTuple(args, "O:get_complex_arr", &py_cmplx_arr))
    {
      ORBIT_MPI_Finalize("ERROR! You have to specify a parameter - array of complex numbers!");
    }
    if(PySequence_Check(py_cmplx_arr) != 1)
    {
      ORBIT_MPI_Finalize("ERROR! You have to specify a parameter - array of complex numbers!");
    }

    int size = PySequence_Size(py_cmplx_arr);
    Py_complex cmplx;
    PyObject* py_cmplx;
    double real, imag;
    for(int n = 0; n < size; n++)
    {
      py_cmplx = PySequence_Fast_GET_ITEM(py_cmplx_arr, n);
      if(!PyComplex_Check(py_cmplx))
      {
        ORBIT_MPI_Finalize("ERROR! No complex numbers!");
      }
      cmplx = PyComplex_AsCComplex(py_cmplx);
      real = cmplx.real;
      imag = cmplx.imag;
      cpp_LImpedance->assignImpedanceValue(n, real, imag);
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //-----------------------------------------------------
  // assignImpedanceValue(int, real, real).
  // Wraps the LImpedance routine assigning an impedance mode
  //-----------------------------------------------------

  static PyObject* LImpedance_assignImpedanceValue(PyObject *self,
                                                   PyObject *args)
  {
    pyORBIT_Object* pyLImpedance = (pyORBIT_Object*) self;
    LImpedance* cpp_LImpedance = (LImpedance*) pyLImpedance->cpp_obj;

    int n = 0;
    double real = 0.0;
    double imag = 0.0;

    if(!PyArg_ParseTuple(args, "idd:arguments", &n, &real, &imag))
    {
      ORBIT_MPI_Finalize("PyLImpedance - assignImpedanceValue(n, real, imag) - constructor needs parameters.");
    }
    cpp_LImpedance->assignImpedanceValue(n, real, imag);

    Py_INCREF(Py_None);
    return Py_None;  
  }

  //-----------------------------------------------------
  //  trackBunchBunch(Bunch* bunch)
  //-----------------------------------------------------

  static PyObject* LImpedance_trackBunch(PyObject *self, PyObject *args)
  {
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pyLImpedance = (pyORBIT_Object*) self;
    LImpedance* cpp_LImpedance = (LImpedance*) pyLImpedance->cpp_obj;
    PyObject* pyBunch;

    if(!PyArg_ParseTuple(args, "O:trackBunch", &pyBunch))
    {
      ORBIT_MPI_Finalize("PyLImpedance.trackBunch(pyBunch) - method needs parameters.");
    }

    PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
    if(!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type))
    {
      ORBIT_MPI_Finalize("PyLImpedance.trackBunch(pyBunch) - pyBunch is not a Bunch.");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;		
    cpp_LImpedance->trackBunch(cpp_bunch);

    Py_INCREF(Py_None);
    return Py_None;  
  }

  //-----------------------------------------------------
  // Destructor for python LImpedance class (__del__ method).
  //-----------------------------------------------------

  static void LImpedance_del(pyORBIT_Object* self)
  {
    LImpedance* cpp_LImpedance = (LImpedance*) self->cpp_obj;
    if(cpp_LImpedance != NULL)
    {
      delete cpp_LImpedance;
    }
      self->ob_type->tp_free((PyObject*)self);
  }

  //-----------------------------------------------------
  // Definition of the methods of the python LImpedance wrapper class
  // They will be vailable from python level
  //-----------------------------------------------------

  static PyMethodDef LImpedanceClassMethods[] =
  {
    {"assignImpedance", LImpedance_assignImpedance, METH_VARARGS,
     "assigns overall impedance - assignImpedance(Z))"},
    {"assignImpedanceValue",  LImpedance_assignImpedanceValue, METH_VARARGS,
     "assigns impedance for the nth mode - assignImpedanceValue(n, real, imag)"},
    {"trackBunch", LImpedance_trackBunch, METH_VARARGS,
     "trackBunch tracks the bunch - trackBunch(pyBunch)"},
    {NULL}
  };

  //-----------------------------------------------------
  // Definition of the members of the python LImpedance wrapper class
  // They will be available from python level
  //-----------------------------------------------------

  static PyMemberDef LImpedanceClassMembers [] =
  {
    {NULL}
  };

  //-----------------------------------------------------
  //new python LImpedance wrapper type definition
  //-----------------------------------------------------

  static PyTypeObject pyORBIT_LImpedance_Type =
  {
    PyObject_HEAD_INIT(NULL)
    0, /*ob_size*/
    "LImpedance", /*tp_name*/
    sizeof(pyORBIT_Object), /*tp_basicsize*/
    0, /*tp_itemsize*/
    (destructor) LImpedance_del, /*tp_dealloc*/
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
    "The LImpedance python wrapper", /* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    LImpedanceClassMethods, /* tp_methods */
    LImpedanceClassMembers, /* tp_members */
    0, /* tp_getset */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    (initproc) LImpedance_init, /* tp_init */
    0, /* tp_alloc */
    LImpedance_new, /* tp_new */
  };

  //-----------------------------------------------------
  // Initialization function of the pyLImpedance class
  // It will be called from LImpedance wrapper initialization
  //-----------------------------------------------------

  void initLImpedance(PyObject* module)
  {
    if (PyType_Ready(&pyORBIT_LImpedance_Type) < 0) return;
    Py_INCREF(&pyORBIT_LImpedance_Type);
    PyModule_AddObject(module, "LImpedance",
                       (PyObject *)&pyORBIT_LImpedance_Type);
  }

#ifdef __cplusplus
}
#endif

//end of namespace wrap_impedances
}

