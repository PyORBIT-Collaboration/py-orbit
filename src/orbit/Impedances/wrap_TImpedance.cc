#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_TImpedance.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "TImpedance.hh"

using namespace OrbitUtils;

namespace wrap_impedances
{

#ifdef __cplusplus
extern "C"
{
#endif

  //---------------------------------------------------------
  // Python TImpedance class definition
  //---------------------------------------------------------

  //-----------------------------------------------------
  // Constructor for python class wrapping TImpedance instance
  // It never will be called directly
  //-----------------------------------------------------

  static PyObject* TImpedance_new(PyTypeObject *type,
                                  PyObject *args,
                                  PyObject *kwds)
  {
    pyORBIT_Object* self;
    self = (pyORBIT_Object *) type->tp_alloc(type, 0);
    self->cpp_obj = NULL;
    // std::cerr << "The TImpedance new has been called!" << std::endl;
    return (PyObject *) self;
  }

  //-----------------------------------------------------
  // Initialization of python TImpedance class
  // This is implementation of the
  // __init__ method TImpedance(double length,
  //                            int nMacrosMin,
  //                            int nBins,
  //                            int useX,
  //                            int useY)
  //-----------------------------------------------------

  static int TImpedance_init(pyORBIT_Object *self,
                             PyObject *args,
                             PyObject *kwds)
  {
    pyORBIT_Object* pyTImpedance = (pyORBIT_Object*) self;
    TImpedance* cpp_TImpedance = (TImpedance*) pyTImpedance->cpp_obj;

    double length  = 1.0;
    int nMacrosMin = 0;
    int nBins      = 0;
    int useX       = 0;
    int useY       = 0;

    if(!PyArg_ParseTuple(args, "diiii:arguments",
                         &length, &nMacrosMin, &nBins, &useX, &useY))
    {
      ORBIT_MPI_Finalize("PyTImpedance - TImpedance(length, nMacrosMin, nBins, useX, useY) - constructor needs parameters.");
    }

    self->cpp_obj = new TImpedance(length, nMacrosMin, nBins, useX, useY);

    ((TImpedance*) self->cpp_obj)->setPyWrapper((PyObject*) self);

    return 0;
  }

  //-----------------------------------------------------
  // assignLatFuncs(double qX, double alphaX, double betaX,
  //                double qY, double alphaY, double betaY)
  // Wraps the TImpedance routine assigning lattice functions
  //-----------------------------------------------------

  static PyObject* TImpedance_assignLatFuncs(PyObject *self,
                                             PyObject *args)
  {
    pyORBIT_Object* pyTImpedance = (pyORBIT_Object*) self;
    TImpedance* cpp_TImpedance = (TImpedance*) pyTImpedance->cpp_obj;

    double qX     = 0.0;
    double alphaX = 0.0;
    double betaX  = 0.0;
    double qY     = 0.0;
    double alphaY = 0.0;
    double betaY  = 0.0;

    if(!PyArg_ParseTuple(args, "dddddd:arguments",
                         &qX, &alphaX, &betaX, &qY, &alphaY, &betaY))
    {
      ORBIT_MPI_Finalize("PyTImpedance - assignLatFuncs(double, double, double, double, double, double) - constructor needs parameters.");
    }
    cpp_TImpedance->assignLatFuncs(qX, alphaX, betaX, qY, alphaY, betaY);

    Py_INCREF(Py_None);
    return Py_None;  
  }

  //-----------------------------------------------------
  // assignImpedance: Routine to import python complex tuples
  // and convert to c++ X or Y impedance arrays
  //-----------------------------------------------------

  static PyObject* assignImpedance(PyObject *self, PyObject *args)
  {
    pyORBIT_Object* pyTImpedance = (pyORBIT_Object*) self;
    TImpedance* cpp_TImpedance = (TImpedance*) pyTImpedance->cpp_obj;
    const char* XorY;
    PyObject* py_cmplx_arrp;
    PyObject* py_cmplx_arrm;

    if(!PyArg_ParseTuple(args, "sOO:get_impedance", &XorY,
                         &py_cmplx_arrp, &py_cmplx_arrm))
    {
      ORBIT_MPI_Finalize("ERROR! You have to specify a parameters - XorY and two arrays of complex numbers!");
    }
    if(PySequence_Check(py_cmplx_arrp) != 1)
    {
      ORBIT_MPI_Finalize("ERROR! You have to specify a parameter - py_cmplx_arrp array of complex numbers!");
    }
    if(PySequence_Check(py_cmplx_arrm) != 1)
    {
      ORBIT_MPI_Finalize("ERROR! You have to specify a parameter - py_cmplx_arrm array of complex numbers!");
    }

    int size = PySequence_Size(py_cmplx_arrp);
    if(size != PySequence_Size(py_cmplx_arrm))
    {
      ORBIT_MPI_Finalize("ERROR! Size of py_cmplx_arrp != size of py_cmplx_arrm!");
    }

    Py_complex cmplx;
    PyObject* py_cmplx;

    int nm;
    double realp, imagp, realm, imagm;
    for(int n = 0; n < size; n++)
    {
      py_cmplx = PySequence_Fast_GET_ITEM(py_cmplx_arrp, n);
      if(!PyComplex_Check(py_cmplx))
      {
        ORBIT_MPI_Finalize("ERROR! py_cmplx_arrp - No complex numbers!");
      }
      cmplx = PyComplex_AsCComplex(py_cmplx);
      realp = cmplx.real;
      imagp = cmplx.imag;

      py_cmplx = PySequence_Fast_GET_ITEM(py_cmplx_arrm, n);
      if(!PyComplex_Check(py_cmplx))
      {
        ORBIT_MPI_Finalize("ERROR! py_cmplx_arrm - No complex numbers!");
      }
      cmplx = PyComplex_AsCComplex(py_cmplx);
      realm = cmplx.real;
      imagm = cmplx.imag;

      nm = 2 * size - n;
      if(XorY == "X")
      {
        cpp_TImpedance->assignImpedanceX(n, realp, imagp, realm, imagm);
        if(nm < 2 * size)
        {
          cpp_TImpedance->assignImpedanceX(nm, -realm, imagm, -realp, imagp);
        }
      }
      else if(XorY == "Y")
      {
        cpp_TImpedance->assignImpedanceY(n, realp, imagp, realm, imagm);
        if(nm < 2 * size)
        {
          cpp_TImpedance->assignImpedanceY(nm, -realm, imagm, -realp, imagp);
        }
      }
      else
      {
        ORBIT_MPI_Finalize("ERROR! XorY wrong value!");
      }
    }
    if(XorY == "X")
    {
      cpp_TImpedance->assignImpedanceX(size, 0.0, 0.0, 0.0, 0.0);
    }
    if(XorY == "Y")
    {
      cpp_TImpedance->assignImpedanceY(size, 0.0, 0.0, 0.0, 0.0);
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //-----------------------------------------------------
  // assignImpedanceX(int n,
  //                  double realp, double imagp,
  //                  double realm, double imagm)
  // Wraps the TImpedance routine assigning horizontal impedance mode
  //-----------------------------------------------------

  static PyObject* TImpedance_assignImpedanceX(PyObject *self,
                                               PyObject *args)
  {
    pyORBIT_Object* pyTImpedance = (pyORBIT_Object*) self;
    TImpedance* cpp_TImpedance = (TImpedance*) pyTImpedance->cpp_obj;

    int n = 0;
    double realp = 0.0;
    double imagp = 0.0;
    double realm = 0.0;
    double imagm = 0.0;

    if(!PyArg_ParseTuple(args, "idddd:arguments",
                         &n, &realp, &imagp, &realm, &imagm))
    {
      ORBIT_MPI_Finalize("PyTImpedance - assignImpedanceX(n, realp, imagp, realm, imagm) - constructor needs parameters.");
    }
    cpp_TImpedance->assignImpedanceX(n, realp, imagp, realm, imagm);

    Py_INCREF(Py_None);
    return Py_None;  
  }

  //-----------------------------------------------------
  // assignImpedanceY(int n,
  //                  double realp, double imagp,
  //                  double realm, double imagm)
  // Wraps the TImpedance routine assigning vertical impedance mode
  //-----------------------------------------------------

  static PyObject* TImpedance_assignImpedanceY(PyObject *self,
                                               PyObject *args)
  {
    pyORBIT_Object* pyTImpedance = (pyORBIT_Object*) self;
    TImpedance* cpp_TImpedance = (TImpedance*) pyTImpedance->cpp_obj;

    int n = 0;
    double realp = 0.0;
    double imagp = 0.0;
    double realm = 0.0;
    double imagm = 0.0;

    if(!PyArg_ParseTuple(args, "idddd:arguments",
                         &n, &realp, &imagp, &realm, &imagm))
    {
      ORBIT_MPI_Finalize("PyTImpedance - assignImpedanceY(n, realp, imagp, realm, imagm) - constructor needs parameters.");
    }
    cpp_TImpedance->assignImpedanceY(n, realp, imagp, realm, imagm);

    Py_INCREF(Py_None);
    return Py_None;  
  }

  //-----------------------------------------------------
  //  trackBunchBunch(Bunch* bunch)
  //-----------------------------------------------------

  static PyObject* TImpedance_trackBunch(PyObject *self, PyObject *args)
  {
    int nVars = PyTuple_Size(args);
    pyORBIT_Object* pyTImpedance = (pyORBIT_Object*) self;
    TImpedance* cpp_TImpedance = (TImpedance*) pyTImpedance->cpp_obj;
    PyObject* pyBunch;

    if(!PyArg_ParseTuple(args, "O:trackBunch", &pyBunch))
    {
      ORBIT_MPI_Finalize("PyTImpedance.trackBunch(pyBunch) - method needs parameters.");
    }

    PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
    if(!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type))
    {
      ORBIT_MPI_Finalize("PyTImpedance.trackBunch(pyBunch) - pyBunch is not a Bunch.");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;		
    cpp_TImpedance->trackBunch(cpp_bunch);

    Py_INCREF(Py_None);
    return Py_None;  
  }

  //-----------------------------------------------------
  // Destructor for python TImpedance class (__del__ method).
  //-----------------------------------------------------

  static void TImpedance_del(pyORBIT_Object* self)
  {
    TImpedance* cpp_TImpedance = (TImpedance*) self->cpp_obj;
    if(cpp_TImpedance != NULL)
    {
      delete cpp_TImpedance;
    }
      self->ob_type->tp_free((PyObject*)self);
  }

  //-----------------------------------------------------
  // Definition of the methods of the python TImpedance wrapper class
  // They will be vailable from python level
  //-----------------------------------------------------

  static PyMethodDef TImpedanceClassMethods[] =
  {
    {"assignImpedance", assignImpedance, METH_VARARGS,
     "assigns overall impedance - assignImpedance(XorY, Zp, Zm))"},
    {"assignImpedanceX",  TImpedance_assignImpedanceX, METH_VARARGS,
     "assigns horizontal impedance for the nth mode - assignImpedanceX(n, realp, imagp, realm, imagm)"},
    {"assignImpedanceY",  TImpedance_assignImpedanceY, METH_VARARGS,
     "assigns vertical impedance for the nth mode - assignImpedanceY(n, realp, imagp, realm, imagm)"},
    {"trackBunch", TImpedance_trackBunch, METH_VARARGS,
     "trackBunch tracks the bunch - trackBunch(pyBunch)"},
    {NULL}
  };

  //-----------------------------------------------------
  // Definition of the members of the python TImpedance wrapper class
  // They will be available from python level
  //-----------------------------------------------------

  static PyMemberDef TImpedanceClassMembers [] =
  {
    {NULL}
  };

  //-----------------------------------------------------
  //new python TImpedance wrapper type definition
  //-----------------------------------------------------

  static PyTypeObject pyORBIT_TImpedance_Type =
  {
    PyObject_HEAD_INIT(NULL)
    0, /*ob_size*/
    "TImpedance", /*tp_name*/
    sizeof(pyORBIT_Object), /*tp_basicsize*/
    0, /*tp_itemsize*/
    (destructor) TImpedance_del, /*tp_dealloc*/
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
    "The TImpedance python wrapper", /* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    TImpedanceClassMethods, /* tp_methods */
    TImpedanceClassMembers, /* tp_members */
    0, /* tp_getset */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    (initproc) TImpedance_init, /* tp_init */
    0, /* tp_alloc */
    TImpedance_new, /* tp_new */
  };

  //-----------------------------------------------------
  // Initialization function of the pyTImpedance class
  // It will be called from TImpedance wrapper initialization
  //-----------------------------------------------------

  void initTImpedance(PyObject* module)
  {
    if (PyType_Ready(&pyORBIT_TImpedance_Type) < 0) return;
    Py_INCREF(&pyORBIT_TImpedance_Type);
    PyModule_AddObject(module, "TImpedance",
                       (PyObject *)&pyORBIT_TImpedance_Type);
  }

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}

