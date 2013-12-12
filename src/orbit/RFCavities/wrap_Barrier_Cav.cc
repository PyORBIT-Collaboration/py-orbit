#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_Barrier_Cav.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "Barrier_Cav.hh"

using namespace OrbitUtils;

namespace wrap_rfcavities
{

#ifdef __cplusplus
extern "C"
{
#endif

//---------------------------------------------------------
// Python Barrier_Cav class definition
//---------------------------------------------------------

//-----------------------------------------------------
// Constructor for python class wrapping Barrier_Cav instance
// It never will be called directly
//-----------------------------------------------------

static PyObject* Barrier_Cav_new(PyTypeObject *type,
                                 PyObject *args,
                                 PyObject *kwds)
{
  pyORBIT_Object* self;
  self = (pyORBIT_Object*) type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject*) self;
}

//-----------------------------------------------------
// Initialization for python Barrier_Cav class
// This is implementation of the __init__ method
//-----------------------------------------------------

static int Barrier_Cav_init(pyORBIT_Object *self,
                            PyObject *args,
                            PyObject *kwds)
{
  double ZtoPhi    = 0.0;
  double RFVoltage = 0.0;
  double RFPhasep  = 0.0;
  double RFPhasem  = 0.0;
  double dRFPhasep = 0.0;
  double dRFPhasem = 0.0;

  if(!PyArg_ParseTuple(args, "dddddd:arguments",
                       &ZtoPhi,
                       &RFVoltage,
                       &RFPhasep,
                       &RFPhasem,
                       &dRFPhasep,
                       &dRFPhasem))
  {
    ORBIT_MPI_Finalize("PyBarrier_Cav - Barrier_Cav_init - cannot parse arguments! They should be (ZtoPhi, RFVoltage, RFPhasep, RFPhasem, dRFPhasep, dRFPhasem)");
  }
  self->cpp_obj = new Barrier_Cav(ZtoPhi,
                                  RFVoltage,
                                  RFPhasep,
                                  RFPhasem,
                                  dRFPhasep,
                                  dRFPhasem);
  ((Barrier_Cav*) self->cpp_obj)->setPyWrapper((PyObject*) self);
  return 0;
}

//-----------------------------------------------------
// Destructor for python Barrier_Cav class (__del__ method)
//-----------------------------------------------------

static void Barrier_Cav_del(pyORBIT_Object* self)
{
  Barrier_Cav* cpp_Barrier_Cav = (Barrier_Cav*) self->cpp_obj;
  delete cpp_Barrier_Cav;
  self->ob_type->tp_free((PyObject*) self);
}

//-----------------------------------------------------
// Accessor routines:
// Sets or returns value depending on the number of arguments
// if nVars == 1 set value
// if nVars == 0 get value:
// variable(value) - sets new value
// variable()      - returns value
// These are implementations of the set(value) and get() methods
//-----------------------------------------------------

static PyObject* Barrier_Cav_ZtoPhi(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyBarrier_Cav = (pyORBIT_Object*) self;
  Barrier_Cav* cpp_Barrier_Cav = (Barrier_Cav*) pyBarrier_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:ZtoPhi", &val))
    {
      ORBIT_MPI_Finalize("PyBarrier_Cav_ZtoPhi(value) - value is needed");
    }
    cpp_Barrier_Cav->setZtoPhi(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Barrier_Cav->getZtoPhi();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyBarrier_Cav_ZtoPhi. You should call ZtoPhi() or ZtoPhi(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Barrier_Cav_RFVoltage(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyBarrier_Cav = (pyORBIT_Object*) self;
  Barrier_Cav* cpp_Barrier_Cav = (Barrier_Cav*) pyBarrier_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:RFVoltage", &val))
    {
      ORBIT_MPI_Finalize("PyBarrier_Cav_RFVoltage(value) - value is needed");
    }
    cpp_Barrier_Cav->setRFVoltage(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Barrier_Cav->getRFVoltage();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyBarrier_Cav_RFVoltage. You should call RFVoltage() or RFVoltage(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Barrier_Cav_RFPhasep(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyBarrier_Cav = (pyORBIT_Object*) self;
  Barrier_Cav* cpp_Barrier_Cav = (Barrier_Cav*) pyBarrier_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:RFPhasep", &val))
    {
      ORBIT_MPI_Finalize("PyBarrier_Cav_RFPhasep(value) - value is needed");
    }
    cpp_Barrier_Cav->setRFPhasep(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Barrier_Cav->getRFPhasep();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyBarrier_Cav_RFPhasep. You should call RFPhasep() or RFPhasep(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Barrier_Cav_RFPhasem(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyBarrier_Cav = (pyORBIT_Object*) self;
  Barrier_Cav* cpp_Barrier_Cav = (Barrier_Cav*) pyBarrier_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:RFPhasem", &val))
    {
      ORBIT_MPI_Finalize("PyBarrier_Cav_RFPhasem(value) - value is needed");
    }
    cpp_Barrier_Cav->setRFPhasem(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Barrier_Cav->getRFPhasem();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyBarrier_Cav_RFPhasem. You should call RFPhasem() or RFPhasem(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Barrier_Cav_dRFPhasep(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyBarrier_Cav = (pyORBIT_Object*) self;
  Barrier_Cav* cpp_Barrier_Cav = (Barrier_Cav*) pyBarrier_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:dRFPhasep", &val))
    {
      ORBIT_MPI_Finalize("PyBarrier_Cav_dRFPhasep(value) - value is needed");
    }
    cpp_Barrier_Cav->setdRFPhasep(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Barrier_Cav->getdRFPhasep();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyBarrier_Cav_dRFPhasep. You should call dRFPhasep() or dRFPhasep(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Barrier_Cav_dRFPhasem(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyBarrier_Cav = (pyORBIT_Object*) self;
  Barrier_Cav* cpp_Barrier_Cav = (Barrier_Cav*) pyBarrier_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:dRFPhasem", &val))
    {
      ORBIT_MPI_Finalize("PyBarrier_Cav_dRFPhasem(value) - value is needed");
    }
    cpp_Barrier_Cav->setdRFPhasem(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Barrier_Cav->getdRFPhasem();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyBarrier_Cav_dRFPhasem. You should call dRFPhasem() or dRFPhasem(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

//-----------------------------------------------------
// trackBunch(Bunch* bunch)
//-----------------------------------------------------

static PyObject* Barrier_Cav_trackBunch(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyBarrier_Cav = (pyORBIT_Object*) self;
  Barrier_Cav* cpp_Barrier_Cav = (Barrier_Cav*) pyBarrier_Cav->cpp_obj;
  PyObject* pyBunch;
  if(!PyArg_ParseTuple(args, "O:trackBunch", &pyBunch))
  {
    ORBIT_MPI_Finalize("PyBarrier_Cav - trackBunch(Bunch* bunch) - parameter is needed.");
  }
  PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if(!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type))
  {
    ORBIT_MPI_Finalize("PyBarrier_Cav - trackBunch(Bunch* bunch) - the parameter should be a Bunch.");
  }
  Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*) pyBunch)->cpp_obj;
  cpp_Barrier_Cav->trackBunch(cpp_bunch);
  Py_INCREF(Py_None);
  return Py_None;
}

//-----------------------------------------------------
// Definition of the methods of the python Barrier_Cav wrapper class
// They will be available at python level
//-----------------------------------------------------

static PyMethodDef Barrier_CavClassMethods[] =
{
  { "ZtoPhi", Barrier_Cav_ZtoPhi, METH_VARARGS, "Set ZtoPhi(value) or get ZtoPhi() the conversion from bunch longitudinal length coordinate to cavity phase"},
  { "RFVoltage", Barrier_Cav_RFVoltage, METH_VARARGS, "Set RFVoltage(value) or get RFVoltage() the RF cavity voltage in GeV"},
  { "RFPhasep", Barrier_Cav_RFPhasep, METH_VARARGS, "Set RFPhasep(value) or get RFPhasep() the RF cavity phase in radians"},
  { "RFPhasem", Barrier_Cav_RFPhasem, METH_VARARGS, "Set RFPhasem(value) or get RFPhasem() the RF cavity phase in radians"},
  { "dRFPhasep", Barrier_Cav_dRFPhasep, METH_VARARGS, "Set dRFPhasep(value) or get dRFPhasep() the RF cavity phase in radians"},
  { "dRFPhasem", Barrier_Cav_dRFPhasem, METH_VARARGS, "Set dRFPhasem(value) or get dRFPhasem() the RF cavity phase in radians"},
  {"trackBunch", Barrier_Cav_trackBunch, METH_VARARGS,"tracks the Bunch through a barrier RF cavity"},
  {NULL}
};

//-----------------------------------------------------
// Definition of the members of the python Barrier_Cav wrapper class
// They will be available at python level
//-----------------------------------------------------

static PyMemberDef Barrier_CavClassMembers [] =
{
  {NULL}
};

//-----------------------------------------------------
//new python Barrier_Cav wrapper type definition
//-----------------------------------------------------

static PyTypeObject pyORBIT_Barrier_Cav_Type =
{
  PyObject_HEAD_INIT(NULL)
  0, /*ob_size*/
  "Barrier_Cav", /*tp_name*/
  sizeof(pyORBIT_Object), /*tp_basicsize*/
  0, /*tp_itemsize*/
  (destructor) Barrier_Cav_del , /*tp_dealloc*/
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
  "The Barrier_Cav python wrapper", /* tp_doc */
  0, /* tp_traverse */
  0, /* tp_clear */
  0, /* tp_richcompare */
  0, /* tp_weaklistoffset */
  0, /* tp_iter */
  0, /* tp_iternext */
  Barrier_CavClassMethods, /* tp_methods */
  Barrier_CavClassMembers, /* tp_members */
  0, /* tp_getset */
  0, /* tp_base */
  0, /* tp_dict */
  0, /* tp_descr_get */
  0, /* tp_descr_set */
  0, /* tp_dictoffset */
  (initproc) Barrier_Cav_init, /* tp_init */
  0, /* tp_alloc */
  Barrier_Cav_new, /* tp_new */
};

//--------------------------------------------------
//Initialization function of the pyBarrier_Cav class
//It will be called from Bunch wrapper initialization
//--------------------------------------------------

void initBarrier_Cav(PyObject* module)
{
  if (PyType_Ready(&pyORBIT_Barrier_Cav_Type) < 0) return;
  Py_INCREF(&pyORBIT_Barrier_Cav_Type);
  PyModule_AddObject(module, "Barrier_Cav", (PyObject*) &pyORBIT_Barrier_Cav_Type);
}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_rfcavities
}
