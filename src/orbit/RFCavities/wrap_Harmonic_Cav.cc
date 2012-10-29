#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_Harmonic_Cav.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "Harmonic_Cav.hh"

using namespace OrbitUtils;

namespace wrap_rfcavities
{

#ifdef __cplusplus
extern "C"
{
#endif

//---------------------------------------------------------
// Python Harmonic_Cav class definition
//---------------------------------------------------------

//-----------------------------------------------------
// Constructor for python class wrapping Harmonic_Cav instance
// It never will be called directly
//-----------------------------------------------------

static PyObject* Harmonic_Cav_new(PyTypeObject *type,
                                  PyObject *args,
                                  PyObject *kwds)
{
  pyORBIT_Object* self;
  self = (pyORBIT_Object*) type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject*) self;
}

//-----------------------------------------------------
// Initializator for python Harmonic_Cav class
// This is implementation of the __init__ method
//-----------------------------------------------------

static int Harmonic_Cav_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
{
  double ZtoPhi    = 0.0;
  double dESync    = 0.0;
  double RFHNum    = 0.0;
  double RFVoltage = 0.0;
  double RFPhase   = 0.0;

  if(!PyArg_ParseTuple(args, "ddddd:arguments",
                       &ZtoPhi, &dESync, &RFHNum, &RFVoltage, &RFPhase))
  {
    ORBIT_MPI_Finalize("PyBunch - addParticle - cannot parse arguments! They should be (ZtoPhi, dESync, RFHNum, RFVoltage, RFPhase)");
  }
  self->cpp_obj = new Harmonic_Cav(ZtoPhi, dESync, RFHNum, RFVoltage, RFPhase);
  ((Harmonic_Cav*) self->cpp_obj)->setPyWrapper((PyObject*) self);
  return 0;
}

//-----------------------------------------------------
// Destructor for python Harmonic_Cav class (__del__ method)
//-----------------------------------------------------

static void Harmonic_Cav_del(pyORBIT_Object* self)
{
  Harmonic_Cav* cpp_Harmonic_Cav = (Harmonic_Cav*) self->cpp_obj;
  delete cpp_Harmonic_Cav;
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

static PyObject* Harmonic_Cav_ZtoPhi(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyHarmonic_Cav = (pyORBIT_Object*) self;
  Harmonic_Cav* cpp_Harmonic_Cav = (Harmonic_Cav*) pyHarmonic_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:ZtoPhi", &val))
    {
      ORBIT_MPI_Finalize("PyHarmonic_Cav_ZtoPhi(value) - value is needed");
    }
    cpp_Harmonic_Cav->setZtoPhi(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Harmonic_Cav->getZtoPhi();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyHarmonic_Cav_ZtoPhi. You should call ZtoPhi() or ZtoPhi(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Harmonic_Cav_dESync(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyHarmonic_Cav = (pyORBIT_Object*) self;
  Harmonic_Cav* cpp_Harmonic_Cav = (Harmonic_Cav*) pyHarmonic_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:dESync", &val))
    {
      ORBIT_MPI_Finalize("PyHarmonic_Cav_dESync(value) - value is needed");
    }
    cpp_Harmonic_Cav->setdESync(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Harmonic_Cav->getdESync();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyHarmonic_Cav_dESync. You should call dESync() or dESync(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Harmonic_Cav_RFHNum(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyHarmonic_Cav = (pyORBIT_Object*) self;
  Harmonic_Cav* cpp_Harmonic_Cav = (Harmonic_Cav*) pyHarmonic_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:RFHNum", &val))
    {
      ORBIT_MPI_Finalize("PyHarmonic_Cav_RFHNum(value) - value is needed");
    }
    cpp_Harmonic_Cav->setRFHNum(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Harmonic_Cav->getRFHNum();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyHarmonic_Cav_RFHNum. You should call RFHNum() or RFHNum(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Harmonic_Cav_RFVoltage(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyHarmonic_Cav = (pyORBIT_Object*) self;
  Harmonic_Cav* cpp_Harmonic_Cav = (Harmonic_Cav*) pyHarmonic_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:RFVoltage", &val))
    {
      ORBIT_MPI_Finalize("PyHarmonic_Cav_RFVoltage(value) - value is needed");
    }
    cpp_Harmonic_Cav->setRFVoltage(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Harmonic_Cav->getRFVoltage();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyHarmonic_Cav_RFVoltage. You should call RFVoltage() or RFVoltage(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Harmonic_Cav_RFPhase(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyHarmonic_Cav = (pyORBIT_Object*) self;
  Harmonic_Cav* cpp_Harmonic_Cav = (Harmonic_Cav*) pyHarmonic_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:RFPhase", &val))
    {
      ORBIT_MPI_Finalize("PyHarmonic_Cav_RFPhase(value) - value is needed");
    }
    cpp_Harmonic_Cav->setRFPhase(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Harmonic_Cav->getRFPhase();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyHarmonic_Cav_RFPhase. You should call RFPhase() or RFPhase(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

//-----------------------------------------------------
// trackBunch(Bunch* bunch)
//-----------------------------------------------------

static PyObject* Harmonic_Cav_trackBunch(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyHarmonic_Cav = (pyORBIT_Object*) self;
  Harmonic_Cav* cpp_Harmonic_Cav = (Harmonic_Cav*) pyHarmonic_Cav->cpp_obj;
  PyObject* pyBunch;
  if(!PyArg_ParseTuple(args, "O:trackBunch", &pyBunch))
  {
    ORBIT_MPI_Finalize("PyHarmonic_Cav - trackBunch(Bunch* bunch) - parameter is needed.");
  }
  PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if(!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type))
  {
    ORBIT_MPI_Finalize("PyHarmonic_Cav - trackBunch(Bunch* bunch) - the parameter should be a Bunch.");
  }
  Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*) pyBunch)->cpp_obj;
  cpp_Harmonic_Cav->trackBunch(cpp_bunch);
  Py_INCREF(Py_None);
  return Py_None;
}

//-----------------------------------------------------
// Definition of the methods of the python Harmonic_Cav wrapper class
// They will be available at python level
//-----------------------------------------------------

static PyMethodDef Harmonic_CavClassMethods[] =
{
  { "ZtoPhi", Harmonic_Cav_ZtoPhi ,METH_VARARGS,"Set ZtoPhi(value) or get ZtoPhi() the conversion from bunch longitudinal length coordinate to cavity phase"},
  { "dESync", Harmonic_Cav_dESync ,METH_VARARGS,"Set dESync(value) or get dESync() the energy change of the synchronous particle"},
  { "RFHNum", Harmonic_Cav_RFHNum ,METH_VARARGS,"Set RFHNum(value) or get RFHNum() the RF cavity harmonic number"},
  { "RFVoltage", Harmonic_Cav_RFVoltage ,METH_VARARGS,"Set RFVoltage(value) or get RFVoltage() the RF cavity voltage in GeV"},
  { "RFPhase", Harmonic_Cav_RFPhase ,METH_VARARGS,"Set RFPhase(value) or get RFPhase() the RF cavity phase in radians"},
  {"trackBunch", Harmonic_Cav_trackBunch, METH_VARARGS, "tracks the Bunch through a harmonic RF cavity"},
  {NULL}
};

//-----------------------------------------------------
// Definition of the members of the python Harmonic_Cav wrapper class
// They will be available at python level
//-----------------------------------------------------

static PyMemberDef Harmonic_CavClassMembers [] =
{
  {NULL}
};

//-----------------------------------------------------
//new python Harmonic_Cav wrapper type definition
//-----------------------------------------------------

static PyTypeObject pyORBIT_Harmonic_Cav_Type =
{
  PyObject_HEAD_INIT(NULL)
  0, /*ob_size*/
  "Harmonic_Cav", /*tp_name*/
  sizeof(pyORBIT_Object), /*tp_basicsize*/
  0, /*tp_itemsize*/
  (destructor) Harmonic_Cav_del , /*tp_dealloc*/
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
  "The Harmonic_Cav python wrapper", /* tp_doc */
  0, /* tp_traverse */
  0, /* tp_clear */
  0, /* tp_richcompare */
  0, /* tp_weaklistoffset */
  0, /* tp_iter */
  0, /* tp_iternext */
  Harmonic_CavClassMethods, /* tp_methods */
  Harmonic_CavClassMembers, /* tp_members */
  0, /* tp_getset */
  0, /* tp_base */
  0, /* tp_dict */
  0, /* tp_descr_get */
  0, /* tp_descr_set */
  0, /* tp_dictoffset */
  (initproc) Harmonic_Cav_init, /* tp_init */
  0, /* tp_alloc */
  Harmonic_Cav_new, /* tp_new */
};

//--------------------------------------------------
//Initialization function of the pyHarmonic_Cav class
//It will be called from Bunch wrapper initialization
//--------------------------------------------------

void initHarmonic_Cav(PyObject* module)
{
  if (PyType_Ready(&pyORBIT_Harmonic_Cav_Type) < 0) return;
  Py_INCREF(&pyORBIT_Harmonic_Cav_Type);
  PyModule_AddObject(module, "Harmonic_Cav", (PyObject*) &pyORBIT_Harmonic_Cav_Type);
}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_rfcavities
}
