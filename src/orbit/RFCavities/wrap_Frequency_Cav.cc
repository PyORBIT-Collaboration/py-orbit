#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_Frequency_Cav.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "Frequency_Cav.hh"

using namespace OrbitUtils;

namespace wrap_rfcavities
{

#ifdef __cplusplus
extern "C"
{
#endif

//---------------------------------------------------------
// Python Frequency_Cav class definition
//---------------------------------------------------------

//-----------------------------------------------------
// Constructor for python class wrapping Frequency_Cav instance
// It never will be called directly
//-----------------------------------------------------

static PyObject* Frequency_Cav_new(PyTypeObject *type,
                                   PyObject *args,
                                   PyObject *kwds)
{
  pyORBIT_Object* self;
  self = (pyORBIT_Object*) type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject*) self;
}

//-----------------------------------------------------
// Initialization for python Frequency_Cav class
// This is implementation of the __init__ method
//-----------------------------------------------------

static int Frequency_Cav_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
{
  double RFFreq  = 0.0;
  double RFE0TL  = 0.0;
  double RFPhase = 0.0;

  if(!PyArg_ParseTuple(args, "ddd:arguments",
                       &RFFreq, &RFE0TL, &RFPhase))
  {
    ORBIT_MPI_Finalize("PyBunch - addParticle - cannot parse arguments! They should be (RFFreq, RFE0TL, RFPhase)");
  }
  self->cpp_obj = new Frequency_Cav(RFFreq, RFE0TL, RFPhase);
  ((Frequency_Cav*) self->cpp_obj)->setPyWrapper((PyObject*) self);
  return 0;
}

//-----------------------------------------------------
// Destructor for python Frequency_Cav class (__del__ method)
//-----------------------------------------------------

static void Frequency_Cav_del(pyORBIT_Object* self)
{
  Frequency_Cav* cpp_Frequency_Cav = (Frequency_Cav*) self->cpp_obj;
  delete cpp_Frequency_Cav;
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

static PyObject* Frequency_Cav_RFFreq(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyFrequency_Cav = (pyORBIT_Object*) self;
  Frequency_Cav* cpp_Frequency_Cav = (Frequency_Cav*) pyFrequency_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:RFFreq", &val))
    {
      ORBIT_MPI_Finalize("PyFrequency_Cav_RFFreq(value) - value is needed");
    }
    cpp_Frequency_Cav->setRFFreq(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Frequency_Cav->getRFFreq();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyFrequency_Cav_RFFreq. You should call RFFreq() or RFFreq(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Frequency_Cav_RFE0TL(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyFrequency_Cav = (pyORBIT_Object*) self;
  Frequency_Cav* cpp_Frequency_Cav = (Frequency_Cav*) pyFrequency_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:RFE0TL", &val))
    {
      ORBIT_MPI_Finalize("PyFrequency_Cav_RFE0TL(value) - value is needed");
    }
    cpp_Frequency_Cav->setRFE0TL(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Frequency_Cav->getRFE0TL();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyFrequency_Cav_RFE0TL. You should call RFE0TL() or RFE0TL(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Frequency_Cav_RFPhase(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyFrequency_Cav = (pyORBIT_Object*) self;
  Frequency_Cav* cpp_Frequency_Cav = (Frequency_Cav*) pyFrequency_Cav->cpp_obj;
  int nVars = PyTuple_Size(args);
  double val = 0.;
  if(nVars == 1)
  {
    if(!PyArg_ParseTuple(args,"d:RFPhase", &val))
    {
      ORBIT_MPI_Finalize("PyFrequency_Cav_RFPhase(value) - value is needed");
    }
    cpp_Frequency_Cav->setRFPhase(val);
    return Py_BuildValue("d", val);
  }
  else if(nVars == 0)
  {
    val = cpp_Frequency_Cav->getRFPhase();
    return Py_BuildValue("d", val);
  }
  else
  {
    ORBIT_MPI_Finalize("PyFrequency_Cav_RFPhase. You should call RFPhase() or RFPhase(value)");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

//-----------------------------------------------------
// trackBunch(Bunch* bunch)
//-----------------------------------------------------

static PyObject* Frequency_Cav_trackBunch(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyFrequency_Cav = (pyORBIT_Object*) self;
  Frequency_Cav* cpp_Frequency_Cav = (Frequency_Cav*) pyFrequency_Cav->cpp_obj;
  PyObject* pyBunch;
  if(!PyArg_ParseTuple(args, "O:trackBunch", &pyBunch))
  {
    ORBIT_MPI_Finalize("PyFrequency_Cav - trackBunch(Bunch* bunch) - parameter is needed.");
  }
  PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if(!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type))
  {
    ORBIT_MPI_Finalize("PyFrequency_Cav - trackBunch(Bunch* bunch) - the parameter should be a Bunch.");
  }
  Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*) pyBunch)->cpp_obj;
  cpp_Frequency_Cav->trackBunch(cpp_bunch);
  Py_INCREF(Py_None);
  return Py_None;
}

//-----------------------------------------------------
// Definition of the methods of the python Frequency_Cav wrapper class
// They will be available at python level
//-----------------------------------------------------

static PyMethodDef Frequency_CavClassMethods[] =
{
  { "RFFreq", Frequency_Cav_RFFreq ,METH_VARARGS,"Set RFFreq(value) or get RFFreq() the RF cavity frequency in Hz"},
  { "RFE0TL", Frequency_Cav_RFE0TL ,METH_VARARGS,"Set RFE0TL(value) or get RFE0TL() the RF cavity E0TL in GeV"},
  { "RFPhase", Frequency_Cav_RFPhase ,METH_VARARGS,"Set RFPhase(value) or get RFPhase() the RF cavity phase in radians"},
  {"trackBunch", Frequency_Cav_trackBunch, METH_VARARGS, "tracks the Bunch through a frequency-specified RF cavity"},
  {NULL}
};

//-----------------------------------------------------
// Definition of the members of the python Frequency_Cav wrapper class
// They will be available at python level
//-----------------------------------------------------

static PyMemberDef Frequency_CavClassMembers [] =
{
  {NULL}
};

//-----------------------------------------------------
//new python Frequency_Cav wrapper type definition
//-----------------------------------------------------

static PyTypeObject pyORBIT_Frequency_Cav_Type =
{
  PyObject_HEAD_INIT(NULL)
  0, /*ob_size*/
  "Frequency_Cav", /*tp_name*/
  sizeof(pyORBIT_Object), /*tp_basicsize*/
  0, /*tp_itemsize*/
  (destructor) Frequency_Cav_del , /*tp_dealloc*/
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
  "The Frequency_Cav python wrapper", /* tp_doc */
  0, /* tp_traverse */
  0, /* tp_clear */
  0, /* tp_richcompare */
  0, /* tp_weaklistoffset */
  0, /* tp_iter */
  0, /* tp_iternext */
  Frequency_CavClassMethods, /* tp_methods */
  Frequency_CavClassMembers, /* tp_members */
  0, /* tp_getset */
  0, /* tp_base */
  0, /* tp_dict */
  0, /* tp_descr_get */
  0, /* tp_descr_set */
  0, /* tp_dictoffset */
  (initproc) Frequency_Cav_init, /* tp_init */
  0, /* tp_alloc */
  Frequency_Cav_new, /* tp_new */
};

//--------------------------------------------------
//Initialization function of the pyFrequency_Cav class
//It will be called from Bunch wrapper initialization
//--------------------------------------------------

void initFrequency_Cav(PyObject* module)
{
  if (PyType_Ready(&pyORBIT_Frequency_Cav_Type) < 0) return;
  Py_INCREF(&pyORBIT_Frequency_Cav_Type);
  PyModule_AddObject(module, "Frequency_Cav", (PyObject*) &pyORBIT_Frequency_Cav_Type);
}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_rfcavities
}
