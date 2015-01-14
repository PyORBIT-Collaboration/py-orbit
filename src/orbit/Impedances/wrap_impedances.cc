#include "orbit_mpi.hh"

#include "wrap_LImpedance.hh"
#include "wrap_TImpedance.hh"

static PyMethodDef impedancesMethods[] = {{NULL,NULL}};


namespace wrap_impedances
{

#ifdef __cplusplus
extern "C"
{
#endif

void initimpedances()
{
  //create new module
  PyObject* module = Py_InitModule("impedances", impedancesMethods);
  wrap_impedances::initLImpedance(module);
  wrap_impedances::initTImpedance(module);
}

PyObject* getImpedanceType(char* name)
{
  PyObject* mod = PyImport_ImportModule("impedances");
  PyObject* pyType = PyObject_GetAttrString(mod, name);
  Py_DECREF(mod);
  Py_DECREF(pyType);
  return pyType;
}

#ifdef __cplusplus
}
#endif // __cplusplus

} // end of namespace wrap_impedances

