#ifndef WRAP_PY_EXTERNAL_EFFECTS_H
#define WRAP_PY_EXTERNAL_EFFECTS_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_trackerrk4_py_external_effects{
    void initPyExternalEffects(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
