#ifndef WRAP_LIMPEDANCE_H
#define WRAP_LIMPEDANCE_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_impedances
  {
    void initLImpedance(PyObject* module);
  }

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // WRAP_LIMPEDANCE_H

