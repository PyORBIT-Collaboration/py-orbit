#ifndef WRAP_TIMPEDANCE_H
#define WRAP_TIMPEDANCE_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_impedances
  {
    void initTImpedance(PyObject* module);
  }

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // WRAP_TIMPEDANCE_H

