#ifndef WRAP_BARRIER_CAV_H
#define WRAP_BARRIER_CAV_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_rfcavities
  {
    void initBarrier_Cav(PyObject* module);
  }

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_BARRIER_CAV_H
