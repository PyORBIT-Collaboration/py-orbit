#ifndef WRAP_HARMONIC_CAV_H
#define WRAP_HARMONIC_CAV_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_rfcavities
  {
    void initHarmonic_Cav(PyObject* module);
  }

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_HARMONIC_CAV_H
