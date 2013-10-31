#ifndef WRAP_FREQUENCY_CAV_H
#define WRAP_FREQUENCY_CAV_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_rfcavities
  {
    void initFrequency_Cav(PyObject* module);
  }

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_FREQUENCY_CAV_H
