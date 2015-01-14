#ifndef WRAP_IMPEDANCES_H
#define WRAP_IMPEDANCES_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_impedances
  {
    void initimpedances(void);
    PyObject* getImpedanceType(char* name);
  }

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_RFCAVITIES_H

