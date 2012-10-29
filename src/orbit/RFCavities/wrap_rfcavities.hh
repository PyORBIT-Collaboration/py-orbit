#ifndef WRAP_RFCAVITIES_H
#define WRAP_RFCAVITIES_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_rfcavities
  {
    void initrfcavities(void);
    PyObject* getRFCavityType(char* name);
  }

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_RFCAVITIES_H

