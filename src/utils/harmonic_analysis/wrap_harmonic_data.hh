#ifndef WRAP_ORBIT_UTILS_HARMONIC_DATA_HH_
#define WRAP_ORBIT_UTILS_HARMONIC_DATA_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_harmonicdata{
    void initHarmonicData(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_UTILS_HARMONIC_DATA_HH_*/
