#ifndef WRAP_ORBIT_APERTURE_HH_
#define WRAP_ORBIT_APERTURE_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_aperture{
    void initAperture(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_APERTURE_HH_*/
