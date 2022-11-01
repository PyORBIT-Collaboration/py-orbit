#ifndef WRAP_ORBIT_CONVEX_APERTURE_SHAPE_HH_
#define WRAP_ORBIT_CONVEX_APERTURE_SHAPE_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_convex_aperture_shape{
    void initConvexApertureShape(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_CONVEX_APERTURE_SHAPE_HH_*/