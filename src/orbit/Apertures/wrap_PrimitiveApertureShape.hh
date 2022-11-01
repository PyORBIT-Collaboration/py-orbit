#ifndef WRAP_ORBIT_PRIMITIVE_APERTURE_SHAPE_HH_
#define WRAP_ORBIT_PRIMITIVE_APERTURE_SHAPE_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_primitive_aperture_shape{
    void initPrimitiveApertureShape(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_PRIMITIVE_APERTURE_SHAPE_HH_*/