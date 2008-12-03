#ifndef WRAP_CPP_EXTERNAL_EFFECTS_HH_
#define WRAP_CPP_EXTERNAL_EFFECTS_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_tracker3dfield_las_strip_external_effects{
    void initLasStripExternalEffects(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_CPP_EXTERNAL_EFFECTS_HH_*/



