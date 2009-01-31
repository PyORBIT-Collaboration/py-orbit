#ifndef WRAP_UTILS_PHASE_VECTOR_H
#define WRAP_UTILS_PHASE_VECTOR_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_utils_phase_vector{
    void initPhaseVector(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
