#ifndef WRAP_QUAD_FIELD_SOURCE_H
#define WRAP_QUAD_FIELD_SOURCE_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_quad_field_source{
    void initQuadFieldSource(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
