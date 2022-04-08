#ifndef WRAP_LOOP_FIELD_SOURCE_H
#define WRAP_LOOP_FIELD_SOURCE_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_loop_field_source{
    void initLoopFieldSource(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif