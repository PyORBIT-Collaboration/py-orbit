#ifndef WRAP_SUPER_FISH_RF_FIELD_SOURCE_H
#define WRAP_SUPER_FISH_RF_FIELD_SOURCE_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_linac{
    void initSuperFishFieldSource(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
