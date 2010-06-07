#ifndef WRAP_BASE_RF_GAP_H
#define WRAP_BASE_RF_GAP_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_linac{
    void initBaseRfGap(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
