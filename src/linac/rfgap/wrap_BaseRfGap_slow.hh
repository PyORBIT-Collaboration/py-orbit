#ifndef WRAP_BASE_RF_GAP_SLOW_H
#define WRAP_BASE_RF_GAP_SLOW_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_linac{
    void initBaseRfGap_slow(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
