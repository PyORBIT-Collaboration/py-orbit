#ifndef WRAP_RF_GAP_TTF_SLOW_H
#define WRAP_RF_GAP_TTF_SLOW_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_linac{
    void initRfGapTTF_slow(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif