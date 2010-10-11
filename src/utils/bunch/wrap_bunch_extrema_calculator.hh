#ifndef WRAP_UTILS_BUNCH_H
#define WRAP_UTILS_BUNCH_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_utils_bunch{
    void initBunchExtremaCalculator(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
