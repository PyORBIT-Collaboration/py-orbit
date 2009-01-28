#ifndef WRAP_SPACE_CHARGE_GRID_2D_H
#define WRAP_SPACE_CHARGE_GRID_2D_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_spacecharge{
    void initGrid2D(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
