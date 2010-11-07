#ifndef WRAP_SPACE_CHARGE_GRID_3D_H
#define WRAP_SPACE_CHARGE_GRID_3D_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_spacecharge{
    void initGrid3D(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
