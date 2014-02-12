#ifndef WRAP_SPACE_CHARGE_GRID_1D_H
#define WRAP_SPACE_CHARGE_GRID_1D_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_spacecharge
  {
    void initGrid1D(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
