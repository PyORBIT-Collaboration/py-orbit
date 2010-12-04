#ifndef WRAP_SPACE_CHARGE_CALC_3D_H
#define WRAP_SPACE_CHARGE_CALC_3D_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_spacecharge{
    void initSpaceChargeCalc3D(PyObject* module);
  }
	
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // WRAP_SPACE_CHARGE_CALC_3D_H
