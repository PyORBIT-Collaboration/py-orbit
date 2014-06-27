#ifndef WRAP_SPACE_CHARGE_CALC_SLICE_BY_SLICE_2D_H
#define WRAP_SPACE_CHARGE_CALC_SLICE_BY_SLICE_2D_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_spacecharge{
    void initSpaceChargeCalcSliceBySlice2D(PyObject* module);
  }
	
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // WRAP_SPACE_CHARGE_CALC_SLICE_BY_SLICE_2D_H
