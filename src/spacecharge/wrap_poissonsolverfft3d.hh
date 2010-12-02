#ifndef WRAP_SPACE_CHARGE_POISSON_SOLVER_FFT_3D_H
#define WRAP_SPACE_CHARGE_POISSON_SOLVER_FFT_3D_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_spacecharge{
    void initPoissonSolverFFT3D(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
