#ifndef WRAP_MATRIX_GENERATOR_H
#define WRAP_MATRIX_GENERATOR_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_teapotbase_matrix_generator{
    void initMatrixGenerator(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
