#ifndef WRAP_ORBIT_UTILS_POLYNOMIAL_HH_
#define WRAP_ORBIT_UTILS_POLYNOMIAL_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_polynomial{
    void initPolynomial(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_UTILS_POLYNOMIAL_HH_*/
