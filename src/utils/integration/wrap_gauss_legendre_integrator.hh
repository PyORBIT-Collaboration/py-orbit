#ifndef WRAP_ORBIT_UTILS_GAUSS_LEGENDRE_INTEGRATOR_HH_
#define WRAP_ORBIT_UTILS_GAUSS_LEGENDRE_INTEGRATOR_HH_


#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_gl_integrator{
    void initGLIntegrator(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_UTILS_GAUSS_LEGENDRE_INTEGRATOR_HH_*/
