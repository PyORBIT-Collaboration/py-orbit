#ifndef WRAP_ORBIT_UTILS_FUNCTION_HH_
#define WRAP_ORBIT_UTILS_FUNCTION_HH_


#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_function{
    void initFunction(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_UTILS_FUNCTION_HH_*/
