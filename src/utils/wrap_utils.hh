#ifndef WRAP_ORBIT_UTILS_H
#define WRAP_ORBIT_UTILS_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

namespace wrap_orbit_utils{
	void initutils(void);
	PyObject* getOrbitUtilsType(const char* name);
}

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_ORBIT_UTILS_H
