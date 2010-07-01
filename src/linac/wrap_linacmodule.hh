#ifndef WRAP_LINAC_H
#define WRAP_LINAC_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif
  namespace wrap_linac{
	  void initlinac(void);
	  PyObject* getLinacType(char* name);
	}
	
#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_LINAC_H
