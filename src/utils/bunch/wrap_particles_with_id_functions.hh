#ifndef WRAP_UTILS_BUNCH_FUNCTIONS_H
#define WRAP_UTILS_BUNCH_FUNCTIONS_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_utils_bunch_functions{
    void initParticlesWithIdFunctions(PyObject* module, const char* part_with_id_module_name);
  }

#ifdef __cplusplus
}
#endif

#endif
