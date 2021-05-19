#ifndef WRAP_DIPOLE_FIELD_SOURCE_H
#define WRAP_DIPOLE_FIELD_SOURCE_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_dipole_field_source{
    void initDipoleFieldSource(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
