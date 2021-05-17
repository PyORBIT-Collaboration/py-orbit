#ifndef WRAP_FIELD_SOURCES_MODULE_H
#define WRAP_FIELD_SOURCES_MODULE_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_field_sources_module{
    void initFieldSourcesModule(PyObject* module, const char* part_with_id_module_name);
  }

#ifdef __cplusplus
}
#endif

#endif
