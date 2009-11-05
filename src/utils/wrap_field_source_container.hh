#ifndef WRAP_FIELD_SOURCE_CONTAINER_HH_
#define WRAP_FIELD_SOURCE_CONTAINER_HH_


#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_field_source_container{
    void initFieldSourceContainer(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_FIELD_SOURCE_CONTAINER_HH_*/
