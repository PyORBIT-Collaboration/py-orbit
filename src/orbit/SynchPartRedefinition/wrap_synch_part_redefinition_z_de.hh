#ifndef WRAP_SYNCH_PARTICLE_REDEFINITION_HH_
#define WRAP_SYNCH_PARTICLE_REDEFINITION_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_synch_part_redefinition{
    void initsynchpartredefinition(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_SYNCH_PARTICLE_REDEFINITION_HH_*/
