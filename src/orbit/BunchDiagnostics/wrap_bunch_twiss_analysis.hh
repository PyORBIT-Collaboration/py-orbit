#ifndef WRAP_ORBIT_BUNCH_TWISS_ANALYSIS_HH_
#define WRAP_ORBIT_BUNCH_TWISS_ANALYSIS_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_bunch_twiss_analysis{
    void initbunchtwissanalysis(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_BUNCH_TWISS_ANALYSIS_HH_*/
