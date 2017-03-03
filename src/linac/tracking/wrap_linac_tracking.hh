#ifndef WRAP_LINAC_TRACKING_H
#define WRAP_LINAC_TRACKING_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

namespace wrap_linac_tracking
{
    void initlinactracking();
    PyObject* getLinacTrackingType(char* name);
}

#ifdef __cplusplus
}
#endif

#endif // WRAP_LINAC_TRACKING_H check

