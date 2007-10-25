//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    pyORBIT_Object.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    06/08/2007
//
// DESCRIPTION
//    The wrapper around any pointer to the ORBIT C++ class instance.
//
///////////////////////////////////////////////////////////////////////////

#ifndef PY_ORBIT_OBJECT_H
#define PY_ORBIT_OBJECT_H

#include "structmember.h"

#ifdef __cplusplus
extern "C" {
#endif

 typedef struct {
   PyObject_HEAD
   void* cpp_obj;
 } pyORBIT_Object;

#ifdef __cplusplus
}
#endif

#endif
