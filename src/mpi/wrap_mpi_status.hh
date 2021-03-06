#ifndef WRAP_ORBIT_MPI_STATUS_H
#define WRAP_ORBIT_MPI_STATUS_H

///////////////////////////////////////////////////////////////////////////
//
// This is a wrapper for the MPI_Status data type from MPI
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_orbit_mpi_status{
    void init_orbit_mpi_status(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
