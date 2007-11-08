#ifndef WRAP_ORBIT_MPI_GROUP_H
#define WRAP_ORBIT_MPI_GROUP_H

///////////////////////////////////////////////////////////////////////////
//
// This is a wrapper for the MPI_Group data type from MPI
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_orbit_mpi_group{
    void init_orbit_mpi_group(PyObject* module);
		
		//The function that will be exposed as C/C++ API for MPI_Group
		pyORBIT_MPI_Group* newMPI_Group();
		void freeMPI_Group(pyORBIT_MPI_Group* pyMPI_Group);
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
