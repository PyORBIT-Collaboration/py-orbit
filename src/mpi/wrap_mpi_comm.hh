#ifndef WRAP_ORBIT_MPI_COMM_H
#define WRAP_ORBIT_MPI_COMM_H

///////////////////////////////////////////////////////////////////////////
//
// This is a wrapper for the MPI_Comm data type from MPI
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_orbit_mpi_comm{
    void init_orbit_mpi_comm(PyObject* module);
	  
		//The function that will be exposed as C/C++ API for MPI_Comm		
		pyORBIT_MPI_Comm* newMPI_Comm();
		void freeMPI_Comm(pyORBIT_MPI_Comm* pyMPI_Comm);
		
		//Returns the type object for the specified name, like MPI_Comm
		PyObject* getMPI_CommType(char* name);
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
