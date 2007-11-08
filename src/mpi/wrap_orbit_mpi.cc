#include "Python.h"
#include "orbit_mpi.hh"

#include <cstring>
#include <iostream>

#include "wrap_orbit_mpi.hh"

//wrappers of mpi objects
#include "wrap_mpi_comm.hh"

#ifdef __cplusplus
extern "C" {
#endif

 namespace wrap_orbit_mpi{

    static PyObject* get_cpu_name(PyObject *self, PyObject *args) {
      char* name = new char[MPI_MAX_PROCESSOR_NAME];
      int len;
      ORBIT_MPI_Get_processor_name(name,&len);
      PyObject* nm = PyString_FromStringAndSize(name, len);
      delete [] name;
      return nm;
    }


    static PyObject* get_comm_world_size(PyObject *self, PyObject *args) {
      int result = 0;
      ORBIT_MPI_Comm_size(MPI_COMM_WORLD,&result);
      return PyInt_FromLong((long)result);
    }


    static PyObject* get_comm_world_rank(PyObject *self, PyObject *args) {
      int result = 0;
      ORBIT_MPI_Comm_rank(MPI_COMM_WORLD,&result);
      return PyInt_FromLong((long)result);
    }

    static PyObject* get_time(PyObject *self, PyObject *args) {
      double result ;
      result = (double) ORBIT_MPI_Wtime();
      return PyFloat_FromDouble(result);
    }


  //Finalizes the execution of program
  //  the action is depended on the number of arguments
  //  () - no message
  //  (message) - will print message
  //this is a wrapper of
  // ORBIT_MPI_Finalize(const char* message)
  static PyObject* finalize(PyObject *self, PyObject *args){
    //if nVars == 0 no message
    //if nVars == 1 stop with message
    int nVars = PyTuple_Size(args);

    const char* message = NULL;

    if(nVars == 0 ||  nVars == 1){

      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!
        //NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(	args,"s:finalize",&message)){
					ORBIT_MPI_Finalize("orbit_mpi - something wrong with error message.");
        }

				ORBIT_MPI_Finalize(message);
      }

    }
    else{
			ORBIT_MPI_Finalize("orbit_mpi. You should call finalize() or finalize(message)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }


    //end of namespace
 }

  static PyMethodDef orbit_mpiMethods[] = {
    { (char *)"cpu",  wrap_orbit_mpi::get_cpu_name,  METH_VARARGS },
    { (char *)"world_size", wrap_orbit_mpi::get_comm_world_size, METH_VARARGS },
    { (char *)"world_rank", wrap_orbit_mpi::get_comm_world_rank, METH_VARARGS },
    { (char *)"time", wrap_orbit_mpi::get_time,      METH_VARARGS },
    { (char *)"finalize", wrap_orbit_mpi::finalize,      METH_VARARGS },
    { NULL, NULL }
  };


  void wrap_orbit_mpi::initorbit_mpi(void) {
    PyObject *m, *d;
    m = Py_InitModule((char*)"orbit_mpi",orbit_mpiMethods);
    d = PyModule_GetDict(m);
		
		//add MPI_Comm class and fields
		wrap_orbit_mpi_comm::init_orbit_mpi_comm(m);
  }


#ifdef __cplusplus
}
#endif
