#include "orbit_mpi.hh"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ctime>

int ORBIT_MPI_Init(int *len, char ***ch){
  int res = 0;

#ifdef USE_MPI
  res = MPI_Init(len,ch);
#else
  res  = 1;
#endif

  return res;
}


int ORBIT_MPI_Initialized(int *init){
  int res = 0;

#ifdef USE_MPI
  res = MPI_Initialized(init);
#else
  init = 0;
#endif

  return res;
}


int ORBIT_MPI_Finalize(void){
  return ORBIT_MPI_Finalize(NULL);
}

int ORBIT_MPI_Finalize(const char* message){
  int res = 0;

  int rank;
  ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int init;
  ORBIT_MPI_Initialized(&init);

#ifdef USE_MPI
  if(init > 0){
    res = MPI_Finalize();
  }
#else
  res  = 1;
#endif

  if(rank == 0){
    if(Py_IsInitialized()){
      PyErr_SetString(PyExc_RuntimeError,"ORBIT_MPI_Finalize.");
      PyErr_Print();
      PyRun_SimpleString("import traceback; traceback.print_stack()");
    }
    if(message != NULL){
      std::cerr<<"Error PyORBIT Message:"<<std::endl;
      std::cerr<<message<<std::endl;
      std::cerr<<"Stop."<<std::endl;
    }
  }

  if(Py_IsInitialized()){
    Py_Exit(1);
  }
  else{
    exit(1);
  }
  return res;
}

int ORBIT_MPI_Get_processor_name(char *name, int* len){
  int res = 0;

#ifdef USE_MPI
  res = MPI_Get_processor_name(name,len);
#else
  res = 1;
  *len = 0;
  name = new char[0];
#endif

  return res;
}


int ORBIT_MPI_Comm_size(MPI_Comm comm, int * size){
  int res = 0;

#ifdef USE_MPI
  res = MPI_Comm_size(comm,size);
#else
  res  = 1;
  *size = 1;
#endif

  return res;
}

int ORBIT_MPI_Comm_rank(MPI_Comm comm, int * rank){
  int res = 0;

#ifdef USE_MPI
  res = MPI_Comm_rank(comm,rank);
#else
  res  = 1;
  *rank = 0;
#endif

  return res;
}

double ORBIT_MPI_Wtime(void){
  double time = 0;

#ifdef USE_MPI
  time = MPI_Wtime();
#else
  time  = ((double) clock())/CLOCKS_PER_SEC;
#endif

  return time;
}

int ORBIT_MPI_Allreduce(void* ar1, void* ar2, int n, MPI_Datatype data, MPI_Op op, MPI_Comm comm){
  int res = 0;

#ifdef USE_MPI
  res = MPI_Allreduce(ar1, ar2, n, data, op, comm);
#else
  res  = 1;
#endif

  return res;
}

int ORBIT_MPI_Send(void* ar, int n1, MPI_Datatype data, int n2, int n3, MPI_Comm comm){
  int res = 0;

#ifdef USE_MPI
  res = MPI_Send(ar, n1, data, n2, n3, comm);
#else
  res  = 1;
#endif

  return res;
}

int ORBIT_MPI_Recv(void* ar, int n1, MPI_Datatype data, int n2, int n3, MPI_Comm comm, MPI_Status * stat){
  int res = 0;

#ifdef USE_MPI
  res = MPI_Recv(ar, n1, data, n2, n3, comm, stat);
#else
  res  = 1;
#endif

  return res;
}

int ORBIT_MPI_Bcast(void* ar, int n1, MPI_Datatype data, int n2, MPI_Comm comm ){
  int res = 0;

#ifdef USE_MPI
  res = MPI_Bcast(ar, n1, data, n2, comm ) ;
#else
  res  = 1;
#endif

  return res;
}
