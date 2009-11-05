#include "orbit_mpi.hh"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ctime>


/** A C wrapper around MPI_Init. */
int ORBIT_MPI_Init(int *len, char ***ch){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Init(len,ch);
#else
  res  = 1;
#endif
  return res;
}

/** A C wrapper around MPI_Initialized. */
int ORBIT_MPI_Initialized(int *init){
 int res = 0;
#if USE_MPI > 0
  res = MPI_Initialized(init);
#else
  *init = 0;
#endif
  return res;
}

/** A C wrapper around MPI_Finalize. */
int ORBIT_MPI_Finalize(){
  return ORBIT_MPI_Finalize(NULL);
}

/** A C wrapper around MPI_Finalize(message). */
int ORBIT_MPI_Finalize(const char* message){
  int res = 0;
  int rank;
  ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int init;
  ORBIT_MPI_Initialized(&init);
#if USE_MPI > 0
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

/** A C wrapper around MPI_Get_processor_name. */
int ORBIT_MPI_Get_processor_name(char *name, int* len){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Get_processor_name(name,len);
#else
  res = 1;
  *len = 0;
  name = "no mpi";
#endif

  return res;
}

/** A C wrapper around MPI_Wtime. */
double ORBIT_MPI_Wtime(void){
  double time = 0;
#if USE_MPI > 0
  time = MPI_Wtime();
#else
  time  = ((double) clock())/CLOCKS_PER_SEC;
#endif
  return time;
}

/** A C wrapper around MPI_Wtick. */
double ORBIT_MPI_Wtick(){  
	double tick = 0;
#if USE_MPI > 0
  tick = MPI_Wtick();
#else
  tick  = 1.0/CLOCKS_PER_SEC;
#endif
  return tick;
}

//--------------------------------------------------------
// MPI functions related to the MPI_Comm manipulations
//--------------------------------------------------------

/** A C wrapper around MPI_Comm_create. */
int ORBIT_MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *comm_out){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_create(comm, group, comm_out);
#else
  res = 1;
#endif

  return res;	
}

/** A C wrapper around . */
int ORBIT_MPI_Comm_group(MPI_Comm comm, MPI_Group *group ){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_group(comm, group);
#else
  res  = 1;
#endif		
	return res;	
}

/** A C wrapper around . */
int ORBIT_MPI_Comm_dup(MPI_Comm comm, MPI_Comm *comm_out){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_dup(comm, comm_out);
#else
  res  = 1;
#endif		
	return res;	
}

/** A C wrapper around . */
int ORBIT_MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *comm_out){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_split(comm, color, key, comm_out);
#else
  res  = 1;
#endif		
	return res;	
}

/** A C wrapper around . */
int ORBIT_MPI_Comm_remote_size(MPI_Comm comm, int *size){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_remote_size(comm, size);
#else
  res  = 1;
	*size = 0;
#endif		
	return res;	
}

/** A C wrapper around . */
int ORBIT_MPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_remote_group(comm, group);
#else
  res  = 1;
#endif		
	return res;
}

/** A C wrapper around MPI_Comm_test_inter. */
int ORBIT_MPI_Comm_test_inter(MPI_Comm comm, int *flag){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_test_inter(comm, flag);
#else
  res  = 1;
	*flag = 0;
#endif		
	return res;
}

/** A C wrapper around MPI_Comm_compare. */
int ORBIT_MPI_Comm_compare(MPI_Comm  comm1, MPI_Comm  comm2, int *result){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_compare(comm1, comm2, result);
#else
  res  = 1;
#endif		
	return res;	
}

/** A C wrapper around MPI_Comm_set_name. */
int ORBIT_MPI_Comm_set_name(MPI_Comm com, char *name){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_set_name(com, name);
#else
  res  = 1;
#endif		
	return res;	
}

/** A C wrapper around MPI_Comm_get_name. */
int ORBIT_MPI_Comm_get_name(MPI_Comm comm, char *namep, int *reslen){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_get_name(comm, namep, reslen);
#else
  res  = 1;
	namep = "no mpi";
	*reslen = 6;
#endif		
	return res;	
}

/** A C wrapper around MPI_Comm_size. */
int ORBIT_MPI_Comm_size(MPI_Comm comm, int * size){
  int res = 0;

#if USE_MPI > 0
  res = MPI_Comm_size(comm,size);
#else
  res  = 1;
  *size = 1;
#endif

  return res;
}

/** A C wrapper around MPI_Comm_rank. */
int ORBIT_MPI_Comm_rank(MPI_Comm comm, int * rank){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_rank(comm,rank);
#else
  res  = 1;
  *rank = 0;
#endif

  return res;
}

/** A C wrapper around MPI_Comm_free. */
int ORBIT_MPI_Comm_free(MPI_Comm* comm){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Comm_free(comm);
#else
  res  = 1;
#endif		
	return res;	
}

//--------------------------------------------------------
// MPI functions related to the MPI_Group manipulations
//--------------------------------------------------------

/** A C wrapper around MPI_Group_incl. */
int ORBIT_MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *group_out ){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Group_incl(group, n, ranks, group_out);
#else
  res  = 1;
#endif		
	return res;		
}
  
/** A C wrapper around MPI_Group_excl. */
int ORBIT_MPI_Group_excl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Group_excl(group, n, ranks, newgroup);
#else
  res  = 1;
#endif		
	return res;		
}

/** A C wrapper around MPI_Group_union. */
int ORBIT_MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *group_out){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Group_union(group1, group2, group_out);
#else
  res  = 1;
#endif		
	return res;	
}

/** A C wrapper around MPI_Group_difference. */
int ORBIT_MPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *group_out){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Group_difference(group1, group2, group_out);
#else
  res  = 1;
	*group_out = MPI_GROUP_EMPTY;
#endif		
	return res;
}

/** A C wrapper around MPI_Group_intersection. */
int ORBIT_MPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *group_out){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Group_intersection(group1, group2, group_out);
#else
  res  = 1;
	*group_out = MPI_GROUP_EMPTY;
#endif		
	return res;
}

/** A C wrapper around MPI_Group_compare. */
int ORBIT_MPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Group_compare(group1, group2, result);
#else
  res  = 1;
	*result = MPI_IDENT;
#endif		
	return res;		
}

/** A C wrapper around MPI_Group_translate_ranks. */
int ORBIT_MPI_Group_translate_ranks(MPI_Group group_a, int n, int *ranks_a, MPI_Group group_b, int *ranks_b){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Group_translate_ranks(group_a, n, ranks_a, group_b, ranks_b);
#else
  res  = 1;
#endif		
	return res;
}

/** A C wrapper around MPI_Group_size. */
int ORBIT_MPI_Group_size(MPI_Group group, int *size){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Group_size(group, size);
#else
  res  = 1;
	*size = 1; 
#endif		
	return res;		
}
	
/** A C wrapper around MPI_Group_rank. */
int ORBIT_MPI_Group_rank(MPI_Group group, int *rank){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Group_rank(group,rank);
#else
  res  = 1;
	*rank = 0; 
#endif		
	return res;		
}

/** A C wrapper around MPI_Group_free(. */
int ORBIT_MPI_Group_free(MPI_Group* group){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Group_free(group);
#else
  res  = 1;
#endif		
	return res;	
}

//--------------------------------------------------------
// MPI functions related to the MPI_Intercomm manipulations
//--------------------------------------------------------		

/** A C wrapper around MPI_Intercomm_create. */
int ORBIT_MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm, 
	int remote_leader, int tag, MPI_Comm *comm_out){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Intercomm_create(local_comm, local_leader, peer_comm, remote_leader, tag, comm_out);
#else
  res  = 1;
#endif		
	return res;	
}

/** A C wrapper around MPI_Intercomm_merge. */
int ORBIT_MPI_Intercomm_merge(MPI_Comm comm, int high, MPI_Comm *comm_out){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Intercomm_merge(comm, high, comm_out);
#else
  res  = 1;
#endif		
	return res;	
}

//--------------------------------------------------------
// MPI functions related to the MPI_Graph manipulations
//--------------------------------------------------------	

/** A C wrapper around MPI_Graph_create. */
int ORBIT_MPI_Graph_create(MPI_Comm comm_old, int nnodes, int *index, int *edges, 
	                         int reorder, MPI_Comm *comm_graph){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Graph_create(comm_old, nnodes, index, edges, reorder, comm_graph);
#else
  res  = 1;
#endif
  return res;
}

/** A C wrapper around MPI_Graphdims_get. */
int ORBIT_MPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Graphdims_get(comm, nnodes, nedges);
#else
  res  = 1;
	*nnodes = 1;
	*nedges = 1;
#endif
  return res;	
}

/** A C wrapper around MPI_Graph_get. */
int ORBIT_MPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int *index, int *edges){
  int res = 0;
#if USE_MPI > 0
  res = ORBIT_MPI_Graph_get(comm, maxindex, maxedges, index, edges);
#else
  res  = 1;
#endif
  return res;	
}

/** A C wrapper around MPI_Graph_map. */
int ORBIT_MPI_Graph_map(MPI_Comm comm_old, int nnodes, int *index, int *edges, int *newrank){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Graph_map(comm_old, nnodes, index, edges, newrank);
#else
  res  = 1;
#endif
  return res;	
}

/** A C wrapper around MPI_Graph_neighbors_count. */
int ORBIT_MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Graph_neighbors_count(comm, rank, nneighbors);
#else
  res  = 1;
#endif
  return res;		
}

/** A C wrapper around MPI_Graph_neighbors. */
int ORBIT_MPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int *neighbors){
	int res = 0;
#if USE_MPI > 0
  res = MPI_Graph_neighbors(comm, rank, maxneighbors, neighbors);
#else
  res  = 1;
#endif
  return res;	
}

//--------------------------------------------------------
// MPI functions related to the Send-Receive operations
//--------------------------------------------------------

/** A C wrapper around MPI_Barrier. */
int ORBIT_MPI_Barrier(MPI_Comm comm){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Barrier(comm);
#else
  res  = 1;
#endif
  return res;
}


/** A C wrapper around MPI_Wait. */
int ORBIT_MPI_Wait(MPI_Request  *request, MPI_Status *status){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Wait(request, status);
#else
  res  = 1;
#endif
  return res;
}

/** A C wrapper around MPI_Allreduce. */
int ORBIT_MPI_Allreduce(void* ar1, void* ar2, int n, MPI_Datatype data, MPI_Op op, MPI_Comm comm){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Allreduce(ar1, ar2, n, data, op, comm);
#else
  res  = 1;
#endif
  return res;
}

/** A C wrapper around MPI_Bcast. */
int ORBIT_MPI_Bcast(void* ar, int n1, MPI_Datatype data, int n2, MPI_Comm comm ){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Bcast(ar, n1, data, n2, comm ) ;
#else
  res  = 1;
#endif
  return res;
}

/** A C wrapper around MPI_Send. */
int ORBIT_MPI_Send(void* ar, int n1, MPI_Datatype data, int n2, int n3, MPI_Comm comm){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Send(ar, n1, data, n2, n3, comm);
#else
  res  = 1;
#endif
  return res;
}

/** A C wrapper around MPI_Recv. */
int ORBIT_MPI_Recv(void* ar, int n1, MPI_Datatype data, int n2, int n3, MPI_Comm comm, MPI_Status * stat){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Recv(ar, n1, data, n2, n3, comm, stat);
#else
  res  = 1;
#endif
  return res;
}

/** A C wrapper around MPI_Probe. */
int ORBIT_MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Probe(source, tag, comm, status);
#else
  res  = 1;
	status->count = 0;
	status->MPI_SOURCE = source;
	status->MPI_TAG = tag;
#endif
  return res;	
}

/** A C wrapper around MPI_Get_count. */
int ORBIT_MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count){
  int res = 0;
#if USE_MPI > 0
  res = MPI_Get_count(status, datatype, count);
#else
  res  = 1;
	status->count = 0;
#endif
  return res;		
}

