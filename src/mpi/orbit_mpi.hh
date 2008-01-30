#include "Python.h"

#ifndef ORBIT_MPI_INCLUDE
#define ORBIT_MPI_INCLUDE

#ifdef USE_MPI
 #include "mpi.h"
#else
//---------------------------------------------------------------
//START the case when USE_MPI is not defined.
//There is no MPI library. We will define everything by ourselves.
//----------------------------------------------------------------

 //Communicators ( handler)
 typedef int MPI_Comm;
 #define MPI_COMM_WORLD 91
 #define MPI_COMM_SELF  92

 //data type
 typedef int MPI_Datatype;
 #define MPI_CHAR           ((MPI_Datatype)1)
 #define MPI_UNSIGNED_CHAR  ((MPI_Datatype)2)
 #define MPI_BYTE           ((MPI_Datatype)3)
 #define MPI_SHORT          ((MPI_Datatype)4)
 #define MPI_UNSIGNED_SHORT ((MPI_Datatype)5)
 #define MPI_INT            ((MPI_Datatype)6)
 #define MPI_UNSIGNED       ((MPI_Datatype)7)
 #define MPI_LONG           ((MPI_Datatype)8)
 #define MPI_UNSIGNED_LONG  ((MPI_Datatype)9)
 #define MPI_FLOAT          ((MPI_Datatype)10)
 #define MPI_DOUBLE         ((MPI_Datatype)11)
 #define MPI_LONG_DOUBLE    ((MPI_Datatype)12)
 #define MPI_LONG_LONG_INT  ((MPI_Datatype)13)

 // Groups
 typedef int MPI_Group;
 #define MPI_GROUP_EMPTY 90

 //Operations
 typedef int MPI_Op;
 #define MPI_MAX    (MPI_Op)(100)
 #define MPI_MIN    (MPI_Op)(101)
 #define MPI_SUM    (MPI_Op)(102)
 #define MPI_PROD   (MPI_Op)(103)
 #define MPI_LAND   (MPI_Op)(104)
 #define MPI_BAND   (MPI_Op)(105)
 #define MPI_LOR    (MPI_Op)(106)
 #define MPI_BOR    (MPI_Op)(107)
 #define MPI_LXOR   (MPI_Op)(108)
 #define MPI_BXOR   (MPI_Op)(109)
 #define MPI_MINLOC (MPI_Op)(110)
 #define MPI_MAXLOC (MPI_Op)(111)

 typedef struct {
    int count;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
 } MPI_Status;

 // Request ( handler)
 typedef int MPI_Request; 
 
 /* Define some null objects */
 #define MPI_COMM_NULL      ((MPI_Comm)0)
 #define MPI_OP_NULL        ((MPI_Op)0)
 #define MPI_GROUP_NULL     ((MPI_Group)0)
 #define MPI_DATATYPE_NULL  ((MPI_Datatype)0)

 /* Topology types */
 #define MPI_GRAPH  1
 #define MPI_CART   2

 /* Results of the compare operations */   
 #define MPI_IDENT     0  
 #define MPI_CONGRUENT 1  
 #define MPI_SIMILAR   2  
 #define MPI_UNEQUAL   3 

 /* Names' lengths */
 #define MPI_MAX_PROCESSOR_NAME 256
 #define MPI_MAX_ERROR_STRING   512
 #define MPI_MAX_NAME_STRING     63

 /* MPI Constants */
 #define MPI_UNDEFINED      (-32766)
 #define MPI_UNDEFINED_RANK MPI_UNDEFINED
 #define MPI_SUCCESS              0
 #define MPI_ANY_SOURCE         (-2)   
 #define MPI_ANY_TAG            (-1)

#endif
//-------------------------------------------------------------
//END the case when USE_MPI is defined or not.
//--------------------------------------------------------------


//--------------------------------------------------------------
//     The PyORBIT MPI classes definitions         START
//--------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif

 typedef struct {
   PyObject_HEAD
   MPI_Comm comm;
 } pyORBIT_MPI_Comm;
 
 typedef struct {
   PyObject_HEAD
   MPI_Group group;
 } pyORBIT_MPI_Group; 
 
 typedef struct {
   PyObject_HEAD
   MPI_Status status;
 } pyORBIT_MPI_Status;
 
 typedef struct {
   PyObject_HEAD
   MPI_Request request;
 } pyORBIT_MPI_Request;  
 
 typedef struct {
   PyObject_HEAD
   MPI_Datatype datatype;
 } pyORBIT_MPI_Datatype;
 
 typedef struct {
   PyObject_HEAD
   MPI_Op op;
 } pyORBIT_MPI_Op; 
 
#ifdef __cplusplus
}
#endif
//--------------------------------------------------------------
//     The PyORBIT MPI classes definitions         STOP
//--------------------------------------------------------------


int ORBIT_MPI_Init(int *len, char ***ch);
int ORBIT_MPI_Initialized(int *init);
int ORBIT_MPI_Finalize();
int ORBIT_MPI_Finalize(const char* message);
int ORBIT_MPI_Get_processor_name(char *name, int* len);
double ORBIT_MPI_Wtime(void);
double ORBIT_MPI_Wtick();

//--------------------------------------------------------
// MPI functions related to the MPI_Comm manipulations
//--------------------------------------------------------	
int ORBIT_MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *comm_out);
int ORBIT_MPI_Comm_group(MPI_Comm comm, MPI_Group *group );
int ORBIT_MPI_Comm_dup(MPI_Comm comm, MPI_Comm *comm_out);
int ORBIT_MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *comm_out);
int ORBIT_MPI_Comm_remote_size(MPI_Comm comm, int *size);	
int ORBIT_MPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group);
int ORBIT_MPI_Comm_test_inter(MPI_Comm comm, int *flag);
int ORBIT_MPI_Comm_compare(MPI_Comm  comm1, MPI_Comm  comm2, int *result);
int ORBIT_MPI_Comm_set_name(MPI_Comm com, char *name);
int ORBIT_MPI_Comm_get_name(MPI_Comm comm, char *namep, int *reslen);
int ORBIT_MPI_Comm_size(MPI_Comm comm, int * size);
int ORBIT_MPI_Comm_rank(MPI_Comm comm, int * rank);
int ORBIT_MPI_Comm_free(MPI_Comm* comm);

//--------------------------------------------------------
// MPI functions related to the MPI_Group manipulations
//--------------------------------------------------------	
int ORBIT_MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *group_out );
int ORBIT_MPI_Group_excl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
int ORBIT_MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *group_out);	
int ORBIT_MPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *group_out);
int ORBIT_MPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *group_out);
int ORBIT_MPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result);	
int ORBIT_MPI_Group_translate_ranks(MPI_Group group_a, int n, int *ranks_a, MPI_Group group_b, int *ranks_b);	
int ORBIT_MPI_Group_size(MPI_Group group, int *size);
int ORBIT_MPI_Group_rank(MPI_Group group, int *rank);
int ORBIT_MPI_Group_free(MPI_Group* group);

//--------------------------------------------------------
// MPI functions related to the MPI_Intercomm manipulations
//--------------------------------------------------------		
int ORBIT_MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm, 
															 int remote_leader, int tag, MPI_Comm *comm_out);	
int ORBIT_MPI_Intercomm_merge(MPI_Comm comm, int high, MPI_Comm *comm_out);

//--------------------------------------------------------
// MPI functions related to the Graph manipulations
//--------------------------------------------------------	
int ORBIT_MPI_Graph_create(MPI_Comm comm_old, int nnodes, int *index, int *edges, int reorder, MPI_Comm *comm_graph);
int ORBIT_MPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges);
int ORBIT_MPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int *index, int *edges);
int ORBIT_MPI_Graph_map(MPI_Comm comm_old, int nnodes, int *index, int *edges, int *newrank);
int ORBIT_MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors);
int ORBIT_MPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int *neighbors);

//--------------------------------------------------------
// MPI functions related to the Send-Receive operations
//--------------------------------------------------------
int ORBIT_MPI_Barrier(MPI_Comm comm);
int ORBIT_MPI_Wait(MPI_Request  *request, MPI_Status *status);
int ORBIT_MPI_Allreduce(void* buf_in, void* buf_out, int count, MPI_Datatype, MPI_Op, MPI_Comm);
int ORBIT_MPI_Bcast(void* buf, int count, MPI_Datatype, int rank, MPI_Comm);	
int ORBIT_MPI_Send(void* buf, int count, MPI_Datatype, int dest,   int tag, MPI_Comm);
int ORBIT_MPI_Recv(void* buf, int count, MPI_Datatype, int source, int tag, MPI_Comm, MPI_Status *);
int ORBIT_MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);
int ORBIT_MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);	


#endif   //end of ---ifndef ORBIT_MPI_INCLUDE---
