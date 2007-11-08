#include "Python.h"

#ifndef ORBIT_MPI_INCLUDE
#define ORBIT_MPI_INCLUDE

#ifdef USE_MPI
 #include "mpi.h"
#else
//-------------------------------------------------------------
//START  #ifdef USE_MPI
//Defenition of the MPI realted types for the case when there 
//is no MPI library
//--------------------------------------------------------------

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

#endif
//-------------------------------------------------------------
//END of  #ifdef USE_MPI
//Defenition of the MPI realted types for the case when there 
//is no MPI library
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
  int ORBIT_MPI_Finalize(void);
  int ORBIT_MPI_Finalize(const char* message);
  int ORBIT_MPI_Get_processor_name(char *name, int* len);
  int ORBIT_MPI_Comm_size(MPI_Comm comm, int * size);
  int ORBIT_MPI_Comm_rank(MPI_Comm comm, int * rank);

	int ORBIT_MPI_Comm_free(MPI_Comm* comm);
	
	int ORBIT_MPI_Group_free(MPI_Group* group);
	
  double ORBIT_MPI_Wtime(void);

  int ORBIT_MPI_Allreduce(void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm);

  int ORBIT_MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
  int ORBIT_MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);

  int ORBIT_MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm );

#endif   //end of ---ifndef ORBIT_MPI_INCLUDE---
