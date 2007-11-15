#include "Python.h"
#include "orbit_mpi.hh"

#include <cstring>
#include <iostream>

//arrays buffers for MPI exchange
#include "BufferStore.hh"

#include "wrap_orbit_mpi.hh"

//wrappers of mpi objects
#include "wrap_mpi_comm.hh"
#include "wrap_mpi_group.hh"
#include "wrap_mpi_status.hh"
#include "wrap_mpi_request.hh"
#include "wrap_mpi_datatype.hh"
#include "wrap_mpi_op.hh"


namespace wrap_orbit_mpi{
	
  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }
	
	#ifdef __cplusplus
	extern "C" {
		#endif
		
		static PyObject* mpi_initialized(PyObject *self, PyObject *args) {	
      int result = 0;
      ORBIT_MPI_Initialized(&result);
			return Py_BuildValue("i",result);		 
		}
		
    static PyObject* mpi_get_processor_name(PyObject *self, PyObject *args) {
      char* name = new char[MPI_MAX_PROCESSOR_NAME];
      int len;
      ORBIT_MPI_Get_processor_name(name,&len);
      //PyObject* nm = PyString_FromStringAndSize(name, len);
			PyObject* nm = Py_BuildValue("s#",name, len);
      delete [] name;
      return nm;
    }

    static PyObject* mpi_wtime(PyObject *self, PyObject *args) {
      double result ;
      result = (double) ORBIT_MPI_Wtime();
      return Py_BuildValue("d",result);
    }
		
    static PyObject* mpi_wtick(PyObject *self, PyObject *args) {
      double result ;
      result = (double) ORBIT_MPI_Wtick();
      return Py_BuildValue("d",result);
    }
	
		//------------------------------------------------------------------
		// Wrappers for MPI functions related to the MPI_Comm manipulations
		//------------------------------------------------------------------		
		#include 	"wrap_orbit_mpi_comm_functions.hh"	
		
		//------------------------------------------------------------------
		// Wrappers for MPI functions related to the MPI_Group manipulations
		//------------------------------------------------------------------		
		#include 	"wrap_orbit_mpi_group_functions.hh"	

		//----------------------------------------------------------------------
		// Wrappers for MPI functions related to the MPI_Intercomm manipulations
		//----------------------------------------------------------------------
		#include 	"wrap_orbit_mpi_intercomm_functions.hh"	
	
		//-------------------------------------------------------------------
		// Wrappers for MPI functions related to the MPI_Graph manipulations
		//-------------------------------------------------------------------		
		#include 	"wrap_orbit_mpi_graph_functions.hh"			
		
		//------------------------------------------------------------------
		// Wrappers for MPI functions related to the Send-Receive operations
		//------------------------------------------------------------------		
		#include 	"wrap_orbit_mpi_send_receive_functions.hh"	
		
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
			if(nVars == 0 ||  nVars == 1){
				if(nVars == 1){
					char* message = NULL;
					//NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
					if(!PyArg_ParseTuple(	args,"s:finalize",&message)){
						ORBIT_MPI_Finalize("orbit_mpi - something wrong with error message.");
					}
					ORBIT_MPI_Finalize(message);
				}
				ORBIT_MPI_Finalize();
			}
			else{
				ORBIT_MPI_Finalize("orbit_mpi. You should call finalize() or finalize(message)");
			}
			
			Py_INCREF(Py_None);
			return Py_None;
		}
		
		static PyMethodDef orbit_mpiMethods[] = {
			{"MPI_Initialized",         mpi_initialized,          METH_VARARGS },			
			{"MPI_Get_processor_name",  mpi_get_processor_name,   METH_VARARGS },
			{"MPI_Wtime",               mpi_wtime,                METH_VARARGS },
			{"MPI_Wtick",               mpi_wtick,                METH_VARARGS },	
			{"finalize",                finalize,                 METH_VARARGS },
			// Wrappers for MPI functions related to the MPI_Comm manipulations
			{"MPI_Comm_create",         mpi_comm_create,          METH_VARARGS },			
			{"MPI_Comm_group",          mpi_comm_group,           METH_VARARGS },			
			{"MPI_Comm_dup",            mpi_comm_dup,             METH_VARARGS },			
			{"MPI_Comm_split",          mpi_comm_split,           METH_VARARGS },			
			{"MPI_Comm_remote_size",    mpi_comm_remote_size,     METH_VARARGS },			
			{"MPI_Comm_remote_group",   mpi_comm_remote_group,    METH_VARARGS },			
			{"MPI_Comm_test_inter",     mpi_comm_test_inter,      METH_VARARGS },			
			{"MPI_Comm_compare",        mpi_comm_compare,         METH_VARARGS },			
			{"MPI_Comm_set_name",       mpi_comm_set_name,        METH_VARARGS },			
			{"MPI_Comm_get_name",       mpi_comm_get_name,        METH_VARARGS },			
			{"MPI_Comm_size",           mpi_comm_size,            METH_VARARGS },
			{"MPI_Comm_rank",           mpi_comm_rank,            METH_VARARGS },
			{"MPI_Comm_free",           mpi_comm_free,            METH_VARARGS },		
			// Wrappers for MPI functions related to the MPI_Group manipulations
			{"MPI_Group_incl",            mpi_group_incl,            METH_VARARGS },		
			{"MPI_Group_excl",            mpi_group_excl,            METH_VARARGS },		
			{"MPI_Group_union",           mpi_group_union,           METH_VARARGS },		
			{"MPI_Group_difference",      mpi_group_difference,      METH_VARARGS },		
			{"MPI_Group_intersection",    mpi_group_intersection,    METH_VARARGS },		
			{"MPI_Group_compare",         mpi_group_compare,         METH_VARARGS },		
			{"MPI_Group_translate_ranks", mpi_group_translate_ranks, METH_VARARGS },		
			{"MPI_Group_size",            mpi_group_size,            METH_VARARGS },		
			{"MPI_Group_rank",            mpi_group_rank,            METH_VARARGS },		
			{"MPI_Group_free",            mpi_group_free,            METH_VARARGS },
			// Wrappers for MPI functions related to the MPI_Intercomm manipulations
			{"MPI_Intercomm_create",      mpi_intercomm_create,      METH_VARARGS },
			{"MPI_Intercomm_merge",       mpi_intercomm_merge,       METH_VARARGS },
			// Wrappers for MPI functions related to the MPI_Graph manipulations			
			{"MPI_Graph_create",          mpi_graph_create,          METH_VARARGS },			
			{"MPI_Graphdims_get",         mpi_graphdims_get,         METH_VARARGS },			
			{"MPI_Graph_get",             mpi_graph_get,             METH_VARARGS },			
			{"MPI_Graph_map",             mpi_graph_map,             METH_VARARGS },			
			{"MPI_Graph_neighbors_count", mpi_graph_neighbors_count, METH_VARARGS },			
			{"MPI_Graph_neighbors",       mpi_graph_neighbors,       METH_VARARGS },	
			// Wrappers for MPI functions related to the Send-Receive operations
			{"MPI_Barrier",               mpi_barrier,               METH_VARARGS },	
			{"MPI_Wait",                  mpi_wait,                  METH_VARARGS },	
			{"MPI_Allreduce",             mpi_allreduce,             METH_VARARGS },	
			{"MPI_Bcast",                 mpi_bcast,                 METH_VARARGS },	
			{"MPI_Send",                  mpi_send,                  METH_VARARGS },	
			{"MPI_Recv",                  mpi_recv,                  METH_VARARGS },	
			{ NULL, NULL }
		};
		
		void initorbit_mpi(void) {
			PyObject *m, *d;
			m = Py_InitModule((char*)"orbit_mpi",orbit_mpiMethods);
			d = PyModule_GetDict(m);

			//add the results of comparisons	constants	
	    PyModule_AddObject(m,"MPI_IDENT", Py_BuildValue("i",(int) MPI_IDENT));
	    PyModule_AddObject(m,"MPI_CONGRUENT", Py_BuildValue("i",(int) MPI_CONGRUENT));
	    PyModule_AddObject(m,"MPI_SIMILAR", Py_BuildValue("i",(int) MPI_SIMILAR));
	    PyModule_AddObject(m,"MPI_UNEQUAL", Py_BuildValue("i",(int) MPI_UNEQUAL));
			
			//add constants
			PyModule_AddObject(m,"MPI_UNDEFINED", Py_BuildValue("i",(int) MPI_UNDEFINED));
			PyModule_AddObject(m,"MPI_UNDEFINED_RANK", Py_BuildValue("i",(int) MPI_UNDEFINED_RANK));
			PyModule_AddObject(m,"MPI_SUCCESS", Py_BuildValue("i",(int) MPI_SUCCESS));
			PyModule_AddObject(m,"MPI_ANY_SOURCE", Py_BuildValue("i",(int) MPI_ANY_SOURCE));
			PyModule_AddObject(m,"MPI_ANY_TAG", Py_BuildValue("i",(int) MPI_ANY_TAG));
	
			//add MPI_Comm class and fields
			wrap_orbit_mpi_comm::init_orbit_mpi_comm(m);
			
			//add MPI_Group class and fields
			wrap_orbit_mpi_group::init_orbit_mpi_group(m);
			
			//add MPI_Status class and fields
			wrap_orbit_mpi_status::init_orbit_mpi_status(m);
			
			//add MPI_Request class and fields
			wrap_orbit_mpi_request::init_orbit_mpi_request(m);
			
			//add MPI_Datatype class and fields
			wrap_orbit_mpi_datatype::init_orbit_mpi_datatype(m);
			
			//add MPI_Op class and fields
			wrap_orbit_mpi_op::init_orbit_mpi_op(m);
		}
		
		
		#ifdef __cplusplus
	}
	#endif
	
	
	//end of namespace
}
