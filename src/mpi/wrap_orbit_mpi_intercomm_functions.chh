//----------------------------------------------------------------------
// Wrappers for MPI functions related to the MPI_Intercomm manipulations
//----------------------------------------------------------------------
static PyObject* mpi_intercomm_create(PyObject *self, PyObject *args){
	PyObject* pyO_0; PyObject* pyO_1; PyObject* pyO_out;
	int local_leader, remote_leader, tag;
	if(!PyArg_ParseTuple(	args,"OiOiiO:mpi_intercomm_create",&pyO_0,&local_leader,&pyO_1,&remote_leader,&tag,&pyO_out)){
		ORBIT_MPI_Finalize("ORBIT_MPI_Intercomm_create(MPI_Comm local,int loc_leader,MPI_Comm local, int rem_leader, int tag, MPI_Comm out) - needs 6 params.");
	}
	pyORBIT_MPI_Comm* pyComm_local = (pyORBIT_MPI_Comm*) pyO_0;
	pyORBIT_MPI_Comm* pyComm_remote = (pyORBIT_MPI_Comm*) pyO_1;
	pyORBIT_MPI_Comm* pyComm_out = (pyORBIT_MPI_Comm*) pyO_out;
	int res = ORBIT_MPI_Intercomm_create(pyComm_local->comm, local_leader, pyComm_remote->comm,remote_leader,tag,&pyComm_out->comm);
	return Py_BuildValue("i",res);
}

static PyObject* mpi_intercomm_merge(PyObject *self, PyObject *args){
	PyObject* pyO_in; PyObject* pyO_out;
	int high;
	if(!PyArg_ParseTuple(	args,"OiO:mpi_intercomm_merge",&pyO_in,&high,&pyO_out)){
		ORBIT_MPI_Finalize("ORBIT_MPI_Intercomm_merge(MPI_Comm in,int high,MPI_Comm out) - needs 3 params.");
	}
	pyORBIT_MPI_Comm* pyComm_in = (pyORBIT_MPI_Comm*) pyO_in;
	pyORBIT_MPI_Comm* pyComm_out = (pyORBIT_MPI_Comm*) pyO_out;
	return Py_BuildValue("i",ORBIT_MPI_Intercomm_merge(pyComm_in->comm, high, &pyComm_out->comm));
}
