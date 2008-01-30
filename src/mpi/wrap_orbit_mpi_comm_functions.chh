//------------------------------------------------------------------
// Wrappers for MPI functions related to the MPI_Comm manipulations
//------------------------------------------------------------------
static PyObject* mpi_comm_create(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 2){
		error("MPI_Comm_create(MPI_Comm, MPI_Group) needs 2 parameters.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
	pyORBIT_MPI_Group* pyGroup = (pyORBIT_MPI_Group*) PyTuple_GetItem(args,1);
	pyORBIT_MPI_Comm* pyComm_out = wrap_orbit_mpi_comm::newMPI_Comm();
	if(ORBIT_MPI_Comm_create(pyComm->comm,pyGroup->group,&pyComm_out->comm) != MPI_SUCCESS){
		error("MPI_Comm_create(MPI_Comm, MPI_Group) - fatal error. STOP.");
	}
	return (PyObject*) pyComm_out;
}

static PyObject* mpi_comm_group(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Comm_group(MPI_Comm) needs one parameter.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
	pyORBIT_MPI_Group* pyGroup = wrap_orbit_mpi_group::newMPI_Group();
	if(ORBIT_MPI_Comm_group(pyComm->comm,&pyGroup->group) != MPI_SUCCESS){
		error("MPI_Comm_group(MPI_Comm) - fatal error. STOP.");
	}
	return (PyObject*)  pyGroup;
}

static PyObject* mpi_comm_dup(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Comm_dup(MPI_Comm) needs one parameter.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
	pyORBIT_MPI_Comm* pyComm_out = wrap_orbit_mpi_comm::newMPI_Comm();
	if(ORBIT_MPI_Comm_dup(pyComm->comm,&pyComm_out->comm) != MPI_SUCCESS){
		error("MPI_Comm_dup(MPI_Comm) - fatal error. STOP.");
	}
	return (PyObject*) pyComm_out;
}

static PyObject* mpi_comm_split(PyObject *self, PyObject *args){
	PyObject* pyO_in;
	int color, key;
	if(!PyArg_ParseTuple(	args,"Oii:mpi_comm_split",&pyO_in,&color,&key)){
		ORBIT_MPI_Finalize("ORBIT_MPI_Comm_split(MPI_Comm in,int color,int key) - needs 3 params.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO_in;
	pyORBIT_MPI_Comm* pyComm_out = wrap_orbit_mpi_comm::newMPI_Comm();
	if(ORBIT_MPI_Comm_split(pyComm->comm,color,key,&pyComm_out->comm)  != MPI_SUCCESS){
		ORBIT_MPI_Finalize("ORBIT_MPI_Comm_split(...) - fatal error. STOP.");
	}
	return (PyObject*) pyComm_out;
}

static PyObject* mpi_comm_remote_size(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Comm_remote_size(MPI_Comm) needs MPI_Comm as a parameter.");
	}			
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
	int result = 0;
	if(ORBIT_MPI_Comm_remote_size(pyComm->comm,&result) != MPI_SUCCESS){
		error("MPI_Comm_remote_size(MPI_Comm) - fatal error. STOP.");
	}
	return Py_BuildValue("i",result);
}

static PyObject* mpi_comm_remote_group(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Comm_remote_group(MPI_Comm) needs one parameter.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
	pyORBIT_MPI_Group* pyGroup = wrap_orbit_mpi_group::newMPI_Group();
	if(ORBIT_MPI_Comm_remote_group(pyComm->comm,&pyGroup->group) != MPI_SUCCESS){
		error("MPI_Comm_remote_group(MPI_Comm) - fatal error. STOP.");
	}
	return (PyObject*) pyGroup;
}

static PyObject* mpi_comm_test_inter(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Comm_test_inter(MPI_Comm) needs MPI_Comm as a parameter.");
	}			
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
	int result = 0;
	if(ORBIT_MPI_Comm_test_inter(pyComm->comm,&result) != MPI_SUCCESS){
		error("MPI_Comm_test_inter(MPI_Comm) - fatal error. STOP.");
	}
	return Py_BuildValue("i",result);
}

static PyObject* mpi_comm_compare(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 2){
		error("MPI_Comm_compare(MPI_Comm, MPI_Comm) needs two parameters.");
	}
	pyORBIT_MPI_Comm* pyComm0 = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
	pyORBIT_MPI_Comm* pyComm1 = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,1);
	int result = 0;
	if(ORBIT_MPI_Comm_compare(pyComm0->comm,pyComm1->comm,&result) != MPI_SUCCESS){
		error("MPI_Comm_compare(MPI_Comm, MPI_Comm) - fatal error. STOP.");
	}
	return Py_BuildValue("i",result);
}

static PyObject* mpi_comm_set_name(PyObject *self, PyObject *args) {
	PyObject* pyO;
	char* comm_name = NULL;
	if(!PyArg_ParseTuple(	args,"Os:mpi_comm_set_name",&pyO,&comm_name)){
		error("MPI_Comm_set_name(MPI_Comm, name) needs 2 parameters.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO;
	return Py_BuildValue("i",ORBIT_MPI_Comm_set_name(pyComm->comm, comm_name));	
}

static PyObject* mpi_comm_get_name(PyObject *self, PyObject *args) {
	PyObject* pyO;
	if(!PyArg_ParseTuple(	args,"O:mpi_comm_get_name",&pyO)){
		error("MPI_Comm_get_name(MPI_Comm) needs 1 parameter.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO;
	char* comm_name = new char[MPI_MAX_NAME_STRING];
	int len = 0;	
	ORBIT_MPI_Comm_get_name(pyComm->comm, comm_name, &len);
	PyObject* name = Py_BuildValue("s#",comm_name, len);
	delete [] comm_name;
	return name;
}

static PyObject* mpi_comm_size(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Comm_size(MPI_Comm) needs MPI_Comm as a parameter.");
	}			
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
	int result = 0;
	ORBIT_MPI_Comm_size(pyComm->comm,&result);
	return Py_BuildValue("i",result);
}

static PyObject* mpi_comm_rank(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Comm_rank(MPI_Comm) needs MPI_Comm as a parameter.");
	}			
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
	int result = 0;
	ORBIT_MPI_Comm_rank(pyComm->comm,&result);
	return Py_BuildValue("i",result);
}

static PyObject* mpi_comm_free(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Comm_free(MPI_Comm) needs MPI_Comm as an input parameter.");
	}	
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
	return Py_BuildValue("i",ORBIT_MPI_Comm_free(&pyComm->comm));
}

		

