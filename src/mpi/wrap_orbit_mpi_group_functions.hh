//------------------------------------------------------------------
// Wrappers for MPI functions related to the MPI_Group manipulations
//------------------------------------------------------------------
static PyObject* mpi_group_incl(PyObject *self, PyObject *args) {
	PyObject* pyO_in; PyObject* pyO_arr;
	if(!PyArg_ParseTuple(args,"OO:mpi_group_incl",&pyO_in,&pyO_arr)){	
		error("MPI_Group_incl(MPI_Group in, [...ranks]) needs 2 parameters.");
	}
	pyORBIT_MPI_Group* pyGroup_in = (pyORBIT_MPI_Group*) pyO_in;
	pyORBIT_MPI_Group* pyGroup_out = wrap_orbit_mpi_group::newMPI_Group();
  if(!PySequence_Check(pyO_arr)){
		error("MPI_Group_incl(MPI_Group in, [...ranks]), [...ranks] should be a sequence.");
	}
	int size = PySequence_Size(pyO_arr);
	int buff_index = 0;
	int* i_arr = BufferStore::getBufferStore()->getFreeIntArr(buff_index,size);
	for(int i = 0; i < size; i++){
		i_arr[i]= (int) PyInt_AsLong(PySequence_Fast_GET_ITEM(pyO_arr, i));
	}
	if(ORBIT_MPI_Group_incl(pyGroup_in->group,size,i_arr,&pyGroup_out->group)  != MPI_SUCCESS){
		error("MPI_Group_incl(MPI_Group in, [...ranks]) - fatal error. STOP.");
	}
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index);
	return (PyObject*) pyGroup_out;
}

static PyObject* mpi_group_excl(PyObject *self, PyObject *args) {
	PyObject* pyO_in; PyObject* pyO_arr;
	if(!PyArg_ParseTuple(args,"OO:mpi_group_excl",&pyO_in,&pyO_arr)){	
		error("MPI_Group_excl(MPI_Group in, [...ranks]) needs 2 parameters.");
	}
	pyORBIT_MPI_Group* pyGroup_in = (pyORBIT_MPI_Group*) pyO_in;
	pyORBIT_MPI_Group* pyGroup_out = wrap_orbit_mpi_group::newMPI_Group();
  if(!PySequence_Check(pyO_arr)){
		error("MPI_Group_excl(MPI_Group in, [...ranks]), [...ranks] should be a sequence.");
	}
	int size = PySequence_Size(pyO_arr);
	int buff_index = 0;
	int* i_arr = BufferStore::getBufferStore()->getFreeIntArr(buff_index,size);
	for(int i = 0; i < size; i++){
		i_arr[i]= (int) PyInt_AsLong(PySequence_Fast_GET_ITEM(pyO_arr, i));
	}
	if(ORBIT_MPI_Group_excl(pyGroup_in->group,size,i_arr,&pyGroup_out->group) != MPI_SUCCESS){
		error("MPI_Group_excl(MPI_Group in, [...ranks]) - fatal error. STOP.");
	}
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index);
	return (PyObject*) pyGroup_out;
}

static PyObject* mpi_group_union(PyObject *self, PyObject *args) {
	PyObject* pyO_in0; PyObject* pyO_in1;
	if(!PyArg_ParseTuple(args,"OO:mpi_group_union",&pyO_in0,&pyO_in1)){	
		error("MPI_Group_union(MPI_Group group_0,MPI_Group group_1) needs 2 parameters.");
	}
	pyORBIT_MPI_Group* pyGroup_in0 = (pyORBIT_MPI_Group*) pyO_in0;
	pyORBIT_MPI_Group* pyGroup_in1 = (pyORBIT_MPI_Group*) pyO_in1;
	pyORBIT_MPI_Group* pyGroup_out = wrap_orbit_mpi_group::newMPI_Group();
	if(ORBIT_MPI_Group_union(pyGroup_in0->group,pyGroup_in1->group,&pyGroup_out->group) != MPI_SUCCESS){
		error("MPI_Group_union(MPI_Group group_0,MPI_Group group_1) - fatal error. STOP.");
	}
	return (PyObject*) pyGroup_out;
}

static PyObject* mpi_group_difference(PyObject *self, PyObject *args) {
	PyObject* pyO_in0; PyObject* pyO_in1;
	if(!PyArg_ParseTuple(args,"OO:mpi_group_difference",&pyO_in0,&pyO_in1)){	
		error("MPI_Group_difference(MPI_Group in0,MPI_Group in1) needs 2 parameters.");
	}
	pyORBIT_MPI_Group* pyGroup_in0 = (pyORBIT_MPI_Group*) pyO_in0;
	pyORBIT_MPI_Group* pyGroup_in1 = (pyORBIT_MPI_Group*) pyO_in1;
	pyORBIT_MPI_Group* pyGroup_out = wrap_orbit_mpi_group::newMPI_Group();
	if(ORBIT_MPI_Group_difference(pyGroup_in0->group,pyGroup_in1->group,&pyGroup_out->group) != MPI_SUCCESS){
		error("MPI_Group_difference(MPI_Group in0,MPI_Group in1) - fatal error. STOP.");
	}
	return (PyObject*) pyGroup_out;
}

static PyObject* mpi_group_intersection(PyObject *self, PyObject *args) {
	PyObject* pyO_in0; PyObject* pyO_in1;
	if(!PyArg_ParseTuple(args,"OO:mpi_group_intersection",&pyO_in0,&pyO_in1)){	
		error("MPI_Group_intersection(MPI_Group in_0,MPI_Group in_1) needs 2 parameters.");
	}
	pyORBIT_MPI_Group* pyGroup_in0 = (pyORBIT_MPI_Group*) pyO_in0;
	pyORBIT_MPI_Group* pyGroup_in1 = (pyORBIT_MPI_Group*) pyO_in1;
	pyORBIT_MPI_Group* pyGroup_out = wrap_orbit_mpi_group::newMPI_Group();
	if(ORBIT_MPI_Group_intersection(pyGroup_in0->group,pyGroup_in1->group,&pyGroup_out->group) != MPI_SUCCESS){
		error("MPI_Group_intersection(MPI_Group in_0,MPI_Group in_1) - fatal error. STOP.");
	}
	return (PyObject*) pyGroup_out;
}

static PyObject* mpi_group_compare(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 2){
		error("MPI_Group_compare(MPI_Group, MPI_Group) needs two parameters.");
	}
	pyORBIT_MPI_Group* pyGroup0 = (pyORBIT_MPI_Group*) PyTuple_GetItem(args,0);
	pyORBIT_MPI_Group* pyGroup1 = (pyORBIT_MPI_Group*) PyTuple_GetItem(args,1);
	int result = 0;
	if(ORBIT_MPI_Group_compare(pyGroup0->group,pyGroup1->group,&result) != MPI_SUCCESS){
		error("MPI_Group_compare(MPI_Group, MPI_Group) fata error. STOP.");
	}
	return Py_BuildValue("i",result);
}

static PyObject* mpi_group_translate_ranks(PyObject *self, PyObject *args) {
	PyObject* pyO_a; PyObject* pyO_arr; PyObject* pyO_b;
	
	if(!PyArg_ParseTuple(args,"OOO:mpi_group_translate_ranks",&pyO_a,&pyO_arr,&pyO_b)){	
		error("MPI_Group_translate_ranks(MPI_Group a, [...ranks],MPI_Group b) needs 3 parameters.");
	}
	pyORBIT_MPI_Group* pyGroup_a = (pyORBIT_MPI_Group*) pyO_a;
	pyORBIT_MPI_Group* pyGroup_b = (pyORBIT_MPI_Group*) pyO_b;
  if(!PySequence_Check(pyO_arr)){
		error("MPI_Group_translate_ranks(MPI_Group in, [...ranks],MPI_Group out), [...ranks] should be a sequence.");
	}
	int size = PySequence_Size(pyO_arr);
	int buff_index0 = 0;
	int buff_index1 = 0;
	int* i_arr_in = BufferStore::getBufferStore()->getFreeIntArr(buff_index0,size);
	int* i_arr_out = BufferStore::getBufferStore()->getFreeIntArr(buff_index1,size);
	for(int i = 0; i < size; i++){
		i_arr_in[i] = (int) PyInt_AsLong(PySequence_Fast_GET_ITEM(pyO_arr, i));
		i_arr_out[i] = i_arr_in[i];
	}
	ORBIT_MPI_Group_translate_ranks(pyGroup_a->group, size, i_arr_in, pyGroup_b->group, i_arr_out);
	PyObject* pyArr = PyList_New(size);
	for(int i = 0; i < size; i++){
		PyList_Append(pyArr, Py_BuildValue("i",i_arr_out[i]));
	}
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index0);
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index1);
	return pyArr;
}

static PyObject* mpi_group_size(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Group_size(MPI_Group) needs MPI_Group as a parameter.");
	}			
	pyORBIT_MPI_Group* pyGroup = (pyORBIT_MPI_Group*) PyTuple_GetItem(args,0);
	int result = 0;
	ORBIT_MPI_Group_size(pyGroup->group,&result);
	return Py_BuildValue("i",result);
}

static PyObject* mpi_group_rank(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Group_rank(MPI_Group) needs MPI_Group as a parameter.");
	}			
	pyORBIT_MPI_Group* pyGroup = (pyORBIT_MPI_Group*) PyTuple_GetItem(args,0);
	int result = 0;
	ORBIT_MPI_Group_rank(pyGroup->group,&result);
	return Py_BuildValue("i",result);
}

static PyObject* mpi_group_free(PyObject *self, PyObject *args) {
	if(PyTuple_Size(args) != 1){
		error("MPI_Group_free(MPI_Group) needs MPI_Group as an input parameter.");
	}	
	pyORBIT_MPI_Group* pyGroup = (pyORBIT_MPI_Group*) PyTuple_GetItem(args,0);
	return Py_BuildValue("i",ORBIT_MPI_Group_free(&pyGroup->group));
}


