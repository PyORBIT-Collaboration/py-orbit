//-------------------------------------------------------------------
// Wrappers for MPI functions related to the MPI_Graph manipulations
//-------------------------------------------------------------------
static PyObject* mpi_graph_create(PyObject *self, PyObject *args) {
	PyObject* pyO_old; PyObject* pyO_indexes; PyObject* pyO_edges; PyObject* pyO_new;
	int reorder;
	if(!PyArg_ParseTuple(args,"OOOiO:mpi_graph_create",&pyO_old,&pyO_indexes,&pyO_edges,&reorder,&pyO_new)){	
		error("MPI_Graph_create(MPI_Comm old, [...indexes],[...egdes],reorder,MPI_Comm graph) needs 6 parameters.");
	}
	pyORBIT_MPI_Comm* pyComm_old = (pyORBIT_MPI_Comm*) pyO_old;
	pyORBIT_MPI_Comm* pyComm_new = (pyORBIT_MPI_Comm*) pyO_new;
  if(!PySequence_Check(pyO_indexes) || !PySequence_Check(pyO_edges)){
		error("MPI_Graph_create(MPI_Comm old, [...indexes],[...egdes],reorder,MPI_Comm graph)- ind, and edges should be sequences.");
	}
	int nnodes = PySequence_Size(pyO_indexes);
	int edge_size = PySequence_Size(pyO_edges);
	int buff_index0 = 0;
	int buff_index1 = 0;
	int* indexes = BufferStore::getBufferStore()->getFreeIntArr(buff_index0,nnodes);
	int* edges = BufferStore::getBufferStore()->getFreeIntArr(buff_index1,edge_size);
	for(int i = 0; i < nnodes; i++){
		indexes[i] = (int) PyInt_AsLong(PySequence_Fast_GET_ITEM(pyO_indexes, i));
		if(indexes[i] < 0){
			error("MPI_Graph_create(MPI_Comm old, [...indexes],[...egdes],reorder,MPI_Comm graph), [...inds] is not good.");
		}
	}
	for(int i = 0; i < edge_size; i++){
		edges[i] = (int) PyInt_AsLong(PySequence_Fast_GET_ITEM(pyO_edges, i));
		if(edges[i] < 0){
			error("MPI_Graph_create(MPI_Comm old, [...indexes],[...egdes],reorder,MPI_Comm graph), [...edges] is not good.");
		}
	}
	int res = ORBIT_MPI_Graph_create(pyComm_old->comm,nnodes,indexes,edges,reorder,&pyComm_new->comm);
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index0);
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index1);	
	return Py_BuildValue("i",res);
}	

//this function will return a tuple with (nnodes, nedges)
static PyObject* mpi_graphdims_get(PyObject *self, PyObject *args) {
	PyObject* pyO;
	if(!PyArg_ParseTuple(args,"O:mpi_graphdims_get",&pyO)){	
		error("MPI_Graphdims_get(MPI_Comm comm) needs 1 parameter.");
	}
	int nnodes, nedges;
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO;
	ORBIT_MPI_Graphdims_get(pyComm->comm, &nnodes, &nedges);
	return Py_BuildValue("(ii)",nnodes,nedges);
}	

//this function will return a tuple with (indexes[], edges[])
static PyObject* mpi_graph_get(PyObject *self, PyObject *args) {
	PyObject* pyO;
	if(!PyArg_ParseTuple(args,"O:mpi_graph_get",&pyO)){	
		error("MPI_Graph_get(MPI_Comm comm) needs 1 parameter.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO;	
	int nnodes, nedges;
	ORBIT_MPI_Graphdims_get(pyComm->comm, &nnodes, &nedges);
	int buff_index0 = 0;
	int buff_index1 = 0;	
	int* indexes = BufferStore::getBufferStore()->getFreeIntArr(buff_index0,nnodes);
	int* edges = BufferStore::getBufferStore()->getFreeIntArr(buff_index1,nedges);	
	ORBIT_MPI_Graph_get(pyComm->comm,nnodes,nedges,indexes,edges);
	PyObject* pyArrInd = PyTuple_New(nnodes);
	for(int i = 0; i < nnodes; i++){
		if(!PyTuple_SetItem(pyArrInd,i,Py_BuildValue("i",indexes[i]))){
			error(" MPI_Graph_get(MPI_Comm comm) cannot create a resulting tuple.");
		}
	}
	PyObject* pyArrEdg = PyTuple_New(nnodes);
	for(int i = 0; i < nedges; i++){
		if(!PyTuple_SetItem(pyArrEdg,i,Py_BuildValue("i",edges[i]))){
			error(" MPI_Graph_get(MPI_Comm comm) cannot create a resulting tuple.");
		}
	}
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index0);
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index1);	
	return Py_BuildValue("(OO)",pyArrInd,pyArrEdg);
}		
	
//this function will return a new rank of the calling process
static PyObject* mpi_graph_map(PyObject *self, PyObject *args) {
	PyObject* pyO_old; PyObject* pyO_indexes; PyObject* pyO_edges;
	if(!PyArg_ParseTuple(args,"OOOiO:mpi_graph_map",&pyO_old,&pyO_indexes,&pyO_edges)){	
		error("MPI_Graph_map(MPI_Comm old, [...indexes],[...egdes]) needs 4 parameters.");
	}
	pyORBIT_MPI_Comm* pyComm_old = (pyORBIT_MPI_Comm*) pyO_old;
  if(!PySequence_Check(pyO_indexes) || !PySequence_Check(pyO_edges)){
		error("MPI_Graph_map(MPI_Comm old, [...indexes],[...egdes])- indexes and edges should be sequences.");
	}
	int nnodes = PySequence_Size(pyO_indexes);
	int edge_size = PySequence_Size(pyO_edges);
	int buff_index0 = 0;
	int buff_index1 = 0;	
	int* indexes = BufferStore::getBufferStore()->getFreeIntArr(buff_index0,nnodes);
	int* edges = BufferStore::getBufferStore()->getFreeIntArr(buff_index1,edge_size);
	for(int i = 0; i < nnodes; i++){
		indexes[i] = (int) PyInt_AsLong(PySequence_Fast_GET_ITEM(pyO_indexes, i));
		if(indexes[i] < 0){
			error("MPI_Graph_map(MPI_Comm old, [...indexes],[...egdes]), - [...inds] is not good.");
		}
	}
	for(int i = 0; i < edge_size; i++){
		edges[i] = (int) PyInt_AsLong(PySequence_Fast_GET_ITEM(pyO_edges, i));
		if(edges[i] < 0){
			error("MPI_Graph_map(MPI_Comm old, [...indexes],[...egdes]), - [...edges] is not good.");
		}
	}
	int newrank = MPI_UNDEFINED;
	if(ORBIT_MPI_Graph_map(pyComm_old->comm,nnodes,indexes,edges,&newrank) != MPI_SUCCESS){
		error("MPI_Graph_map(MPI_Comm old, [...indexes],[...egdes]), fatal error STOP.");
	}
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index0);
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index1);	
	return Py_BuildValue("i",newrank);
}		

//this function will return a number of neighbors
static PyObject* mpi_graph_neighbors_count(PyObject *self, PyObject *args) {
	PyObject* pyO;
	int rank, nneighbors;
	if(!PyArg_ParseTuple(args,"Oi:mpi_graph_neighbors_count",&pyO,&rank)){	
		error("MPI_Graph_neighbors_count(MPI_Comm comm, rank) needs 2 parameters.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO;
	if(ORBIT_MPI_Graph_neighbors_count(pyComm->comm, rank, &nneighbors) != MPI_SUCCESS){
		error("MPI_Graph_neighbors_count(MPI_Comm comm, rank) fatal error. STOP.");
	}
	return Py_BuildValue("i",nneighbors);
}

//this function will return a tuple with ranks of neighbors
static PyObject* mpi_graph_neighbors(PyObject *self, PyObject *args) {
	PyObject* pyO;
	int rank, nneighbors;
	if(!PyArg_ParseTuple(args,"Oi:mpi_graph_neighbors",&pyO,&rank)){	
		error("MPI_Graph_neighbors(MPI_Comm comm, rank) needs 2 parameters.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO;
	ORBIT_MPI_Graph_neighbors_count(pyComm->comm, rank, &nneighbors);
	int buff_index = 0;
	int* neighbors = BufferStore::getBufferStore()->getFreeIntArr(buff_index,nneighbors);
	ORBIT_MPI_Graph_neighbors(pyComm->comm, rank, nneighbors, neighbors);
	PyObject* pyArr = PyTuple_New(nneighbors);
	for(int i = 0; i < nneighbors; i++){
		if(!PyTuple_SetItem(pyArr,i,Py_BuildValue("i",neighbors[i]))){
			error(" MPI_Graph_neighbors(MPI_Comm comm, rank) cannot create a resulting tuple.");
		}
	}	
	BufferStore::getBufferStore()->setUnusedIntArr(buff_index);
	return pyArr;
}

