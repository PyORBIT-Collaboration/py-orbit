//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticlesWithIdFunctions.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    09/01/2015
//
// DESCRIPTION
//    A set of functions for bunches with the ParticleIdNumber attribute
//
///////////////////////////////////////////////////////////////////////////

#include <algorithm>    // std::sort
#include <vector>       // std::vector

#include "ParticlesWithIdFunctions.hh"
#include "ParticleIdNumber.hh"

#include "BufferStore.hh"
#include "MatrixOperations.hh"

#include "orbit_mpi.hh"

namespace OrbitUtils{
	
	struct comparator {
		ParticleIdNumber* partAttr;
		bool operator() (int i,int j) {
			int id_i = partAttr->getIdNumber(i);
			int id_j = partAttr->getIdNumber(j);
		  return (id_i<id_j);
		}
	} mycomparator;	
	
	/** A function that will sort bunch according to Id.*/
	void bunch_sort_id(Bunch* bunch){
		int n_parts = bunch->getSize();
		if(n_parts == 0) return;		
		int size_MPI,rank_MPI;
		ORBIT_MPI_Comm_size(bunch->getMPI_Comm_Local()->comm, &size_MPI);
		ORBIT_MPI_Comm_rank(bunch->getMPI_Comm_Local()->comm, &rank_MPI);
		if(bunch->hasParticleAttributes("ParticleIdNumber") == 0){
			if(rank_MPI == 0){
				std::cerr << "OrbitUtils::bunch_utils_functions::bunch_sort_id(Bunch* bunch) function"<< std::endl;
				std::cerr << "There is no ParticleAttributes in Bunch with this name."<< std::endl;
				std::cerr << "name:"<<" ParticleIdNumber "<< std::endl;
			}
			ORBIT_MPI_Finalize();		
		}
		bunch->compress();
		ParticleIdNumber* partAttr = (ParticleIdNumber*) bunch->getParticleAttributes("ParticleIdNumber");
		std::vector<int> numb_vector(n_parts);
		for(int i=0; i<n_parts; ++i){
			numb_vector[i] = i;
		}
		mycomparator.partAttr = partAttr;
		std::sort(numb_vector.begin(),numb_vector.end(),mycomparator);		
		std::vector<double*> partVect(n_parts);
		std::vector<int> idVect(n_parts);
		double** coordArr = bunch->coordArr();
		for(int i=0; i<n_parts; ++i){
			partVect[i] = coordArr[i];
			idVect[i] =partAttr->getIdNumber(i); 
		}
		int ind = 0;
		for(int i=0; i<n_parts; ++i){
			ind = numb_vector[i];
			coordArr[i] = partVect[ind];
			partAttr->setIdNumber(i,idVect[ind]);
		}
	}
	
	
	/** A function analyzes two bunches assuming the vectors of coordinates 
	    transformation x_out = A*x_in + b, where x_in the initial coordinates,
			and x_out final. The results are matrix A (7x7) with the last column as a vector b.
			We assume that bunch_in and bunch_out are sorted according to Id PartAttr, and
			bunch_out has less particles than bunch_in (some of them could be lost).
			It returns 0 if unsuccessful or the size of the statistics otherwise.
	*/
	int transport_mtrx(Bunch* bunch_in, Bunch* bunch_out, Matrix* A_mtr){
		
		int size_MPI,rank_MPI;
		ORBIT_MPI_Comm_size(bunch_in->getMPI_Comm_Local()->comm, &size_MPI);
		ORBIT_MPI_Comm_rank(bunch_in->getMPI_Comm_Local()->comm, &rank_MPI);		
		if(bunch_in->hasParticleAttributes("ParticleIdNumber") == 0 || bunch_out->hasParticleAttributes("ParticleIdNumber") == 0){
			if(rank_MPI == 0){
				std::cerr << "OrbitUtils::bunch_utils_functions::transport_mtrx(...) function"<< std::endl;
				std::cerr << "There is no ParticleAttributes  in Bunch with this name."<< std::endl;
				std::cerr << "name:"<<" ParticleIdNumber "<< std::endl;
			}
			ORBIT_MPI_Finalize();		
		}
		if(A_mtr->rows() != 7 || A_mtr->columns() != 7){
			if(rank_MPI == 0){
				std::cerr << "OrbitUtils::bunch_utils_functions::transport_mtrx(...) function"<< std::endl;
				std::cerr << "Matix have wrong size (not 7x7)!"<< std::endl;
				std::cerr << "Matrix A is:"<<A_mtr->rows()<<"x"<<A_mtr->columns()<< std::endl;
			}
			ORBIT_MPI_Finalize();				
		}		
		if(bunch_in->getMPI_Comm_Local() != bunch_out->getMPI_Comm_Local()){
			if(rank_MPI == 0){
				std::cerr << "OrbitUtils::bunch_utils_functions::transport_mtrx(...) function"<< std::endl;
				std::cerr << "Bunches In and Out have different MPI communicators!"<< std::endl;
				std::cerr << "That is WRONG!"<< std::endl;
			}
			ORBIT_MPI_Finalize();			
		}
		//---------------fill out the temporary bunches
		Bunch* b_in_tmp = new Bunch();
		Bunch* b_out_tmp = new Bunch();
		bunch_in->copyEmptyBunchTo(b_in_tmp);
		bunch_out->copyEmptyBunchTo(b_out_tmp);
		bunch_sort_id(bunch_in);
		bunch_sort_id(bunch_out);
		int n_parts_in = bunch_in->getSize();
		int n_parts_out = bunch_out->getSize();
		ParticleIdNumber* partAttr_in = (ParticleIdNumber*) bunch_in->getParticleAttributes("ParticleIdNumber");		
		ParticleIdNumber* partAttr_out = (ParticleIdNumber*) bunch_out->getParticleAttributes("ParticleIdNumber");
		int ind_start_in = 0;
		for(int ind_out = 0; ind_out < n_parts_out; ind_out++){
			int id_out =  partAttr_out->getIdNumber(ind_out);
			for(int ind_in = ind_start_in; ind_in < n_parts_in; ind_in++){
				if(partAttr_out->getIdNumber(ind_in) == id_out){
					double* arr = bunch_in->coordArr()[ind_in];
					b_in_tmp->addParticle(arr[0],arr[1],arr[2],arr[3],arr[4],arr[5]);
					arr = bunch_out->coordArr()[ind_out];
					b_out_tmp->addParticle(arr[0],arr[1],arr[2],arr[3],arr[4],arr[5]);
					ind_start_in = ind_in;
					break;
				}
			}
		}
		b_in_tmp->compress();
		b_out_tmp->compress();
		//-------------analysis of the temporary bunches
		A_mtr->zero();
		int n_parts =  b_in_tmp->getSize();
		int n_parts_global = b_in_tmp->getSizeGlobal();
		if(n_parts_global == 0) return 0;
		
		int buff_index0 = 0;
		int buff_index1 = 0;
		int buff_index2 = 0;
		int buff_index3 = 0;
		int buff_index4 = 0;
		int buff_index5 = 0;		
		double* arr_avg_in  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0,6);
		double* arr_avg_out = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1,6);
		double* arr_avg_in_mpi  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index2,6);
		double* arr_avg_out_mpi = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index3,6);	
		double* mtrx_arr  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index4,36);
		double* mtrx_arr_mpi = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index5,36);				
		
		for (int i = 0; i < 6; i++){
			arr_avg_in[i] = 0.; arr_avg_out[i] = 0.;
		}
		for(int ind = 0; ind < n_parts; ind++){
			for (int i = 0; i < 6; i++){
				arr_avg_in[i] += b_in_tmp->coordArr()[ind][i]; 
				arr_avg_out[i] += b_out_tmp->coordArr()[ind][i];
			}			
		}
		
		ORBIT_MPI_Allreduce(arr_avg_in,arr_avg_in_mpi,6,MPI_DOUBLE,MPI_SUM,b_in_tmp->getMPI_Comm_Local()->comm);
		ORBIT_MPI_Allreduce(arr_avg_out,arr_avg_out_mpi,6,MPI_DOUBLE,MPI_SUM,b_in_tmp->getMPI_Comm_Local()->comm);
		
		for (int i = 0; i < 6; i++){
			arr_avg_in_mpi[i] /= n_parts_global; 
			arr_avg_out_mpi[i] /= n_parts_global;
		}			
		//--- now we have avg values for coordinates in and out
		
		for(int ind = 0; ind < n_parts; ind++){
			for (int i = 0; i < 6; i++){
				b_in_tmp->coordArr()[ind][i] -= arr_avg_in_mpi[i]; 
				b_out_tmp->coordArr()[ind][i] -= arr_avg_out_mpi[i];
			}			
		}	
		
		// A = (X^T * X)^-1 * (X^T *Y)  start
		Matrix* XTXmtrx = new Matrix(6,6);
		Matrix* XTYmtrx = new Matrix(6,6);
		Matrix* Amtrx = new Matrix(6,6);
		XTXmtrx->zero();
		XTYmtrx->zero();
		Amtrx->zero();
		for (int i = 0; i < 6; i++){
			for (int j = i; j < 6; j++){
				for(int ind = 0; ind < n_parts; ind++){
					 XTXmtrx->getArray()[i][j] += b_in_tmp->coordArr()[ind][i]*b_in_tmp->coordArr()[ind][j];
				}
			}
		}
		for (int i = 0; i < 6; i++){
			for (int j = 0; j < i; j++){
				XTXmtrx->getArray()[i][j] = XTXmtrx->getArray()[j][i];
			}
		}
		for (int i = 0; i < 6; i++){
			for (int j = 0; j < 6; j++){
				for(int ind = 0; ind < n_parts; ind++){
					XTYmtrx->getArray()[i][j] += b_in_tmp->coordArr()[ind][i]*b_out_tmp->coordArr()[ind][j];
				}
			}
		}	
		
		//----- sum XTXmtrx and XTYmtrx over MPI
		int count = 0;
		for (int i = 0; i < 6; i++){
			for (int j = 0; j < 6; j++){
				mtrx_arr[count] = XTXmtrx->getArray()[i][j];
				count++;
			}
		}
		ORBIT_MPI_Allreduce(mtrx_arr,mtrx_arr_mpi,36,MPI_DOUBLE,MPI_SUM,b_in_tmp->getMPI_Comm_Local()->comm);
		count = 0;
		for (int i = 0; i < 6; i++){
			for (int j = 0; j < 6; j++){
				XTXmtrx->getArray()[i][j] = mtrx_arr_mpi[count];
				count++;
			}
		}		
		
		count = 0;
		for (int i = 0; i < 6; i++){
			for (int j = 0; j < 6; j++){
				mtrx_arr[count] = XTYmtrx->getArray()[i][j];
				count++;
			}
		}
		ORBIT_MPI_Allreduce(mtrx_arr,mtrx_arr_mpi,36,MPI_DOUBLE,MPI_SUM,b_in_tmp->getMPI_Comm_Local()->comm);
		count = 0;
		for (int i = 0; i < 6; i++){
			for (int j = 0; j < 6; j++){
				XTYmtrx->getArray()[i][j] = mtrx_arr_mpi[count];
				count++;
			}
		}		
		
		// A = (X^T * X)^-1 * (X^T *Y)
		MatrixOperations::invert(Amtrx);		
		XTXmtrx->copyTo(Amtrx);
		MatrixOperations::invert(Amtrx);
		
		Amtrx->mult(XTYmtrx);
		for (int i = 0; i < 6; i++){
			A_mtr->getArray()[i][6] = 0.;
			for (int j = 0; j < 6; j++){
				A_mtr->getArray()[i][j] = Amtrx->getArray()[i][j];
				A_mtr->getArray()[i][6] -= Amtrx->getArray()[i][j]*arr_avg_in_mpi[j];
			}
			A_mtr->getArray()[i][6] += arr_avg_out_mpi[i];
		}
		A_mtr->getArray()[6][6] = 1.0;
		
		delete b_in_tmp;
		delete b_out_tmp;
		delete XTXmtrx;
		delete XTYmtrx;
		delete Amtrx;
		OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
		OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);			
		OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index2);
		OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index3);			
		OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index4);
		OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index5);			
		
		return n_parts_global;
	}

}
