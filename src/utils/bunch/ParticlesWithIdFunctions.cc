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
	void bunch_sort_id(Bunch* bunch) {
		int size_MPI,rank_MPI;
		ORBIT_MPI_Comm_size(bunch->getMPI_Comm_Local()->comm, &size_MPI);
		ORBIT_MPI_Comm_rank(bunch->getMPI_Comm_Local()->comm, &rank_MPI);
		if(bunch->hasParticleAttributes("ParticleIdNumber") == 0){
			if(rank_MPI == 0){
				std::cerr << "OrbitUtils::bunch_sort_id(Bunch* bunch) function"<< std::endl;
				std::cerr << "There is no ParticleAttributes with this name."<< std::endl;
				std::cerr << "name:"<<" ParticleIdNumber "<< std::endl;
			}
			ORBIT_MPI_Finalize();		
		}
		bunch->compress();
		ParticleIdNumber* partAttr = (ParticleIdNumber*) bunch->getParticleAttributes("ParticleIdNumber");
		int n_parts = bunch->getSize();
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
	

}
