//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticlesWithIdFunctions.cc
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
#include "ParticleMacroSize.hh"

#include "BufferStore.hh"
#include "MatrixOperations.hh"
#include "BunchTwissAnalysis.hh"

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
	
	
	int transport_mtrx(Bunch* bunch_in, Bunch* bunch_out, Matrix* A_mtr){
		return transport_mtrx(bunch_in,bunch_out,A_mtr,0,0,0);
	}
	
	
	/** A function analyzes two bunches assuming the vectors of coordinates 
	    transformation x_out = A*x_in + b, where x_in the initial coordinates,
			and x_out final. The results are matrix A (7x7) with the last column as a vector b.
			We assume that bunch_in and bunch_out are sorted according to Id PartAttr, and
			bunch_out has less particles than bunch_in (some of them could be lost).
			It returns 0 if unsuccessful or the size of the statistics otherwise.
	*/
	int transport_mtrx(Bunch* bunch_in, Bunch* bunch_out, Matrix* A_mtr, int appl_twiss_x, int appl_twiss_y, int appl_twiss_z){
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
		ParticleMacroSize* partMacroSizeAttr_in = NULL;
		ParticleMacroSize* partMacroSizeAttr_in_tmp = NULL;
		ParticleMacroSize* partMacroSizeAttr_out_tmp = NULL;
		if(bunch_in->hasParticleAttributes("macrosize") != 0 && bunch_out->hasParticleAttributes("macrosize") != 0){
			partMacroSizeAttr_in = (ParticleMacroSize*) bunch_in->getParticleAttributes("macrosize");
			partMacroSizeAttr_in_tmp  = (ParticleMacroSize*) b_in_tmp->getParticleAttributes("macrosize");
			partMacroSizeAttr_out_tmp  = (ParticleMacroSize*) b_out_tmp->getParticleAttributes("macrosize");
		}			
		int ind_start_in = 0;
		int count = 0;
		for(int ind_out = 0; ind_out < n_parts_out; ind_out++){
			int id_out =  partAttr_out->getIdNumber(ind_out);
			for(int ind_in = ind_start_in; ind_in < n_parts_in; ind_in++){
				if(partAttr_in->getIdNumber(ind_in) == id_out){
					double* arr = bunch_in->coordArr()[ind_in];
					b_in_tmp->addParticle(arr[0],arr[1],arr[2],arr[3],arr[4],arr[5]);
					arr = bunch_out->coordArr()[ind_out];
					b_out_tmp->addParticle(arr[0],arr[1],arr[2],arr[3],arr[4],arr[5]);
					ind_start_in = ind_in;
					if(partMacroSizeAttr_in != NULL){
						double m_size = partMacroSizeAttr_in->macrosize(ind_in);
						partMacroSizeAttr_in_tmp->macrosize(count) = m_size;
						partMacroSizeAttr_out_tmp->macrosize(count) = m_size;
					}		
					count++;
					break;
				}
			}
		}
		
		b_in_tmp->compress();
		b_out_tmp->compress();
		int n_parts =  b_in_tmp->getSize();
		int n_parts_global = b_in_tmp->getSizeGlobal();		
		double total_macrosize = 1.0*n_parts;
		
		//apply Twiss Gaussian weights to microsize. It will add macrosize Attr. if it does not exist
		apply_twiss_weghts(b_in_tmp, b_out_tmp,appl_twiss_x,appl_twiss_y,appl_twiss_z);
		
		partMacroSizeAttr_in_tmp = NULL;
		if(b_in_tmp->hasParticleAttributes("macrosize") != 0){
			partMacroSizeAttr_in_tmp  = (ParticleMacroSize*) b_in_tmp->getParticleAttributes("macrosize");
			//apply weights if the Macrosize particle attr. is defined 
			total_macrosize = 0.;
			for(int ind = 0; ind < n_parts; ind++){
				double m_size = partMacroSizeAttr_in_tmp->macrosize(ind);
				total_macrosize += m_size;
				for (int i = 0; i < 6; i++){
					b_in_tmp->coordArr()[ind][i] *= m_size; 
					b_out_tmp->coordArr()[ind][i] *= m_size;
				}			
			}
		}		

		//-------------analysis of the temporary bunches
		A_mtr->zero();
		if(n_parts_global > 6){
			
			int buff_index0 = -1;
			int buff_index1 = -1;
			int buff_index2 = -1;
			int buff_index3 = -1;
			int buff_index4 = -1;
			int buff_index5 = -1;		
			double* arr_avg_in  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0,6);
			double* arr_avg_out = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1,6);
			double* arr_avg_in_mpi  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index2,6);
			double* arr_avg_out_mpi = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index3,6);	
			double* mtrx_arr  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index4,36);
			double* mtrx_arr_mpi = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index5,36);	
			double total_macrosize_mpi = 0.;
			
			for (int i = 0; i < 6; i++){
				arr_avg_in[i] = 0.; arr_avg_out[i] = 0.;
			}
			for(int ind = 0; ind < n_parts; ind++){
				for (int i = 0; i < 6; i++){
					arr_avg_in[i] += b_in_tmp->coordArr()[ind][i]; 
					arr_avg_out[i] += b_out_tmp->coordArr()[ind][i];
				}			
			}
			
			ORBIT_MPI_Allreduce(&total_macrosize,&total_macrosize_mpi,1,MPI_DOUBLE,MPI_SUM,b_in_tmp->getMPI_Comm_Local()->comm);
			ORBIT_MPI_Allreduce(arr_avg_in,arr_avg_in_mpi,6,MPI_DOUBLE,MPI_SUM,b_in_tmp->getMPI_Comm_Local()->comm);
			ORBIT_MPI_Allreduce(arr_avg_out,arr_avg_out_mpi,6,MPI_DOUBLE,MPI_SUM,b_in_tmp->getMPI_Comm_Local()->comm);
			
			total_macrosize = total_macrosize_mpi;
			
			//std::cout<<"debug total_macrosize="<<total_macrosize<<std::endl;
			
			for (int i = 0; i < 6; i++){
				arr_avg_in_mpi[i] /= total_macrosize; 
				arr_avg_out_mpi[i] /= total_macrosize;
			}			
			//--- now we have avg values for coordinates in and out
			
			for(int ind = 0; ind < n_parts; ind++){
				double m_size = 1.0;
				if(partMacroSizeAttr_in != NULL) m_size = partMacroSizeAttr_in_tmp->macrosize(ind);
				for (int i = 0; i < 6; i++){
					b_in_tmp->coordArr()[ind][i] -= m_size*arr_avg_in_mpi[i]; 
					b_out_tmp->coordArr()[ind][i] -= m_size*arr_avg_out_mpi[i];
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
			
			delete XTXmtrx;
			delete XTYmtrx;
			delete Amtrx;
			OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
			OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);			
			OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index2);
			OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index3);			
			OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index4);
			OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index5);	
		}	
		
		delete b_in_tmp;
		delete b_out_tmp;		
		
		return n_parts_global;
	}
	
	/** A function analyzes two bunches assuming that they are already
	    sorted and synchronized according to the macro-particles Id. 
	    Coordinates of macro-particles in "in" and "out" bunches will be 
	    multiplied by the same numbers wx*wy*wz where
	    wx = exp(-(x^2+(alphax*x+betax*x')^2)/(2*(betax*emittancex))
	    etc.
	    Alpha, beta, emittance are the Twiss parameters for the corresponding 
	    plane.
	*/
	void apply_twiss_weghts(Bunch* bunch_in, Bunch* bunch_out,int appl_x,int appl_y,int appl_z){	
		if(appl_x == 0 && appl_y == 0 && appl_z == 0) return;
		
		int n_parts =  bunch_in->getSize();
		ParticleMacroSize* partMacroSizeAttr_in = NULL;
		if(bunch_in->hasParticleAttributes("macrosize") == 0){
			std::map<std::string,double> params_dict;
			bunch_in->addParticleAttributes("macrosize",params_dict);
			partMacroSizeAttr_in = (ParticleMacroSize*) bunch_in->getParticleAttributes("macrosize");
			for(int ind = 0; ind < n_parts; ind++){
				partMacroSizeAttr_in->macrosize(ind) = 1.;
			}		
			
		}
		ParticleMacroSize* partMacroSizeAttr_out = NULL;
		if(bunch_out->hasParticleAttributes("macrosize") == 0){
			std::map<std::string,double> params_dict;
			bunch_out->addParticleAttributes("macrosize",params_dict);
			partMacroSizeAttr_out = (ParticleMacroSize*) bunch_out->getParticleAttributes("macrosize");
			for(int ind = 0; ind < n_parts; ind++){
				partMacroSizeAttr_out->macrosize(ind) = 1.;
			}					
		}	
		
		partMacroSizeAttr_in = (ParticleMacroSize*) bunch_in->getParticleAttributes("macrosize");
		partMacroSizeAttr_out = (ParticleMacroSize*) bunch_out->getParticleAttributes("macrosize");

		BunchTwissAnalysis* twissAnalysis_in = new BunchTwissAnalysis();
		twissAnalysis_in->analyzeBunch(bunch_in);
		BunchTwissAnalysis* twissAnalysis_out = new BunchTwissAnalysis();
		twissAnalysis_out->analyzeBunch(bunch_out);	
		
		for(int ic = 0; ic < 3; ic++){
			if(ic == 0 && appl_x == 0) continue;
			if(ic == 1 && appl_y == 0) continue;
			if(ic == 2 && appl_z == 0) continue;
			double emitt_in = twissAnalysis_in->getEffectiveEmittance(ic);
			double gamma_in = twissAnalysis_in->getEffectiveGamma(ic);
			double alpha_in = twissAnalysis_in->getEffectiveAlpha(ic);
			double beta_in = twissAnalysis_in->getEffectiveBeta(ic);
			double emitt_out = twissAnalysis_out->getEffectiveEmittance(ic);
			double gamma_out = twissAnalysis_out->getEffectiveGamma(ic);
			double alpha_out = twissAnalysis_out->getEffectiveAlpha(ic);
			double beta_out = twissAnalysis_out->getEffectiveBeta(ic);			
			double w_in = 0.;
			double w_out = 0.;
			double x=0.;
			double xp = 0.;
			int ic2 = ic*2;
			int ic21 = ic*2+1;
			for(int ind = 0; ind < n_parts; ind++){
				x = bunch_in->coordArr()[ind][ic2];
				xp = bunch_in->coordArr()[ind][ic21];
				w_in = -(gamma_in*x*x+2*alpha_in*x*xp+beta_in*xp*xp)/(2*emitt_in);
				if(fabs(w_in) > 36.)  w_in = 0.;
				w_in = exp(w_in);
				w_out = -(gamma_out*x*x+2*alpha_out*x*xp+beta_out*xp*xp)/(2*emitt_out);
				if(fabs(w_out) > 36.)  w_out = 0.;
				w_out = exp(w_out);
				partMacroSizeAttr_out->macrosize(ind) *= w_in*w_out;
				partMacroSizeAttr_in->macrosize(ind) *= w_in*w_out;
			}
		}
		
		delete twissAnalysis_in;
		delete twissAnalysis_out;
		
	}
	
	

}
