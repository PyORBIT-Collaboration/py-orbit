//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   InitialCoordsAttrFunctions.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    03/27/2017
//
// DESCRIPTION
//    A set of functions for bunches with the ParticleInitialCoordinates attributes.
//    At this moment there are following functions:
//
//    1. A function that will copy coordinates of the particles to 
//    the 6D initial coordinates Attribute (ParticleInitialCoordinates).
//    
//    void copyCoordsToInitCoordsAttr(Bunch* bunch);
//
//    2. Function that swaps init and existing 6D coordinates:
//
//    void swapInitCoordsAttrAndCoords(Bunch* bunch);
//
//    3. transport_mtrx_from_init_coords(bunch_in, Matrix* A_mtr)
//    4. transport_mtrx_from_init_coords(bunch_in, Matrix* A_mtr, int appl_twiss_x, int appl_twiss_y, int appl_twiss_z)
//    These functions will fill out the A_mtr that will be a transport matrix between
//    initial coordinates and 6D coordinates.
//    
//    5. This function will apply weights for macro-size of each particle according to the 
//    value wx = exp(-(x^2+(alphax*x+betax*x')^2)/(2*(betax*emittancex)) etc. This function
//    is used inside the function 4.
//
//    void apply_twiss_weights_for_init_coords(Bunch* bunch,int appl_x,int appl_y,int appl_z);	
//
//    Before using these functions it is recommended to remove particles far away from the center
//    of the phase space by using functions from TwissFilteringFunctions.cc
//    Function 3 assumes the equal weights for all macro-particles, and 4 uses weights according
//    wx = exp(-(x^2+(alphax*x+betax*x')^2)/(2*(betax*emittancex)) etc.
//
///////////////////////////////////////////////////////////////////////////

#include <algorithm>    // std::sort
#include <vector>       // std::vector

#include "InitialCoordsAttrFunctions.hh"
#include "ParticleInitialCoordinates.hh"
#include "ParticleMacroSize.hh"

#include "BufferStore.hh"
#include "MatrixOperations.hh"
#include "BunchTwissAnalysis.hh"

#include "orbit_mpi.hh"

namespace OrbitUtils{
	
	/** 
	  A function that will copy coordinates of the particles to 
	  the 6D initial coordinates Attribute (ParticleInitialCoordinates).
	  */	
	void copyCoordsToInitCoordsAttr(Bunch* bunch){
		if(bunch->hasParticleAttributes("ParticleInitialCoordinates") == 0){
			std::map<std::string,double> params_dict;
			bunch->addParticleAttributes("ParticleInitialCoordinates",params_dict);
		}
		ParticleInitialCoordinates* partAttr = (ParticleInitialCoordinates*) bunch->getParticleAttributes("ParticleInitialCoordinates");
		int n_parts = bunch->getSize();
		double** coordArr = bunch->coordArr();
		for(int i=0; i<n_parts; ++i){
			for(int j=0; j < 6; ++j){
				partAttr->attValue(i,j) = coordArr[i][j];
			}
		}
	}
	
	/** 
	  A function that will swap the initial coordinates Attribute 
	  (ParticleInitialCoordinates) and the 6D coordinates of the particles.
	  */
	void swapInitCoordsAttrAndCoords(Bunch* bunch){
		if(bunch->hasParticleAttributes("ParticleInitialCoordinates") == 0){
			std::map<std::string,double> params_dict;
			bunch->addParticleAttributes("ParticleInitialCoordinates",params_dict);
		}
		ParticleInitialCoordinates* partAttr = (ParticleInitialCoordinates*) bunch->getParticleAttributes("ParticleInitialCoordinates");
		int n_parts = bunch->getSize();
		double** coordArr = bunch->coordArr();
		double val = 0.;
		for(int i=0; i<n_parts; ++i){
			for(int j=0; j < 6; ++j){
				val = coordArr[i][j];
				coordArr[i][j] = partAttr->attValue(i,j);
				partAttr->attValue(i,j) = val;
			}
		}		
	}
		
	/** 
	  Calculates transport matrix A: x_out = A*x_in +b. A is a 7x7 matrix. 
		The last column of A is the b vector.
	  It returns 0 if unsuccessful or the size of the statistics otherwise.
	  The function uses the initial coordinates attributes as x_in and the usual
	  coordinates as x_out.
	*/	
	int transport_mtrx_from_init_coords(Bunch* bunch, Matrix* A_mtr){
		return transport_mtrx_from_init_coords(bunch,A_mtr,0,0,0);
	}
	
	/** A function analyzes the bunch with initail coords attributes assuming the vectors of coordinates 
	    transformation x_out = A*x_in + b, where x_in the initial coordinates,
			and x_out final. The results are matrix A (7x7) with the last column as a vector b.
			It returns 0 if unsuccessful or the size of the statistics otherwise.
	*/
	int transport_mtrx_from_init_coords(Bunch* bunch, Matrix* A_mtr, int appl_twiss_x, int appl_twiss_y, int appl_twiss_z){
		
		int size_MPI,rank_MPI;
		ORBIT_MPI_Comm_size(bunch->getMPI_Comm_Local()->comm, &size_MPI);
		ORBIT_MPI_Comm_rank(bunch->getMPI_Comm_Local()->comm, &rank_MPI);	
		
		if(bunch->hasParticleAttributes("ParticleInitialCoordinates") == 0){
			if(rank_MPI == 0){
				std::cerr << "OrbitUtils::bunch_utils_functions::transport_mtrx_from_init_coords(...) function"<< std::endl;
				std::cerr << "There is no ParticleAttributes  in Bunch with this name."<< std::endl;
				std::cerr << "Attr. name:"<<" ParticleInitialCoordinates "<< std::endl;
			}
			ORBIT_MPI_Finalize();		
		}
		
		if(A_mtr->rows() != 7 || A_mtr->columns() != 7){
			if(rank_MPI == 0){
				std::cerr << "OrbitUtils::bunch_utils_functions::transport_mtrx_init_coords(...) function"<< std::endl;
				std::cerr << "Matix have wrong size (not 7x7)!"<< std::endl;
				std::cerr << "Matrix A is:"<<A_mtr->rows()<<"x"<<A_mtr->columns()<< std::endl;
			}
			ORBIT_MPI_Finalize();				
		}		
		
		//---------------fill out the temporary bunches
		Bunch* b_tmp = new Bunch();
		bunch->copyBunchTo(b_tmp); 
		ParticleInitialCoordinates* partInitCoordsAttr = (ParticleInitialCoordinates*) b_tmp->getParticleAttributes("ParticleInitialCoordinates");
		
		int n_parts =  b_tmp->getSize();
		int n_parts_global = b_tmp->getSizeGlobal();		
		double total_macrosize = 1.0*n_parts;
		
		//apply Twiss Gaussian weights to microsize. It will add macrosize Attr. if it does not exist
		apply_twiss_weights_for_init_coords(b_tmp,appl_twiss_x,appl_twiss_y,appl_twiss_z);
		
		ParticleMacroSize* partMacroSizeAttr_tmp = NULL;
		if(b_tmp->hasParticleAttributes("macrosize") != 0){
			partMacroSizeAttr_tmp  = (ParticleMacroSize*) b_tmp->getParticleAttributes("macrosize");
			//apply weights if the Macrosize particle attr. is defined 
			total_macrosize = 0.;
			for(int ind = 0; ind < n_parts; ind++){
				double m_size = partMacroSizeAttr_tmp->macrosize(ind);
				total_macrosize += m_size;
				for (int i = 0; i < 6; i++){
					b_tmp->coordArr()[ind][i] *= m_size;
					partInitCoordsAttr->attValue(ind,i) *= m_size;
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
			double* arr_avg  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0,6);
			double* arr_avg_mpi  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1,6);
			double* arr_init_avg  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index2,6);
			double* arr_init_avg_mpi  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index3,6);			
			double* mtrx_arr  = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index4,36);
			double* mtrx_arr_mpi = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index5,36);	
			double total_macrosize_mpi = 0.;
			
			for (int i = 0; i < 6; i++){
				arr_avg[i] = 0.;
				arr_init_avg[i] = 0.;
			}
			for(int ind = 0; ind < n_parts; ind++){
				for (int i = 0; i < 6; i++){
					arr_avg[i] += b_tmp->coordArr()[ind][i];
					arr_init_avg[i] += partInitCoordsAttr->attValue(ind,i);
				}			
			}
			
			ORBIT_MPI_Allreduce(&total_macrosize,&total_macrosize_mpi,1,MPI_DOUBLE,MPI_SUM,b_tmp->getMPI_Comm_Local()->comm);
			ORBIT_MPI_Allreduce(arr_avg,arr_avg_mpi,6,MPI_DOUBLE,MPI_SUM,b_tmp->getMPI_Comm_Local()->comm);
			ORBIT_MPI_Allreduce(arr_init_avg,arr_init_avg_mpi,6,MPI_DOUBLE,MPI_SUM,b_tmp->getMPI_Comm_Local()->comm);
			
			total_macrosize = total_macrosize_mpi;
			
			//std::cout<<"debug total_macrosize="<<total_macrosize<<std::endl;
			
			for (int i = 0; i < 6; i++){
				arr_avg_mpi[i] /= total_macrosize; 
				arr_init_avg_mpi[i] /= total_macrosize;
			}			
			//--- now we have avg values for coordinates in and out
			
			for(int ind = 0; ind < n_parts; ind++){
				double m_size = 1.0;
				if(partMacroSizeAttr_tmp != NULL) m_size = partMacroSizeAttr_tmp->macrosize(ind);
				for (int i = 0; i < 6; i++){
					b_tmp->coordArr()[ind][i] -= m_size*arr_avg_mpi[i]; 
					partInitCoordsAttr->attValue(ind,i) -= m_size*arr_init_avg_mpi[i];
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
						double val = partInitCoordsAttr->attValue(ind,i)*partInitCoordsAttr->attValue(ind,j);
						XTXmtrx->getArray()[i][j] += val;
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
						XTYmtrx->getArray()[i][j] += partInitCoordsAttr->attValue(ind,i)*b_tmp->coordArr()[ind][j];
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
			ORBIT_MPI_Allreduce(mtrx_arr,mtrx_arr_mpi,36,MPI_DOUBLE,MPI_SUM,b_tmp->getMPI_Comm_Local()->comm);
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
			ORBIT_MPI_Allreduce(mtrx_arr,mtrx_arr_mpi,36,MPI_DOUBLE,MPI_SUM,b_tmp->getMPI_Comm_Local()->comm);
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
					A_mtr->getArray()[i][6] -= Amtrx->getArray()[i][j]*arr_init_avg_mpi[j];
				}
				A_mtr->getArray()[i][6] += arr_avg_mpi[i];
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
		
		delete b_tmp;
		
		return n_parts_global;
	}
	
	/** A function analyzes 6D coordinates and initial coords.
	    Macrosize of macro-particles in the bunch will be 
	    multiplied by wx*wy*wz where
	    wx = exp(-(x^2+(alphax*x+betax*x')^2)/(2*(betax*emittancex))
	    etc.
	    Alpha, beta, emittance are the Twiss parameters for the corresponding 
	    plane.
	*/
	void apply_twiss_weights_for_init_coords(Bunch* bunch,int appl_x,int appl_y,int appl_z){	
		if(appl_x == 0 && appl_y == 0 && appl_z == 0) return;
		
		int n_parts =  bunch->getSize();
		ParticleMacroSize* partMacroSizeAttr = NULL;
		if(bunch->hasParticleAttributes("macrosize") == 0){
			std::map<std::string,double> params_dict;
			bunch->addParticleAttributes("macrosize",params_dict);
			partMacroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
			for(int ind = 0; ind < n_parts; ind++){
				partMacroSizeAttr->macrosize(ind) = 1.;
			}		
		}

		partMacroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");
		BunchTwissAnalysis* twissAnalysis = new BunchTwissAnalysis();	
		
		for(int iswap = 0; iswap < 2; ++iswap){
	    if(iswap == 1){
	    	swapInitCoordsAttrAndCoords(bunch);
	    }

			twissAnalysis->analyzeBunch(bunch);
			
			for(int ic = 0; ic < 3; ic++){
				if(ic == 0 && appl_x == 0) continue;
				if(ic == 1 && appl_y == 0) continue;
				if(ic == 2 && appl_z == 0) continue;
				double emitt = twissAnalysis->getEffectiveEmittance(ic);
				double gamma = twissAnalysis->getEffectiveGamma(ic);
				double alpha = twissAnalysis->getEffectiveAlpha(ic);
				double beta = twissAnalysis->getEffectiveBeta(ic);		
				double w = 0.;
				double x=0.;
				double xp = 0.;
				int ic2 = ic*2;
				int ic21 = ic*2+1;
				for(int ind = 0; ind < n_parts; ind++){
					x = bunch->coordArr()[ind][ic2];
					xp = bunch->coordArr()[ind][ic21];
					w = -(gamma*x*x+2*alpha*x*xp+beta*xp*xp)/(2*emitt);
					if(fabs(w) > 36.)  w = 0.;
					w = exp(w);
					partMacroSizeAttr->macrosize(ind) *= w;
				}
			}
		}
		delete twissAnalysis;	
		swapInitCoordsAttrAndCoords(bunch);
	}
}
