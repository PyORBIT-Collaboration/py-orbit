//////////////////////////////// -*- C++ -*- /////////////////////////////
//
// FILE NAME
//   ParticleAttributesFactory.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    07/16/2005
//
// DESCRIPTION
//    A factory class for particle attributes classes.
//    Usually it will be used from a Bunch class instance.
//
///////////////////////////////////////////////////////////////////////////

#include "ParticleAttributesFactory.hh"

#include "ParticleMacroSize.hh"
#include "WaveFunctionAmplitudes.hh"
#include "AtomPopulations.hh"
#include "pq_coordinates.hh"
#include "part_time.hh"
#include "Evolution.hh"
#include "LostParticleAttributes.hh"
#include "ParticlePhaseAttributes.hh"
#include "ParticleIdNumber.hh"

ParticleAttributesFactory::ParticleAttributesFactory()
{
}

ParticleAttributesFactory::~ParticleAttributesFactory()
{
}

ParticleAttributes* ParticleAttributesFactory::getParticleAttributesInstance(
	const string name, 
	std::map<std::string,double> params_dict,
	Bunch* bunch)
{
	//for MPI --- start
	int rank_MPI = 0;
	int size_MPI = 1;
	int iMPIini  = 0;
	MPI_Comm MPI_COMM_Local = bunch->getMPI_Comm_Local()->comm;
	ORBIT_MPI_Initialized(&iMPIini);
	
	if(iMPIini > 0){
		ORBIT_MPI_Comm_size(MPI_COMM_Local, &size_MPI);
		ORBIT_MPI_Comm_rank(MPI_COMM_Local, &rank_MPI);
	}	
	//for MPI --- stop
	
	ParticleAttributes* part_atrs = NULL;
	if(name == "empty"){
		part_atrs = new ParticleAttributes(bunch,0);
	}
	
	if(name == "macrosize"){
		part_atrs = new ParticleMacroSize(bunch);
	}
	
	if(name == "ParticleIdNumber"){
		part_atrs = new ParticleIdNumber(bunch);
	}
  
	if(name == "LostParticleAttributes"){
		part_atrs = new LostParticleAttributes(bunch);
	}

	if(name == "ParticlePhaseAttributes"){
		part_atrs = new ParticlePhaseAttributes(bunch); 
	}
  
	if(name == "Amplitudes"){
		if(params_dict.size() == 0){
			cout<<"dictionary Amplitudes(dict) should be defined "<<"\n";
		} else {
			if(params_dict.count("size") == 1){
				part_atrs = new WaveFunctionAmplitudes(bunch,(int) params_dict["size"]);
			} else {
				if(rank_MPI == 0){
					std::cerr << "ParticleAttributesFactory::getParticleAttributesInstance(name,dict)"<< std::endl;
					std::cerr << "MPI Communicator="<< MPI_COMM_Local << std::endl;
					std::cerr << "MPI size="<< size_MPI << std::endl;
					std::cerr << "MPI rank="<< rank_MPI << std::endl;
					std::cerr << "attr. name:"<< name << std::endl;					
					std::cerr << "There is no <size> specification in the dict. "<< std::endl;
				}				
				ORBIT_MPI_Finalize("ParticleAttributesFactory::getParticleAttributesInstance. Stop.");
			}
		}
	}
  
  
	if(name == "Populations"){
		if(params_dict.size() == 0){
			cout<<"dictionary AtomPopulations(dict) should be defined "<<"\n";
		} else {
			if(params_dict.count("size") == 1){
				part_atrs = new AtomPopulations(bunch,(int) params_dict["size"]);
			} else {
				if(rank_MPI == 0){
					std::cerr << "ParticleAttributesFactory::getParticleAttributesInstance(name,dict)"<< std::endl;
					std::cerr << "MPI Communicator="<< MPI_COMM_Local << std::endl;
					std::cerr << "MPI size="<< size_MPI << std::endl;
					std::cerr << "MPI rank="<< rank_MPI << std::endl;
					std::cerr << "attr. name:"<< name << std::endl;					
					std::cerr << "There is no <size> specification in the dict. "<< std::endl;
				}				
				ORBIT_MPI_Finalize("ParticleAttributesFactory::getParticleAttributesInstance. Stop.");
			}
		}
	}
  
	if(name == "pq_coords"){
		if(params_dict.size() == 0){
			cout<<"dictionary pq_coords(dict) should be defined "<<"\n";
		} else {
			if(params_dict.count("size") == 1){
				part_atrs = new pq_coordinates(bunch,(int) params_dict["size"]);
			} else {
				if(rank_MPI == 0){
					std::cerr << "ParticleAttributesFactory::getParticleAttributesInstance(name,dict)"<< std::endl;
					std::cerr << "MPI Communicator="<< MPI_COMM_Local << std::endl;
					std::cerr << "MPI size="<< size_MPI << std::endl;
					std::cerr << "MPI rank="<< rank_MPI << std::endl;
					std::cerr << "attr. name:"<< name << std::endl;					
					std::cerr << "There is no <size> specification in the dict. "<< std::endl;
				}				
				ORBIT_MPI_Finalize("ParticleAttributesFactory::getParticleAttributesInstance. Stop.");
			}
		}
	}
  
	if(name == "part_time"){
		if(params_dict.size() == 0){
			cout<<"dictionary prf_time(dict) should be defined "<<"\n";
		} else {
			if(params_dict.count("size") == 1){
				part_atrs = new part_time(bunch, (int)params_dict["size"]);
			} else {
				if(rank_MPI == 0){
					std::cerr << "ParticleAttributesFactory::getParticleAttributesInstance(name,dict)"<< std::endl;
					std::cerr << "MPI Communicator="<< MPI_COMM_Local << std::endl;
					std::cerr << "MPI size="<< size_MPI << std::endl;
					std::cerr << "MPI rank="<< rank_MPI << std::endl;
					std::cerr << "attr. name:"<< name << std::endl;					
					std::cerr << "There is no <size> specification in the dict. "<< std::endl;
				}				
				ORBIT_MPI_Finalize("ParticleAttributesFactory::getParticleAttributesInstance. Stop.");
			}
		}
	}
  
	if(name == "Evolution"){
		if(params_dict.size() == 0){
			cout<<"dictionary Evolution(dict) should be defined "<<"\n";
		} else {
			if(params_dict.count("size") == 1){
				part_atrs = new Evolution(bunch, (int) params_dict["size"]);
			} else {
				if(rank_MPI == 0){
					std::cerr << "ParticleAttributesFactory::getParticleAttributesInstance(name,dict)"<< std::endl;
					std::cerr << "MPI Communicator="<< MPI_COMM_Local << std::endl;
					std::cerr << "MPI size="<< size_MPI << std::endl;
					std::cerr << "MPI rank="<< rank_MPI << std::endl;
					std::cerr << "attr. name:"<< name << std::endl;					
					std::cerr << "There is no <size> specification in the dict. "<< std::endl;
				}				
				ORBIT_MPI_Finalize("ParticleAttributesFactory::getParticleAttributesInstance. Stop.");
			}
		}
	}
  	
  
	if(part_atrs == NULL) {
		if(rank_MPI == 0){
			std::cerr << "ParticleAttributesFactory::getParticleAttributesInstance(const string name, Bunch* bunch)"<< std::endl;
			std::cerr << "MPI Communicator="<< MPI_COMM_Local << std::endl;
			std::cerr << "MPI size="<< size_MPI << std::endl;
			std::cerr << "MPI rank="<< rank_MPI << std::endl;
			std::cerr << "There is not a particle attirubutes class with such name in the Factory."<< std::endl;
			std::cerr << "attr. name:"<< name << std::endl;
		}
		ORBIT_MPI_Finalize("ParticleAttributesFactory::getParticleAttributesInstance. Stop.");
		return part_atrs;
	}
	
	//copy the particle attributes dictionary
	part_atrs->parameterDict = params_dict;
	part_atrs->parameterDict["size"] = part_atrs->getAttSize();
	
	return part_atrs;
}

void ParticleAttributesFactory::getParticleAttributesNames(std::vector<string>& names){
	names.clear();
	names.push_back("macrosize");
	names.push_back("ParticleIdNumber");
	names.push_back("Amplitudes");
	names.push_back("Populations");
	names.push_back("pq_coords");
	names.push_back("part_time");
	names.push_back("Evolution");
	names.push_back("LostParticleAttributes");
	names.push_back("ParticlePhaseAttributes");
}



