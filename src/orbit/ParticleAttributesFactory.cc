//////////////////////////////// -*- C++ -*- //////////////////////////////
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

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "ParticleAttributesFactory.hh"

#include "ParticleMacroSize.hh"
///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

ParticleAttributesFactory::ParticleAttributesFactory()
{
}

ParticleAttributesFactory::~ParticleAttributesFactory()
{

}

ParticleAttributes* ParticleAttributesFactory::getParticleAttributesInstance(const string name, Bunch* bunch){
  ParticleAttributes* part_atrs = NULL;
  if(name == "empty"){
    part_atrs = new ParticleAttributes(bunch);
  }
  if(name == "macrosize"){
    part_atrs = new ParticleMacroSize(bunch);
  }
  
  if(part_atrs == NULL) {

    //for MPI
    int rank_MPI = 0;
    int size_MPI = 1;
    int iMPIini  = 0;
		MPI_Comm MPI_COMM_Local = bunch->getMPI_Comm_Local()->comm;
    ORBIT_MPI_Initialized(&iMPIini);

    if(iMPIini > 0){
      ORBIT_MPI_Comm_size(MPI_COMM_Local, &size_MPI);
      ORBIT_MPI_Comm_rank(MPI_COMM_Local, &rank_MPI);
    }

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

  return part_atrs;
}

void ParticleAttributesFactory::getParticleAttributesNames(std::vector<string>& names){
  names.clear();
  names.push_back("macrosize");
}



