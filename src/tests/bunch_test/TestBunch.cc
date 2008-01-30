///////////////////////////////////////////////////////////////////////////
//
// Orbit Bunch class test
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"

#include <string>
#include <iostream>

#include "Bunch.hh"
#include "ParticleAttributesFactory.hh"
#include "ParticleMacroSize.hh"

int main (int argc, char **argv)
{

  ORBIT_MPI_Init(&argc,&argv);

  int rank = 0;
  int size = 0;

  ORBIT_MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  ORBIT_MPI_Comm_size(MPI_COMM_WORLD,&size);

  if(rank == 0){std::cout<<"======START Bunch Class Test================"<<std::endl;}

  Bunch* bunch = new Bunch();

  std::cout<<"Bunch attr. =time="<< bunch->getBunchAttributes()->doubleVal("time") <<std::endl;

  ParticleAttributes* partAttr = ParticleAttributesFactory::getParticleAttributesInstance("macrosize",bunch);

  bunch->addParticleAttributes(partAttr);

  ParticleMacroSize* macroSizeAttr = (ParticleMacroSize*) bunch->getParticleAttributes("macrosize");

  int nParts = 10000;
  for(int i=0; i < nParts; i++){
    bunch->addParticle(i*100.+1+10*rank+0.001,
		       i*100.+2+10*rank+0.001,
		       i*100.+3+10*rank+0.001,
		       i*100.+4+10*rank+0.001,
		       i*100.+5+10*rank+0.001,
		       i*100.+6+10*rank+0.001);

    macroSizeAttr->marcosize(i) = 1.1+0.0022222;

		macroSizeAttr->attArr(i)[0] = 2.2 + 0.0033333;

    if(rank == 0 && i%1000 == 0){std::cout<<"add i = "<< i <<std::endl;}
  }

  for(int i=0; i < nParts; i+= 2){
    bunch->deleteParticleFast(i);
    if(rank == 0 && i%1000 == 0){std::cout<<"remove i = "<< i <<std::endl;}
  }

  std::cout<<"Bunch->size()="<< bunch->getSize() <<std::endl;

  bunch->compress();

  std::cout<<"Bunch->size()="<< bunch->getSize() <<std::endl;

  //if(rank == 0) std::cout<<"debug removeParticleAttributes(macrosize)"<<std::endl;
  //bunch->removeParticleAttributes("macrosize");


  if(rank == 0) std::cout<<"debug start writing file"<<std::endl;
  bunch->print("bunch_test.dat");


  if(rank == 0) std::cout<<"debug deleteAllParticles()"<<std::endl;
  bunch->deleteAllParticles();


  //bunch->addParticleAttributes(macroSizeAttr);


  std::vector<string> nms;
  bunch->readParticleAttributesNames("bunch_test.dat",nms);
  for(int i=0, n= nms.size(); i < n; i++){
    std::cout<<"==== readParticleAttributesNames i="<<i<<" name="<<nms[i]<<std::endl;
  }

  bunch->readBunchCoords("bunch_test.dat");

  bunch->print("bunch_test_new.dat");

  std::cout<<"rank="<< rank <<" Bunch->size()="<< bunch->getSize() <<std::endl;

  delete bunch;

  if(rank == 0){std::cout<<"======STOP  Bunch Class Test================"<<std::endl;}

  return 0;
}
