//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    Bunch.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    06/22/2005
//
// DESCRIPTION
//    Specification and inline functions for a container for macro particles.
//
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <string>
#include <set>
#include <map>
#include <vector>

#include "ParticleAttributes.hh"
#include "SyncPart.hh"

//from utils
#include "AttributesBucket.hh"
#include "CppPyWrapper.hh"

using namespace std;

#ifndef BUNCH_H
#define BUNCH_H
///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    Bunch
//
///////////////////////////////////////////////////////////////////////////

class  Bunch: public OrbitUtils::CppPyWrapper
{
public:
  //--------------------------------------
  //the public methods of the Bunch class
  //--------------------------------------

  Bunch();
  virtual ~Bunch();
  
  double& x(int index);
  double& y(int index);
  double& z(int index);

  double& px(int index);
  double& py(int index);
  double& pz(int index);
  double& dE(int index);

  double& xp(int index);
  double& yp(int index);

  //only flag == 0 means that particle is dead
  int& flag(int index);

  //returns the pointer to the 6D coordinates array
  //the values order is : x,px,y,py,z,pz
	//This can speed up operations with coordinates
  double* coordPartArr(int index);

  //returns the pointer to the [NumbOfPart][6D] coordinates array
  //the values order is : [index][x,px,y,py,z,pz]
	//This can speed up operations with coordinates
  double** coordArr();

	//wrap longitudinal coordinates assuming the certain ring length
	void ringwrap(double ring_length);

  //adds macro-particle
  //returns index of the new particle in the bunch
  int addParticle(double x, double px, double y, double py,
                  double z, double pz_or_dE);

  //removes a macro-particle from a bunch
  //you need to compress the bunch after one or +several delete operations
  void deleteParticleFast(int index);

  //removes a macro-particle from a bunch
  //You do not need to call compress method
  //The number of macro-particles will be changed inside this method
  void deleteParticle(int index);

  void compress();

  double getMass();                // GeV
  double getClassicalRadius();     // m
  double getCharge();              // sign and value in abs(e-charge) only
  double getMacroSize();

  double setMass(double mass);                // GeV
  double setClassicalRadius(double clR);      // m
  double setCharge(double chrg);              // sign and value in abs(e-charge) only
  double setMacroSize(double mcrsz);

  //returns the number of macro-particles in this CPU
  int getSize();

  //returns the number of macro-particles in all CPUs
  //it uses communications between CPUs over internal "Local" MPI communicator
  int getSizeGlobal();

  //returns the number of macro-particles in all CPUs
  //it uses latest results of communications in getSizeGlobal() method
	//it does not call any MPI functions, so it is fast
  int getSizeGlobalFromMemory();

  //returns total number of macro particles, alive and dead
  int getTotalCount();

  //return the capacity of the container
  int getCapacity();

  void print(std::ostream& Out);
  void print(const char* fileName);

  //these methods return the number of actual macro-particles that were read
  int readBunchCoords(const char* fileName, int nParts);
  int readBunchCoords(const char* fileName);
  int readParticleAttributesNames(const char* fileName,
		                              std::vector<std::string>& attr_names,
																	std::map<std::string,std::map<std::string,double> >& part_attr_dicts);
  void readParticleAttributes(const char* fileName);

  void deleteAllParticles();

  //methods related to the attribute buckets
  void addParticleAttributes(ParticleAttributes* attr);
  void addParticleAttributes(const std::string att_name,std::map<std::string,double> part_attr_dict);
  int  hasParticleAttributes(const std::string att_name);
  void removeParticleAttributes(const std::string name);
  void removeAllParticleAttributes();
  ParticleAttributes* getParticleAttributes(const std::string name);
  void getParticleAttributesNames(std::vector<std::string>& names);

  //this can be used for reading and writing the coordinates
  //the attributes will be initialized with default values
  void clearAllParticleAttributesAndMemorize();
  void restoreAllParticleAttributesFromMemory();

	//methods related to the sync. particle
	SyncPart* getSyncPart();

  //methods for the bunch attributes
  //user will get the reference and can use it
  OrbitUtils::AttributesBucket* getBunchAttributes();
  double getBunchAttributeDouble(const std::string att_name);
  int getBunchAttributeInt(const std::string att_name);
  void setBunchAttribute(const std::string att_name, double att_val);
  void setBunchAttribute(const std::string att_name, int att_val);
  void getIntBunchAttributeNames(std::vector<std::string>& names);
  void getDoubleBunchAttributeNames(std::vector<std::string>& names);
  void initBunchAttributes(const char* fileName);

	//copy methods

	//copy only bunch attributes and syncPart
	void copyEmptyBunchTo(Bunch* bunch);
	//copy all structure and macro-particles
	void copyBunchTo(Bunch* bunch);
	//copy particles and particles attributes
	void addParticlesTo(Bunch* bunch);

	//Parallel case
	pyORBIT_MPI_Comm* getMPI_Comm_Local();
	void setMPI_Comm_Local(pyORBIT_MPI_Comm* pyComm_Local);
	int getMPI_Size();
	int getMPI_Rank();
	
protected:

  //Initializes the different data that are the same for the all bunches.
  virtual void init();

private:
  //---------------------------------------
  //the private methods of the Bunch class
  //---------------------------------------

  friend class ParticleAttributes;

  void resize();
  void FinalizeExecution();
  void attrInit(int particle_index);

  //remove particle attributes without deleting it
  //It is used inside the memorize and restore particle attributes routines
  ParticleAttributes* removeParticleAttributesWithoutDelete(const std::string name);

  //this method provides access to the particles' attributes
  //array from ParticleAttributes class instance.
  //User is not supposed to use this method directly.
  double& getParticleAttributeVal(int ind, int attr_ind);

protected:

  double** arrAttr;
  double mass;
  double charge;
  double classicalRadius;

  //kinetic energy of the particle in GeV
  double energy;

  double macroSizeForAll;

  //---------------------------------------
  //the private members of the Bunch class
  //---------------------------------------

  int nDim;
  int nTotalSize;
  int nSize;
  int nNew;
  int nChunk;
  int nChunkMin;
  int sizeGlobal;

  int* arrFlag;
  double** arrCoord;

  //need of compress
  int needOfCompress;

  //----------------------------------------------
  //data members related to the ParticleAttributes
  //-----------------------------------------------
  std::map<std::string,ParticleAttributes*> attrCntrMap;
  std::map<std::string,int> attrCntrSizeMap;
  int attrCntrSize;
  int attributesSize;

  //inclusive Low and exclusive Upp indexes
  std::map<std::string,int> attrCntrLowIndMap;
  std::map<std::string,int> attrCntrUppIndMap;


  std::map<std::string,ParticleAttributes*> attrCntrMapTemp;
  std::vector<ParticleAttributes*> attrCntrVect;

  //bunch attributes
  OrbitUtils::AttributesBucket* bunchAttr;


	//synch. particle
	SyncPart* syncPart;

  //for MPI
  int iMPIini;
  int rank_MPI;
  int size_MPI;

	pyORBIT_MPI_Comm* pyComm_Local;
	
	//reference to the python wrapping class instance
	PyObject* py_wrapper;
	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
