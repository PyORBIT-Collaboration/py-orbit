//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    Bunch.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    06/22/2005
//
// DESCRIPTION
//    Source code for the class "Bunch" It is a container for macro particles.
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "Bunch.hh"

#include "ParticleAttributesFactory.hh"
#include "OrbitConst.hh"
#include "StringUtils.hh"
#include "BufferStore.hh"

#include <iomanip>
#include <string>

///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

Bunch::Bunch()
{

  //initialization all necessary variables and attributes
  //such as mass, charge, classical radius.
  //macro-size for all macro-particles
  init();

  //  x,px   y,py   z,pz
  nDim = 6;

  nChunkMin = 10;
  nChunk = nChunkMin;
  nTotalSize = nChunk;
  nNew = nSize = 0;
  sizeGlobal = 0;

  arrFlag = new int[nTotalSize];
  arrCoord = new double*[nTotalSize];
  for(int i=0; i < nTotalSize; i++){
    arrCoord[i] = new double[nDim];
  }

  //for MPI
	MPI_COMM_Local = MPI_COMM_WORLD;
  rank_MPI = 0;
  size_MPI = 1;
  iMPIini  = 0;
  ORBIT_MPI_Initialized(&iMPIini);

  if(iMPIini > 0){
    ORBIT_MPI_Comm_size(MPI_COMM_Local, &size_MPI);
    ORBIT_MPI_Comm_rank(MPI_COMM_Local, &rank_MPI);
  }

  //data members related to the ParticleAttributes
  attrCntrSize = 0;
  attributesSize = 0;
  arrAttr = NULL;

  //we do not need compress in the beginning
  needOfCompress = 0;
}

Bunch::~Bunch()
{
  for(int i=0; i < nTotalSize; i++){
    delete [] arrCoord[i];
  }
  delete [] arrFlag;
  delete [] arrCoord;

  if(attrCntrSize > 0 && attributesSize>0 ){
		for(int i=0; i < nTotalSize; i++){
			delete [] arrAttr[i];
		}
		delete [] arrAttr;
  }

  //delete controllers of particle attributes
  std::set<ParticleAttributes*>::iterator pos;
  for(pos = attrCntrSet.begin(); pos != attrCntrSet.end(); ++pos){
    delete  *pos;
  }

  //delete bunch attributes
  delete bunchAttr;

	//delete synchronous particle instance
	delete syncPart;
}

AttributesBucket* Bunch::getBunchAttributes(){
  return bunchAttr;
}

double Bunch::getBunchAttributeDouble(const std::string att_name){
  double att_val = 0.;
  if(bunchAttr->hasDoubleAttribute(att_name) > 0){
    return bunchAttr->doubleVal(att_name);
  }
  else{
    if(rank_MPI == 0){
      std::cerr << "The Bunch::getBunchAttributeDouble(const std::string att_name)"<< std::endl;
      std::cerr << "No attribute with name:"<< att_name<< std::endl;
      std::cerr << "Fatal error!"<< std::endl;
    }
    ORBIT_MPI_Finalize();
  }
  return att_val;
}

int Bunch::getBunchAttributeInt(const std::string att_name){
  int att_val = 0;
  if(bunchAttr->hasIntAttribute(att_name) > 0){
    return bunchAttr->intVal(att_name);
  }
  else{
    if(rank_MPI == 0){
      std::cerr << "The Bunch::getBunchAttributeInt(const std::string att_name)"<< std::endl;
      std::cerr << "No attribute with name:"<< att_name<< std::endl;
      std::cerr << "Fatal error!"<< std::endl;
    }
    ORBIT_MPI_Finalize();
  }
  return att_val;
}

void Bunch::setBunchAttribute(const std::string att_name, double att_val){
  bunchAttr->doubleVal(att_name) = att_val;
  if(att_name == "mass"){
    mass = att_val;
  }
  if(att_name == "charge"){
    charge = att_val;
  }
  if(att_name == "classical_radius"){
    classicalRadius = att_val;
  }
  if(att_name == "macro_size"){
    macroSizeForAll = att_val;
  }
}

void Bunch::setBunchAttribute(const std::string att_name, int att_val){
  bunchAttr->intVal(att_name) = att_val;
}

void Bunch::getIntBunchAttributeNames(std::vector<std::string>& names){
  bunchAttr->getIntAttributeNames(names);
}

void Bunch::getDoubleBunchAttributeNames(std::vector<std::string>& names){
  bunchAttr->getDoubleAttributeNames(names);
}

void Bunch::initBunchAttributes(const char* fileName){

  std::vector<std::string> attr_names;
  attr_names.clear();

  ifstream is;

  int error_ind = 0;
  if(rank_MPI == 0){
    is.open(fileName, std::ios::in);
    if (is.bad()){
      std::cerr << "The Bunch::initBunchAttributes(char* fileName)"<< std::endl;
      std::cerr << "Can not open file:"<< fileName<< std::endl;
      error_ind = 1;
    }
  }

  if(size_MPI > 1){
    ORBIT_MPI_Bcast (&error_ind,1, MPI_INT,    0, MPI_COMM_Local );
  }

  if(error_ind > 0){
    ORBIT_MPI_Finalize("The Bunch::initBunchAttributes. Stop.");
  }

  std::string  str;
  std::vector<string> v_str;

  int stop_ind = 0;
  int def_found_ind = 0;

  while(stop_ind == 0){

    def_found_ind = 0;

    if(rank_MPI == 0){
      if(!is.eof()){
        getline(is,str);
        if(strlen(str.c_str()) > 0  && str.c_str()[0] == '%'){
          int nT = StringUtils::Tokenize(str,v_str);
          if(nT > 3 && v_str[1] == "BUNCH_ATTRIBUTE_INT"){
            def_found_ind = 1;
          }
          if(nT > 3 && v_str[1] == "BUNCH_ATTRIBUTE_DOUBLE"){
            def_found_ind = 1;
          }
        }
        else{
          stop_ind = 1;
        }
      }
      else{
        stop_ind = 1;
      }
    }


    if(size_MPI > 1){
      ORBIT_MPI_Bcast (&stop_ind,1, MPI_INT,    0, MPI_COMM_Local );
      ORBIT_MPI_Bcast (&def_found_ind,1, MPI_INT,    0, MPI_COMM_Local );
    }

    if(stop_ind == 0 && def_found_ind == 1){
      if(size_MPI > 1){
        int strLength = strlen(str.c_str());
        ORBIT_MPI_Bcast ( &strLength,1, MPI_INT,    0, MPI_COMM_Local );
        char* char_tmp = BufferStore::getBufferStore()->getCharArr(0,strLength+1);
        strcpy(char_tmp, str.c_str());
        ORBIT_MPI_Bcast ( char_tmp,  strLength+1, MPI_CHAR,    0, MPI_COMM_Local );
        std::string str_new(char_tmp);
        StringUtils::Tokenize(str_new,v_str);
      }
      if(v_str.size() > 3 && v_str[1] == "BUNCH_ATTRIBUTE_INT"){
        int val = 0;
        sscanf( v_str[3].c_str(),"%d",&val);
        setBunchAttribute(v_str[2],val);
      }
      if(v_str.size() > 3 && v_str[1] == "BUNCH_ATTRIBUTE_DOUBLE"){
        double val = 0.;
        sscanf( v_str[3].c_str(),"%lf",&val);
        setBunchAttribute(v_str[2],val);
      }
    }

  }

  if(rank_MPI == 0){
    is.close();
  }

	//set synchronous particle parameters from file
	syncPart->readSyncPart(fileName);

}

void Bunch::copyEmptyBunchTo(Bunch* bunch){
	bunch->setMPI_Comm_Local(this->getMPI_Comm_Local());

	bunch->deleteAllParticles();

	//copy bunch attributes
	//the Attribute backet direct copy cannot be used
	//because mass,charge, class. radius, macro-size_for_all are
	//located not only in bunch attribute, but in the class Bunch
	//data also

	std::vector<std::string> names;

	this->getIntBunchAttributeNames(names);
	for(int i = 0, n = names.size(); i < n; i++){
		bunch->setBunchAttribute(names[i],this->getBunchAttributeInt(names[i]));
	}

	names.clear();

	this->getDoubleBunchAttributeNames(names);
	for(int i = 0, n = names.size(); i < n; i++){
		bunch->setBunchAttribute(names[i],this->getBunchAttributeDouble(names[i]));
	}

	//copy syncPart
	SyncPart* source = this->getSyncPart();
	SyncPart* target = bunch->getSyncPart();
	target->setTime(source->getTime());
	target->setFrequency(source->getFrequency());
	target->setXYZ(source->getX(),source->getY(),source->getZ());
	target->setPXYZ(source->getPX(),source->getPY(),source->getPZ());


	//copy particles attributes
	bunch->removeAllParticleAttributes();

	names.clear();

	this->getParticleAttributesNames(names);
	for(int i = 0, n = names.size(); i < n; i++){
		bunch->addParticleAttributes(names[i]);
	}
}

void Bunch::copyBunchTo(Bunch* bunch){
	this->copyEmptyBunchTo(bunch);
	this->addParticlesTo(bunch);
}

void Bunch::addParticlesTo(Bunch* bunch){
	for(int i = 0, n = this->getSize(); i < n; i++){
		//add particles coordinates
		bunch->addParticle(this->x(i),this->y(i),this->z(i),
			this->px(i),this->py(i),this->pz(i));
	}

	//add particles attributes
	std::vector<std::string> names_source;
	std::vector<std::string> names_target;
	this->getParticleAttributesNames(names_source);
	bunch->getParticleAttributesNames(names_target);

	if(names_source.size() != names_target.size()){
		if(rank_MPI == 0){
			std::cerr << "ParticleAttributes* Bunch::addParticlesTo(const Bunch* bunch)"<< std::endl;
			std::cerr << "The particles attributes have different structures for two bunches."<< std::endl;
			std::cerr << "===========   Target bunch:"<< std::endl;
			for(int i = 0, in = names_target.size(); i < in; i++){
				std::cerr << "Particle Attrubute:"<< names_target[i] << std::endl;
			}
			std::cerr << "===========   Source bunch:"<< std::endl;
			for(int i = 0, in = names_source.size(); i < in; i++){
				std::cerr << "Particle Attrubute:"<< names_source[i] << std::endl;
			}
		}
	}

	for(int j = 0, jn = names_source.size(); j < jn; j++){
		int found = 0;
		for(int i = 0, in = names_target.size(); i < in; i++){
			if(names_source[j] == names_target[i] ) found = 1;
		}
		if(found != 1){
			if(rank_MPI == 0){
				std::cerr << "ParticleAttributes* Bunch::addParticlesTo(const Bunch* bunch)"<< std::endl;
				std::cerr << "The particles attributes have different structures for two bunches."<< std::endl;
				std::cerr << "===========   Target bunch:"<< std::endl;
				for(int i = 0, in = names_target.size(); i < in; i++){
					std::cerr << "Particle Attrubute:"<< names_target[i] << std::endl;
				}
				std::cerr << "===========   Source bunch:"<< std::endl;
				for(int i = 0, in = names_source.size(); i < in; i++){
					std::cerr << "Particle Attrubute:"<< names_source[i] << std::endl;
				}
			}
			ORBIT_MPI_Finalize();
		}
	}

	for(int j = 0, jn = names_source.size(); j < jn; j++){
		ParticleAttributes* att_source = this->getParticleAttributes(names_source[j]);
		ParticleAttributes* att_target = bunch->getParticleAttributes(names_source[j]);
		int attSize = att_source->getAttSize();
		//loop over particles
		for(int i = 0, in = this->getSize(); i < in; i++){
			//loop over attribute values
			for(int k = 0, kn = attSize; k < kn; k++){
				att_target->attValue(i,k) = att_source->attValue(i,k);
			}
		}
	}
}

void Bunch::init(){
  bunchAttr = new AttributesBucket();

  //mass in [MeV/c^2]
  //charge of in the charge of electron
  //classical radius in m
  mass            = OrbitConst::mass_proton;
  charge          = OrbitConst::charge_proton;
  classicalRadius = OrbitConst::classicalRadius_proton;

  macroSizeForAll = 0.;

  bunchAttr->doubleVal("mass") = mass;
  bunchAttr->doubleVal("charge") = charge;
  bunchAttr->doubleVal("classical_radius") = classicalRadius;
  bunchAttr->doubleVal("macro_size") = macroSizeForAll;

	syncPart = new SyncPart(this);
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   Bunch::macroSize, Bunch::x, Bunch::y, Bunch::z,
//   Bunch::px, Bunch::py, Bunch::pz, Bunch::pp, Bunch::pt,
//   Bunch::flag
//
// DESCRIPTION
//   Each returns an element of the given index-th macroparticle
//
//   getMacroSize:  returns MacroSize of each macroparticle
//   x: y: z:       returns coordinates of each macroparticle
//   px: py: pz:    returns momentum of non-drifting macroparticle
//   pp: pt:        returns momentum of drifting  macroparticle
//   flag:          returns the status of macroparticle; dead or alive
//                  flag(index)=0 means index-th macroparticle is dead
//                  flag(index)=1 means index-th macroparticle is alive
//
// RETURNS
//   arrCoord[index][number of corresponding element]
//
///////////////////////////////////////////////////////////////////////////
double& Bunch::x(int index){    return arrCoord[index][0];}
double& Bunch::px(int index){   return arrCoord[index][1];}

double& Bunch::y(int index){    return arrCoord[index][2];}
double& Bunch::py(int index){   return arrCoord[index][3];}

double& Bunch::z(int index){    return arrCoord[index][4];}
double& Bunch::pz(int index){   return arrCoord[index][5];}

double& Bunch::phi(int index){  return arrCoord[index][4];}
double& Bunch::dE(int index){   return arrCoord[index][5];}

double& Bunch::xp(int index){   return arrCoord[index][1];}
double& Bunch::yp(int index){   return arrCoord[index][3];}

int & Bunch::flag(int index){   return arrFlag[index];}

double* Bunch::coordPartArr(int index){ return arrCoord[index];}
double** Bunch::coordArr(){ return arrCoord;}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   Bunch::resize
//
// DESCRIPTION
//   Expands the size of arrFlag,arrCoord,
//   when the number of macroparticles reaches nTotalSize. It should be
//   used in Bunch::addParticle.
//
// RETURNS
//   Nothing.
//
///////////////////////////////////////////////////////////////////////////

void Bunch::resize()
{

  if (nNew <= (nTotalSize - nChunk/3) && nNew >= (nTotalSize - nChunk)) return;

  //chunk should be big enough to avoid frequently changing size
  nChunk = (int) (nNew*0.2);
  if(nChunk < nChunkMin) nChunk = nChunkMin;

  int nOldTotalSize = nTotalSize;
  nTotalSize = (((int)(nNew/nChunk)) + 1)*nChunk;

  int* tmp_arrFlag     = new int    [nTotalSize];
  double** tmp_arrCoord    = new double*[nTotalSize];

  if(nOldTotalSize <= nTotalSize){
    for(int i = nOldTotalSize; i < nTotalSize; i++){
      tmp_arrCoord[i]    = new double[nDim];
    }

    for(int i=0; i < nOldTotalSize; i++){
      tmp_arrFlag[i]     = arrFlag[i];
      tmp_arrCoord[i]    = arrCoord[i];
    }
  }
  else{
    for(int i=0; i < nTotalSize; i++){
      tmp_arrFlag[i] = arrFlag[i];
      tmp_arrCoord[i]    = arrCoord[i];
    }
    for(int i = nTotalSize; i < nOldTotalSize; i++){
      delete [] arrCoord[i];
    }
  }

  delete [] arrFlag;
  delete [] arrCoord;

  arrFlag     = tmp_arrFlag;
  arrCoord    = tmp_arrCoord;

  //attributes resize
  if(attrCntrSize > 0 && attributesSize > 0){

    double** tmp_arrAttr    = new double*[nTotalSize];
    if(nOldTotalSize <= nTotalSize){
      for(int i = nOldTotalSize; i < nTotalSize; i++){
				tmp_arrAttr[i]    = new double[nDim];
      }
      for(int i=0; i < nOldTotalSize; i++){
				tmp_arrAttr[i]    = arrAttr[i];
      }
      delete [] arrAttr;
      arrAttr   = tmp_arrAttr;

      std::map<std::string,ParticleAttributes*>::iterator pos;
      for (pos = attrCntrMap.begin(); pos != attrCntrMap.end(); ++pos) {
				string name = pos->first;
				ParticleAttributes* attrCntrl = pos->second;
				for(int i = nOldTotalSize; i < nTotalSize; i++){
					attrCntrl->init(i);
				}
      }

    }
    else{
      for(int i=0; i < nTotalSize; i++){
				tmp_arrAttr[i]    = arrAttr[i];
      }
      for(int i = nTotalSize; i < nOldTotalSize; i++){
				delete [] arrAttr[i];
      }

      delete [] arrAttr;
      arrAttr   = tmp_arrAttr;
    }
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    Bunch::addParticle
//
// DESCRIPTION
//    adds a macro_particle with unique macro-size.
//    It is needed to call Bunch::compress to activate this new
//    macro-particle in calculation
//
// RETURNS
//    Nothing.
//
///////////////////////////////////////////////////////////////////////////

int Bunch::addParticle(double x, double px,
	double y, double py,
	double z_phi, double pz_dE)
{
  int n = nNew;
  resize();

  arrCoord[n][0] = x;
  arrCoord[n][1] = px;

  arrCoord[n][2] = y;
  arrCoord[n][3] = py;

  arrCoord[n][4] = z_phi;
  arrCoord[n][5] = pz_dE;

  arrFlag[n] = 1; //alive

  //attributes initialization
  attrInit(n);

  nNew = nNew + 1;

  if(needOfCompress == 0){
    nSize = nNew;
  }
  return n;
}

void Bunch::attrInit(int particle_index){
  if(attrCntrSize > 0 && attributesSize > 0){
    std::map<std::string,ParticleAttributes*>::iterator pos;
    for (pos = attrCntrMap.begin(); pos != attrCntrMap.end(); ++pos) {
      ParticleAttributes* attrCntrl = pos->second;
      attrCntrl->init(particle_index);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    Bunch::deleteParticleFast
//
// DESCRIPTION
//    deletes a macro_particle with declaring its flag is dead
//    It is needed to call Bunch::compress to remove this dead
//    macroparticle from calculation
//
// RETURNS
//    Nothing.
//
///////////////////////////////////////////////////////////////////////////

void Bunch::deleteParticleFast(int index)
{
  arrFlag[index] = 0; //dead

  //we need compress in the future
  needOfCompress = 1;
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    Bunch::deleteParticle
//
// DESCRIPTION
//    deletes a macro_particle and call Bunch::compress inside,
//    so the number of macro-particles changes after this method
//    This method is slower than deleteParticleFast
//
// RETURNS
//    Nothing.
//
///////////////////////////////////////////////////////////////////////////

void Bunch::deleteParticle(int index)
{
  deleteParticleFast(index);
  compress();
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    Bunch::compress
//
// DESCRIPTION
//    compresses the workspace to the number of "alive" macroparticles
//    with removing "dead" flagged macroparticles and renumbering "alive"
//    flagged macroparticles
//
// RETURNS
//    Nothing.
//
///////////////////////////////////////////////////////////////////////////

void Bunch::compress()
{

  if(needOfCompress == 0) return;

  int lowInd = 0;
  int uppInd = lowInd;
  int count = 0;

  double* tmp;
  int tmp_flag = 0;

  int lowIndChanged = 0;

  while( uppInd < nNew){

    lowIndChanged = 0;

    while( arrFlag[lowInd] != 0 && lowInd < nNew){
      count++;
      lowInd++;
      lowIndChanged = 1;
    }
    if(lowInd == (nNew -1)) break;

    if(lowIndChanged > 0){
      uppInd = lowInd + 1;
    }
    else{
      uppInd++;
    }

    while( arrFlag[uppInd] == 0 && uppInd < nNew){
      uppInd++;
    }

    if(uppInd < nNew){
      count++;
      tmp = arrCoord[lowInd];
      arrCoord[lowInd] = arrCoord[uppInd];
      arrCoord[uppInd] = tmp;
      tmp_flag = arrFlag[lowInd];
      arrFlag[lowInd] = arrFlag[uppInd];
      arrFlag[uppInd] = tmp_flag;
      if(attributesSize > 0){
        tmp = arrAttr[lowInd];
        arrAttr[lowInd] = arrAttr[uppInd];
        arrAttr[uppInd] = tmp;
      }
      lowInd++;
    }
  }
  nSize=count;
  nNew=count;

  //compression is done
  needOfCompress = 0;

  resize();
}

///////////////////////////////////////////////////////////////////////////
//
// getMass - mass of a particle in MeV
//
///////////////////////////////////////////////////////////////////////////

double Bunch::getMass(){   return mass;}
double Bunch::getCharge(){ return charge;}
double Bunch::getClassicalRadius(){ return classicalRadius;}
double  Bunch::getMacroSize(){ return macroSizeForAll;}

double Bunch::setMass(double val){
  mass = val;
  bunchAttr->doubleVal("mass") = val;
  return mass;
}

double Bunch::setCharge(double val){
  charge = val;
  bunchAttr->doubleVal("charge") = val;
  return charge;
}

double Bunch::setClassicalRadius(double val){
  classicalRadius = val;
  bunchAttr->doubleVal("classical_radius") = val;
  return classicalRadius;
}

double  Bunch::setMacroSize(double val){
  macroSizeForAll = val;
  bunchAttr->doubleVal("macro_size") = val;
  return macroSizeForAll;
}
///////////////////////////////////////////////////////////////////////////
//
// NAME
//    Bunch::getSize
//
// DESCRIPTION
//    returns the number of macroparticles
//    In order to reflect the results of the latest deleteParticle(),
//    you have to use this right after compress()
//
// RETURNS
//    the number of alive macro_particles
//
///////////////////////////////////////////////////////////////////////////

int Bunch::getSize(){  return nSize;}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    Bunch::getSizeGlobal
//
// DESCRIPTION
//    returns the number of macroparticles
//    In order to reflect the results of the latest deleteParticle(),
//    you have to use this right after compress()
//
// RETURNS
//    the number of alive macro_particles for all CPUs
//
///////////////////////////////////////////////////////////////////////////

int Bunch::getSizeGlobal()
{
  if(size_MPI == 1) {
    sizeGlobal = nSize;
    return sizeGlobal;
  }
  else{
    ORBIT_MPI_Allreduce(&nSize,&sizeGlobal,1,
		  MPI_INT,MPI_SUM,MPI_COMM_Local);
    return sizeGlobal;
  }
}


//returns the number of alive macro_particles for all CPUs from memory
//after Bunch::getSizeGlobal() called
int Bunch::getSizeGlobalFromMemory()
{
  return sizeGlobal;
}


//returns total number of macro particles, alive and dead
int Bunch::getTotalCount()
{
  return nNew;
}

//return the capacity of the container
int Bunch::getCapacity(){
  return nTotalSize;
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    Bunch::print
//
// DESCRIPTION
//    print every coordinates' components of macroparticles
//    and all attributes
//
///////////////////////////////////////////////////////////////////////////

void Bunch::print(std::ostream& Out)
{
  int N_flush = 10000;

  //single CPU case
  if(rank_MPI == 0){

    //print particle attributes controllers names
    Out << "% PARTICLE_ATTRIBUTES_CONTROLLERS_NAMES ";
    std::map<std::string,ParticleAttributes*>::iterator pos;
    for (pos = attrCntrMap.begin(); pos != attrCntrMap.end(); ++pos) {
      ParticleAttributes* attrCntr = pos->second;
      Out << attrCntr->name()<<" ";
    }
    Out << std::endl;

    //print bunch attributes
    std::vector<std::string> bunch_attr_names;
    bunchAttr->getIntAttributeNames(bunch_attr_names);
    for(int i = 0, n = bunch_attr_names.size(); i < n; i++){
      Out << "% BUNCH_ATTRIBUTE_INT "<<bunch_attr_names[i]<<"   ";
      Out << bunchAttr->intVal(bunch_attr_names[i]) <<" "<< std::endl;
    }
    bunch_attr_names.clear();
    bunchAttr->getDoubleAttributeNames(bunch_attr_names);
    for(int i = 0, n = bunch_attr_names.size(); i < n; i++){
      Out << "% BUNCH_ATTRIBUTE_DOUBLE "<<bunch_attr_names[i]<<"   ";
      Out << bunchAttr->doubleVal(bunch_attr_names[i]) <<" "<< std::endl;
    }

		//print synchronous particle parameters to the stream
	  syncPart->print(Out);

    Out << "% x[m] px[rad] y[m] py[rad] (z_or_phi or pz_or_dE) ";

    for (pos = attrCntrMap.begin(); pos != attrCntrMap.end(); ++pos) {
      ParticleAttributes* attrCntr = pos->second;
      Out << attrCntr->attrDescription()<<" ";
    }

    Out << std::endl;

    Out <<std::setprecision(8); //<< std::setiosflags(ios::scientific);

    for(int i = 0; i < nSize; i++){

      if(flag(i) > 0){

        Out << x(i)  <<" ";
        Out << px(i) <<" ";

        Out << y(i)  <<" ";
        Out << py(i) <<" ";

        Out << z(i)  <<" ";
        Out << pz(i) <<" ";

        for(int j = 0; j < attributesSize; j++){
          Out << getParticleAttributeVal(i,j)  <<" ";
        }

        Out <<std::endl;
      }

      if (i % N_flush == 0){Out.flush();}
    }
    Out.flush();
    if(size_MPI == 1) return;
  }


  //parallel case                                   ===== MPI start =====
  MPI_Status statusMPI ;
  int* nSizeArr     = BufferStore::getBufferStore()->getIntArr(0,size_MPI);
  int* nSizeArr_MPI = BufferStore::getBufferStore()->getIntArr(1,size_MPI);

  for(int i=0; i<size_MPI;i++){
    nSizeArr[i]=0;
    if(i==rank_MPI){nSizeArr[i]=nSize;}
  }
  ORBIT_MPI_Allreduce(nSizeArr,nSizeArr_MPI,size_MPI,
		MPI_INT,MPI_SUM,MPI_COMM_Local);

  //at this point all CPUs know about number of macro-particles on each CPU

  //max number of macro-particles that can be sent by CPU at once
  int nSizeChank = 1000;

  int nDimAndAttr = nDim+attributesSize+1;


  double* dump_arr = BufferStore::getBufferStore()->getDoubleArr(0,nSizeChank*(nDimAndAttr));

  //sending and receiving coordinates and properties of macroparticles
  for( int i = 1; i < size_MPI; i++){
    int nLoops = nSizeArr_MPI[i]/nSizeChank;
    if(nSizeArr_MPI[i] % nSizeChank != 0) nLoops++;
    for(int k = 0; k < nLoops; k++){

      int j_start = k*nSizeChank;
      int j_stop = j_start + nSizeChank;
      if(j_stop > nSizeArr_MPI[i]) j_stop = nSizeArr_MPI[i];

      if( i == rank_MPI ){
        int j_count = 0;
        for( int j = j_start; j < j_stop; j++){
          dump_arr[(nDimAndAttr)*j_count + 0] = x(j);
          dump_arr[(nDimAndAttr)*j_count + 1] = px(j);
          dump_arr[(nDimAndAttr)*j_count + 2] = y(j);
          dump_arr[(nDimAndAttr)*j_count + 3] = py(j);
          dump_arr[(nDimAndAttr)*j_count + 4] = z(j);
          dump_arr[(nDimAndAttr)*j_count + 5] = pz(j);
          dump_arr[(nDimAndAttr)*j_count + 6] = (double) flag(j);
          for(int k = 0; k < attributesSize; k++){
            dump_arr[(nDimAndAttr)*j_count + 7 + k] = getParticleAttributeVal(j,k);
          }
          j_count++;
        }
        ORBIT_MPI_Send(dump_arr, (nDimAndAttr)*nSizeChank, MPI_DOUBLE, 0,
					1111, MPI_COMM_Local);
      }
      if(rank_MPI == 0){
        ORBIT_MPI_Recv(dump_arr, (nDimAndAttr)*nSizeChank, MPI_DOUBLE, i,
					1111, MPI_COMM_Local, &statusMPI);
        for( int j = 0; j < (j_stop - j_start); j++){
          int flg = (int) dump_arr[(nDimAndAttr)*j + 6];
          if(flg > 0){
            Out<< dump_arr[(nDimAndAttr)*j + 0] << "  "
            << dump_arr[(nDimAndAttr)*j + 1] << "  "
            << dump_arr[(nDimAndAttr)*j + 2] << "  "
            << dump_arr[(nDimAndAttr)*j + 3] << "  "
            << dump_arr[(nDimAndAttr)*j + 4] << "  "
            << dump_arr[(nDimAndAttr)*j + 5] << "  ";

            for(int k = 0; k < attributesSize; k++){
              Out<<   dump_arr[(nDimAndAttr)*j + 7 + k] << "  ";
            }

            Out << std::endl;
          }

          if (j % N_flush == 0){ Out.flush();}
        }
      }
    }
  }

  if(rank_MPI == 0){Out.flush();}
  // ===== MPI end =====
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    Bunch::print
//
// DESCRIPTION
//    print every components of macroparticles with requesting a file name
//    to print for rank0 CPU. The data will be added to the existing file.
//    If you want to clear the file you have to take care about it outside
//    this class.
//
// REMAEKS
//    In order to reflect the results of the latest addParticle() ( if
//    the deleteParticle() method has been called) you have to use compress()
//    this right after all operations addParticle()
//
///////////////////////////////////////////////////////////////////////////

void Bunch::print(const char* fileName)
{
	ofstream F_dump;

	if(rank_MPI == 0){
		F_dump.open (fileName, ios::out);
	}

	print(F_dump);

	if(rank_MPI == 0){
		F_dump.close();
	}
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    Bunch::readBunch(char* fileName, int nPart)
//
// DESCRIPTION
//    reads every components of the given number of macroparticles
//    with assignning them to the concerned processes almost
//    (for the remindar) impartially.
//
// REMARK
//    needs Bunch::compress(); see Bunch::addParticle
//
// RETURNS
//    Nothing
//
///////////////////////////////////////////////////////////////////////////

int Bunch::readBunch(const char* fileName, int nParts)
{
	double x,y,z, px,py,pz;

	ifstream is;

	int error_ind = 0;

	if(rank_MPI == 0){
		is.open(fileName, std::ios::in);
		if (is.bad()){
			std::cerr << "The Bunch::readBunch(char* fileName, int nParts)"<< std::endl;
			std::cerr << "Can not open file:"<< fileName<< std::endl;
			error_ind = 1;
		}
	}

	if(size_MPI > 1){
		ORBIT_MPI_Bcast (&error_ind,1, MPI_INT,    0, MPI_COMM_Local );
	}

	if(error_ind > 0){
		ORBIT_MPI_Finalize();
	}

	//define chunk size for reading particles' coordinates
	//n_c number of particles in each cpu
	int chunk_size = 1000;
	int n_c = chunk_size/size_MPI;
	if(n_c < 1) n_c = 1;
	chunk_size = n_c * size_MPI;

	int nDimAndAttr = nDim+attributesSize;

	double* arr_0 = BufferStore::getBufferStore()->getDoubleArr(0,chunk_size*(nDimAndAttr));;

	int info_stop = 0;
	int nP = 0;
	int nT = 0;
	int indLine = 0;

	std::string  str;
	std::vector<string> v_str;

	if(rank_MPI == 0){
		getline(is,str);
		indLine++;
	}

	while(info_stop != 1){

		int nn = 0;

		if(rank_MPI == 0){

			while(!is.eof() && (nn < chunk_size) && (info_stop != 1)){
				if(strlen(str.c_str()) > 0  && str.c_str()[0] != '%'){
					//here we have the string with coordinates
					nT = StringUtils::Tokenize(str,v_str);
					if(nT != (nDimAndAttr) && nT != 6){
						std::cerr << "The Bunch::readBunch(char* fileName, int nParts)"<< std::endl;
						std::cerr << "File:"<< fileName << std::endl;
						std::cerr << "Error in data tructure in line :"<< indLine << std::endl;
						std::cerr << "line="<< str << std::endl;
						info_stop = 1;
						error_ind = 1;
						break;
					}
					for(int i = 0; i < nT; i++){
						sscanf( v_str[i].c_str(),"%lf",&arr_0[(nDimAndAttr)*nn + i]);
					}
					nn++;
				}

				if(!is.eof()) {
					getline(is,str);
					indLine++;
				}
			}

			if((nP+nn) > nParts) {
				info_stop = 1;
				nn = nParts - nP;
			}
			nP += nn;
			if(is.eof()) info_stop = 1;
		}

		if(size_MPI > 1){
			ORBIT_MPI_Bcast ( &nT,         1, MPI_INT,    0, MPI_COMM_Local );
			ORBIT_MPI_Bcast ( &error_ind,  1, MPI_INT,    0, MPI_COMM_Local );
			ORBIT_MPI_Bcast ( &nn,         1, MPI_INT,    0, MPI_COMM_Local );
			ORBIT_MPI_Bcast ( arr_0, chunk_size*(nDimAndAttr) , MPI_DOUBLE, 0, MPI_COMM_Local );
		}

		if(error_ind > 0){
			ORBIT_MPI_Finalize();
		}

		int i_start = rank_MPI;
		int i_stop  = nn;

		for(int i = i_start; i < i_stop; i += size_MPI){
			x = arr_0[0+i*nDimAndAttr];
			px = arr_0[1+i*nDimAndAttr];

			y = arr_0[2+i*nDimAndAttr];
			py = arr_0[3+i*nDimAndAttr];

			z = arr_0[4+i*nDimAndAttr];
			pz = arr_0[5+i*nDimAndAttr];

			int part_index = addParticle(x, px, y, py, z, pz);

			for(int k = 0; k < attributesSize; k++){
				if(nT != 6){
					getParticleAttributeVal(part_index,k) =  arr_0[(nDimAndAttr)*i + 6 + k];
				}
			}
		}

		if(size_MPI > 1){
			ORBIT_MPI_Bcast ( &info_stop, 1, MPI_INT, 0, MPI_COMM_Local  );
		}

	}

	if(rank_MPI == 0){
		is.close();
	}

	return getSizeGlobal();
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    Bunch::readBunch(char* fileName)
//
// DESCRIPTION
//    reads every components of the given number of macroparticles
//    with assignning them to the concerned processes almost
//    (for the remindar) impartially.
//
// RETURNS
//    Nothing
//
///////////////////////////////////////////////////////////////////////////

int Bunch::readBunch(const char* fileName)
{
	int nMax = 2000000000;
	nMax = readBunch(fileName, nMax);
	return nMax;
}

void Bunch::deleteAllParticles()
{
	for(int i = 0; i < nNew; i++){
		arrFlag[i] = 0; //dead
	}

	nSize = 0;
	nNew = 0;

	//compression is not needed
	needOfCompress = 0;
	resize();
}

void Bunch::initParticleAttributes(const char* fileName){
	removeAllParticleAttributes();
	std::vector<string> attr_names;
	readParticleAttributesNames(fileName,attr_names);
	for(int i = 0, n = attr_names.size(); i < n; i++){
		ParticleAttributes* partAttr = ParticleAttributesFactory::getParticleAttributesInstance(attr_names[i],this);
		addParticleAttributes(partAttr);
	}
}

int Bunch::readParticleAttributesNames(const char* fileName, std::vector<string>& attr_names){

	attr_names.clear();

	ifstream is;

	int error_ind = 0;

	if(rank_MPI == 0){
		is.open(fileName, std::ios::in);
		if (is.bad()){
			std::cerr << "The Bunch::readParticleAttributesNames(const char* fileName, std::vector<string>& attr_names)"<< std::endl;
			std::cerr << "Can not open file:"<< fileName<< std::endl;
			error_ind = 1;
		}
	}

	if(size_MPI > 1){
		ORBIT_MPI_Bcast (&error_ind,1, MPI_INT,    0, MPI_COMM_Local );
	}

	if(error_ind > 0){
		ORBIT_MPI_Finalize();
	}

	std::string  str;
	std::vector<string> v_str;

	if(rank_MPI == 0){
		while(!is.eof()){
			getline(is,str);
			if(strlen(str.c_str()) > 0  && str.c_str()[0] == '%'){
				int nT = StringUtils::Tokenize(str,v_str);
				if(nT > 2 && v_str[1] == "PARTICLE_ATTRIBUTES_CONTROLLERS_NAMES"){
					for(int i = 2, n = v_str.size(); i < n; i++){
						attr_names.push_back(v_str[i]);
					}
					break;
				}
			}
			else{
				break;
			}
		}
		is.close();
	}

	if(rank_MPI == 0){
	}

	int nTypes = attr_names.size();
	int strLength = strlen(str.c_str());

	if(size_MPI > 1){
		ORBIT_MPI_Bcast ( &nTypes,   1, MPI_INT,    0, MPI_COMM_Local );
		ORBIT_MPI_Bcast ( &strLength,1, MPI_INT,    0, MPI_COMM_Local );
	}
	else{
		return nTypes;
	}

	if(nTypes == 0) return 0;

	char* char_tmp = BufferStore::getBufferStore()->getCharArr(0,strLength+1);
	strcpy(char_tmp, str.c_str());
	ORBIT_MPI_Bcast ( char_tmp,  strLength+1, MPI_CHAR,    0, MPI_COMM_Local );
	std::string str_new(char_tmp);
	StringUtils::Tokenize(str_new,v_str);

	attr_names.clear();

	for(int i = 2, n = v_str.size(); i < n; i++){
		attr_names.push_back(v_str[i]);
	}

	return attr_names.size();
}

//--------------------------------------------------
//methods related to the sync. particle
//--------------------------------------------------
SyncPart* Bunch::getSyncPart(){
	return syncPart;
}

//--------------------------------------------------
//methods related to the particle attribute buckets
//--------------------------------------------------
double& Bunch::getParticleAttributeVal(int ind, int attr_ind){
	return arrAttr[ind][attr_ind];
}

void Bunch::addParticleAttributes(const std::string att_name){
	if(attrCntrSizeMap.count(att_name) > 0) return;
	ParticleAttributes* attr = ParticleAttributesFactory::getParticleAttributesInstance(att_name,this);
	addParticleAttributes(attr);
}

int Bunch::hasParticleAttributes(const std::string att_name){
	return attrCntrSizeMap.count(att_name);
}

void Bunch::addParticleAttributes(ParticleAttributes* attr){

	if(attr == NULL) {
		if(rank_MPI == 0){
			std::cerr << "void Bunch::addParticleAttributes(ParticleAttributes* attr)"<< std::endl;
			std::cerr << "You are trying to add to the bunch EMPTY particle attributes"<< std::endl;
		}
		ORBIT_MPI_Finalize("Bunch::addParticleAttributes. Stop.");
		return;
	}

	if(attrCntrSizeMap.count(attr->name()) > 0 || attr->getAttSize() == 0) {
		if(rank_MPI == 0){
			std::cerr << "void Bunch::addParticleAttributes(ParticleAttributes* attr)"<< std::endl;
			std::cerr << "These Attributes already inside the bunch. Name ="<< attr->name() << std::endl;
			std::cerr << "Or the attributes size equals to = "<< attr->getAttSize() << std::endl;
		}
		ORBIT_MPI_Finalize("Bunch::addParticleAttributes. Stop.");
		return;
	}

	if(this != attr->bunch()) {
		if(rank_MPI == 0){
			std::cerr << "void Bunch::addParticleAttributes(ParticleAttributes* attr)"<< std::endl;
			std::cerr << "The Attributes Name = " << attr->name() << std::endl;
			std::cerr << "For this attributes the bunch is different from this!" << std::endl;
		}
		ORBIT_MPI_Finalize("Bunch::addParticleAttributes. Stop.");
		return;
	}

	attrCntrMap[attr->name()] = attr;
	int attr_length = attr->getAttSize();
	if(attrCntrSize == 0){
		arrAttr =  new double*[nTotalSize];
		for(int i=0; i < nTotalSize; i++){
			arrAttr[i] = new double[attr_length];
		}
	}
	else{
		double* tmp_arr = NULL;
		for(int i=0; i < nTotalSize; i++){
			tmp_arr = arrAttr[i];
			arrAttr[i] = new double[attr_length+attributesSize];
			for(int j = 0; j < attributesSize; j++){
				arrAttr[i][j] = tmp_arr[j];
			}
			delete [] tmp_arr;
		}
	}

	attr->setAttrShift(attributesSize);
	attrCntrSizeMap[attr->name()] = attr->getAttSize();
	attrCntrLowIndMap[attr->name()] = attributesSize;
	attributesSize += attr_length;
	attrCntrUppIndMap[attr->name()] = attributesSize;
	attrCntrSize++;
	for(int i=0; i < nTotalSize; i++){
		attr->init(i);
	}

	//memorize attributes in the set.
	//This set will be used in the destructor.
	attrCntrSet.insert(attr);
}

void Bunch::removeAllParticleAttributes(){
	std::vector<string> names;
	getParticleAttributesNames(names);
	for(int i = 0, n= names.size(); i < n; i++){
		removeParticleAttributes(names[i]);
	}
}

void Bunch::removeParticleAttributes(const string name){
	ParticleAttributes* attr = removeParticleAttributesWithoutDelete(name);

	//delete Particle attribute itself
	if(attr != NULL){
		delete attr;
	}
}

ParticleAttributes* Bunch::removeParticleAttributesWithoutDelete(const string name){
	if(attrCntrSizeMap.count(name) == 0) return NULL;

	ParticleAttributes* attr = attrCntrMap[name];
	int attr_length = attr->getAttSize();
	int lowInd = attrCntrLowIndMap[name];
	int uppInd = attrCntrUppIndMap[name];

	int newAttributesSize = attributesSize - attr_length;

	if( newAttributesSize < 0){
		if(rank_MPI == 0){
			std::cerr << "======The structure of particle attributes is wrong ====="<< std::endl;
			std::cerr << "Bunch::removeParticleAttributesWithoutDelete(const string name)"<< std::endl;
			std::cerr << "attributesSize ="<<attributesSize << std::endl;
			std::cerr << "attr_length = "<< attr_length << std::endl;
			std::cerr << "name:"<< name<< std::endl;
		}
		ORBIT_MPI_Finalize();
	}

	if(newAttributesSize > 0){
		double* tmp_arr = NULL;
		int aInd = 0;
		for(int i=0; i < nTotalSize; i++){
			tmp_arr = arrAttr[i];
			arrAttr[i] = new double[newAttributesSize];
			aInd = 0;
			for(int j = 0; j < lowInd; j++){
				arrAttr[i][aInd] = tmp_arr[j];
				aInd++;
			}
			for(int j = uppInd; j < attributesSize; j++){
				arrAttr[i][aInd] = tmp_arr[j];
				aInd++;
			}
			delete [] tmp_arr;
		}
	}
	else{
		for(int i=0; i < nTotalSize; i++){
			delete [] arrAttr[i];
		}
		delete [] arrAttr;
	}

	attributesSize = newAttributesSize;
	attrCntrSize--;

	attrCntrMap.erase(name);
	attrCntrSizeMap.erase(name);
	attrCntrLowIndMap.erase(name);
	attrCntrUppIndMap.erase(name);

	//remove attributes in the set.
	attrCntrSet.erase(attr);

	std::map<std::string,int>::iterator pos;
	for (pos = attrCntrLowIndMap.begin(); pos != attrCntrLowIndMap.end(); ++pos) {
		string name = pos->first;
		int ind = pos->second;
		if(ind >= lowInd) {
			attrCntrLowIndMap[name] =  attrCntrLowIndMap[name] - attr_length;
			attrCntrUppIndMap[name] =  attrCntrUppIndMap[name] - attr_length;
		}
	}

	return attr;
}


ParticleAttributes* Bunch::getParticleAttributes(const string name){
	if(attrCntrSizeMap.count(name) == 0) {
		if(rank_MPI == 0){
			std::cerr << "ParticleAttributes* Bunch::getParticleAttributes(const string name)"<< std::endl;
			std::cerr << "There is no ParticleAttributes with this name."<< std::endl;
			std::cerr << "name:"<< name<< std::endl;
		}
		ORBIT_MPI_Finalize();
		return NULL;
	}
	return attrCntrMap[name];
}

void Bunch::getParticleAttributesNames(std::vector<string>& names){
	names.clear();
	std::map<std::string,ParticleAttributes*>::iterator pos;
	for (pos = attrCntrMap.begin(); pos != attrCntrMap.end(); ++pos) {
		string name = pos->first;
		names.push_back(name);
	}
}

void Bunch::clearAllParticleAttributesAndMemorize(){
	std::map<std::string,ParticleAttributes*>::iterator pos;
	for (pos = attrCntrMap.begin(); pos != attrCntrMap.end(); ++pos) {
		std::string name = pos->first;
		if(attrCntrMapTemp.count(name) == 0){
			ParticleAttributes* attrCntrl = pos->second;
			attrCntrMapTemp[name] = attrCntrl;
		}
		else{
			removeParticleAttributes(name);
		}
	}

	for (pos = attrCntrMap.begin(); pos != attrCntrMap.end(); ++pos) {
		string name = pos->first;
		removeParticleAttributesWithoutDelete(name);
	}
}

void Bunch::restoreAllParticleAttributesFromMemory(){

	std::map<std::string,ParticleAttributes*>::iterator pos;
	for (pos = attrCntrMap.begin(); pos != attrCntrMap.end(); ++pos) {
		string name = pos->first;
		removeParticleAttributes(name);
	}

	for (pos = attrCntrMapTemp.begin(); pos != attrCntrMapTemp.end(); ++pos) {
		std::string name = pos->first;
		ParticleAttributes* attrCntrl = pos->second;
		addParticleAttributes(attrCntrl);
	}
	attrCntrMapTemp.clear();
}

MPI_Comm  Bunch::getMPI_Comm_Local(){
	return MPI_COMM_Local;
}

void  Bunch::setMPI_Comm_Local(MPI_Comm MPI_COMM_Local){
	this->MPI_COMM_Local = MPI_COMM_Local;
}

int Bunch::getMPI_Size(){
	return size_MPI;
}

int Bunch::getMPI_Rank(){
	return rank_MPI;
}

