//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    SyncPart.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    12/14/2005
//
// DESCRIPTION
//    Source code for the synchronous particle class. It keeps info
//    about energy, momentum etc. of the synchronous macro-particle
//    The synchronous particle bears the coordinate system with z-axis 
//    directing along pz vector, and x-axis perpendicular to z. The y-axis
//    ort direction is defined as ny = [nz*nx].
//
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "Python.h"

#include "SyncPart.hh"

#include "Bunch.hh"
#include "OrbitConst.hh"
#include "StringUtils.hh"
#include "BufferStore.hh"

#include <iomanip>
#include <string>

//the limit p/m when we will use non/relativistic formulas
#define P_OVER_M_LIMIT 1.0e-5

using namespace OrbitUtils;

///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

SyncPart::SyncPart(Bunch* bunchIn)
{

	 bunch = bunchIn;

  //initialization all necessary variables and attributes
  //such as energy, time etc.
  init();

}

SyncPart::~SyncPart()
{
	Py_XDECREF(py_wrapper);
}

void SyncPart::setPyWrapper(PyObject* py_wrapper_In){
	if(py_wrapper != NULL) Py_XDECREF(py_wrapper);
	Py_INCREF(py_wrapper_In);
	py_wrapper = py_wrapper_In;
}

PyObject* SyncPart::getPyWrapper(){
	return py_wrapper;
}

//initialization all necessary variables and attributes
void SyncPart::init()
{
	py_wrapper = NULL;

	energy = 0.;
	time = 0.;

  p_abs = 0.;
	beta = 0.;
	gamma = 1.0;

	xyz[0] = 0.;
	xyz[1] = 0.;
	xyz[2] = 0.;

	pxyz[0] = 0.;
	pxyz[1] = 0.;
	pxyz[2] = 0.;
	
	x_normal[0] = 1.;
	x_normal[1] = 0.;
	x_normal[2] = 0.;
}

// Kinetic energy in GeV
double SyncPart::getEnergy(){
	return energy;
}

void SyncPart::setTime(double time){
	this->time = time;
}

double SyncPart::getTime(){
	return time;
}

void SyncPart::setXYZ(const double* xyz){
	this->xyz[0] = xyz[0];
	this->xyz[1] = xyz[1];
	this->xyz[2] = xyz[2];
}

void SyncPart::setXYZ(double x, double y, double z){
	this->xyz[0] = x;
	this->xyz[1] = y;
	this->xyz[2] = z;
}

void SyncPart::setX(double x){
	this->xyz[0] = x;
}

void SyncPart::setY(double y){
	this->xyz[1] = y;
}

void SyncPart::setZ(double z){
	this->xyz[2] = z;
}

void SyncPart::setMomentum(double p){
	if(p_abs == 0.){
		this->pxyz[0] = 0.;
		this->pxyz[1] = 0.;
		this->pxyz[2] = p;
	} else {
		double coef = p/p_abs;
		this->pxyz[0] = coef*this->pxyz[0];
		this->pxyz[1] = coef*this->pxyz[1];
		this->pxyz[2] = coef*this->pxyz[2];
	}
	updateKinematics();
}

double SyncPart::getX(){
	return xyz[0];
}

double SyncPart::getY(){
	return xyz[1];
}

double SyncPart::getZ(){
	return xyz[2];
}

void SyncPart::setPXYZ(const double* pxyz){
	this->pxyz[0] = pxyz[0];
	this->pxyz[1] = pxyz[1];
	this->pxyz[2] = pxyz[2];
	updateKinematics();
}

void SyncPart::setPXYZ(double px, double py, double pz){
	this->pxyz[0] = px;
	this->pxyz[1] = py;
	this->pxyz[2] = pz;
	updateKinematics();
}

void SyncPart::setPX(double px){
	this->pxyz[0] = px;
	updateKinematics();
}

void SyncPart::setPY(double py){
	this->pxyz[1] = py;
	updateKinematics();
}

void SyncPart::setPZ(double pz){
	this->pxyz[2] = pz;
	updateKinematics();
}

double SyncPart::getPX(){
	return pxyz[0];
}

double SyncPart::getPY(){
	return pxyz[1];
}

double SyncPart::getPZ(){
	return pxyz[2];
}

double SyncPart::getMomentum(){
	return p_abs;
}

double SyncPart::getBeta(){
	return beta;
}

double SyncPart::getGamma(){
	return gamma;
}

void SyncPart::setNormalX(double x, double y, double z){
	x_normal[0] = x;
	x_normal[1] = y;
	x_normal[2] = z;	
	updateKinematics();
}

double SyncPart::getNormalXX(){
	return x_normal[0];
}
	
double SyncPart::getNormalXY(){
	return x_normal[1];
}
	
double SyncPart::getNormalXZ(){
	return x_normal[2];
}

double SyncPart::getNormalYX(){
	return y_normal[0];
}
	
double SyncPart::getNormalYY(){
	return y_normal[1];
}
	
double SyncPart::getNormalYZ(){
	return y_normal[2];
}

void SyncPart::updateKinematics(){
	double p2 = pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[2]*pxyz[2];
	p_abs = sqrt(p2);
	//create x-axis ort
	if(p_abs != 0.){
		int iter = 0;
		double x_normal_abs = 0.;
		do {
			if(iter >= 1){
				x_normal[0] = x_normal[0] + iter*1.0;
				if((iter % 2) == 0) x_normal[1] = x_normal[1] + iter*1.0;
				if((iter % 3) == 0) x_normal[2] = x_normal[2] + iter*1.0;
			}
			double coef = x_normal[0]*pxyz[0] + x_normal[1]*pxyz[1] + x_normal[2]*pxyz[2];
			coef = coef/p2;
			x_normal[0] = x_normal[0] - coef*pxyz[0];
			x_normal[1] = x_normal[1] - coef*pxyz[1];
			x_normal[2] = x_normal[2] - coef*pxyz[2];
			x_normal_abs = sqrt(x_normal[0]*x_normal[0] + x_normal[1]*x_normal[1] + x_normal[2]*x_normal[2]);
			iter = iter + 1;
		} while(x_normal_abs == 0.);
		x_normal[0] = x_normal[0]/x_normal_abs;
		x_normal[1] = x_normal[1]/x_normal_abs;
		x_normal[2] = x_normal[2]/x_normal_abs;
	} else {
		double x_normal_abs = sqrt(x_normal[0]*x_normal[0] + x_normal[1]*x_normal[1] + x_normal[2]*x_normal[2]);
		if(x_normal_abs == 0.){
			x_normal[0] = 1.;
			x_normal[1] = 0.;
			x_normal[2] = 0.;
		} else {
			x_normal[0] = x_normal[0]/x_normal_abs;
			x_normal[1] = x_normal[1]/x_normal_abs;
			x_normal[2] = x_normal[2]/x_normal_abs;
		}
	}
	//create y-axis ort
	double y_normal_abs = 0.;
	if(p_abs != 0.){
		y_normal[0] = pxyz[1]*x_normal[2] - pxyz[2]*x_normal[1];
		y_normal[1] = -(pxyz[0]*x_normal[2] - pxyz[2]*x_normal[0]);
		y_normal[0] = pxyz[0]*x_normal[1] - pxyz[1]*x_normal[0];
		y_normal_abs = sqrt(y_normal[0]*y_normal[0] + y_normal[1]*y_normal[1] + y_normal[2]*y_normal[2]);
	} else {
		int iter = 0;
		do {
			if(iter >= 1){
				y_normal[0] = y_normal[0] + iter*1.0;
				if((iter % 2) == 0) y_normal[1] = y_normal[1] + iter*1.0;
				if((iter % 3) == 0) y_normal[2] = y_normal[2] + iter*1.0;
			}
			double coef = y_normal[0]*x_normal[0] + y_normal[1]*x_normal[1] + y_normal[2]*x_normal[2];
			y_normal[0] = y_normal[0] - coef*x_normal[0];
			y_normal[1] = y_normal[1] - coef*x_normal[1];
			y_normal[2] = y_normal[2] - coef*x_normal[2];
			y_normal_abs = sqrt(y_normal[0]*y_normal[0] + y_normal[1]*y_normal[1] + y_normal[2]*y_normal[2]);
			iter = iter + 1;
		} while(y_normal_abs == 0.);
	}
	y_normal[0] = y_normal[0]/y_normal_abs;
	y_normal[1] = y_normal[1]/y_normal_abs;
	y_normal[2] = y_normal[2]/y_normal_abs;
	//--- particle kinematics
	double m = bunch->getMass();
	double m2 = m*m;
	double w2 = m2+p2;
	double w = sqrt(w2);
	//relativistic or non-relativistic approach
	if(p_abs/m < P_OVER_M_LIMIT){
		energy = p2/(2.0*m);
	}
	else{
		energy = w - m;
	}
	beta = p_abs/w;
	gamma = w/m;
}

double SyncPart::momentumToEnergy(double p){
	double m = bunch->getMass();
	double ek = 0.;
	//relativistic or non-relativistic approach
	if(p/m < P_OVER_M_LIMIT){
		ek = p*p/(2.0*m);
	}
	else{
		ek = sqrt(m*m+p*p)-m;
	}
	return ek;
}

double SyncPart::energyToMomentum(double ek){
	double m = bunch->getMass();
	return sqrt(ek*(ek+2.0*m));
}

double SyncPart::getMass(){
	return bunch->getMass();
}

void SyncPart::readSyncPart(const char* fileName){

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

  std::vector<std::string> attr_names;
  attr_names.clear();

  ifstream is;

  int error_ind = 0;
  if(rank_MPI == 0){
    is.open(fileName, std::ios::in);
    if (is.bad()){
      std::cerr << "The SyncPart::initSyncPart(char* fileName)"<< std::endl;
      std::cerr << "Can not open file:"<< fileName<< std::endl;
      error_ind = 1;
    }
  }

  if(size_MPI > 1){
    ORBIT_MPI_Bcast (&error_ind,1, MPI_INT,    0, MPI_COMM_Local );
  }

  if(error_ind > 0){
    ORBIT_MPI_Finalize("The SyncPart::initSyncPart. Stop.");
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
          if(nT > 5 &&
						( v_str[1] == "SYNC_PART_COORDS" ||
							v_str[1] == "SYNC_PART_MOMENTUM" || 
							v_str[1] == "SYNC_PART_X_AXIS" )
						){
            def_found_ind = 1;
          }
          if(nT > 3 && v_str[1] == "SYNC_PART_TIME"){
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
				int buff_index = 0;
        char* char_tmp = BufferStore::getBufferStore()->getFreeCharArr(buff_index,strLength+1);
        strcpy(char_tmp, str.c_str());
        ORBIT_MPI_Bcast ( char_tmp,  strLength+1, MPI_CHAR,    0, MPI_COMM_Local );
        std::string str_new(char_tmp);
        StringUtils::Tokenize(str_new,v_str);
				BufferStore::getBufferStore()->setUnusedCharArr(buff_index);
      }

			//set coordinates
      if(v_str.size() > 5 && v_str[1] == "SYNC_PART_COORDS"){
        double val_arr[3];
				val_arr[0] = 0.;
				val_arr[1] = 0.;
				val_arr[2] = 0.;
        sscanf( v_str[3].c_str(),"%lf",&val_arr[0]);
				sscanf( v_str[4].c_str(),"%lf",&val_arr[1]);
				sscanf( v_str[5].c_str(),"%lf",&val_arr[2]);
				setXYZ(val_arr);
      }

			//set momentum
      if(v_str.size() > 5 && v_str[1] == "SYNC_PART_MOMENTUM"){
        double val_arr[3];
				val_arr[0] = 0.;
				val_arr[1] = 0.;
				val_arr[2] = 0.;
        sscanf( v_str[3].c_str(),"%lf",&val_arr[0]);
				sscanf( v_str[4].c_str(),"%lf",&val_arr[1]);
				sscanf( v_str[5].c_str(),"%lf",&val_arr[2]);
				setPXYZ(val_arr);
      }

			//set x-axis ort
      if(v_str.size() > 5 && v_str[1] == "SYNC_PART_X_AXIS"){
        double val_arr[3];
				val_arr[0] = 0.;
				val_arr[1] = 0.;
				val_arr[2] = 0.;
        sscanf( v_str[3].c_str(),"%lf",&val_arr[0]);
				sscanf( v_str[4].c_str(),"%lf",&val_arr[1]);
				sscanf( v_str[5].c_str(),"%lf",&val_arr[2]);
				setNormalX(val_arr[0],val_arr[1],val_arr[2]);
      }
			
			//set time
      if(v_str.size() > 3 && v_str[1] == "SYNC_PART_TIME"){
        double val = 0.;
        sscanf( v_str[3].c_str(),"%lf",&val);
        setTime(val);
      }
    }

  }

  if(rank_MPI == 0){
    is.close();
  }
}

void SyncPart::print(std::ostream& Out)
{

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

  //single CPU case
  if(rank_MPI == 0){

		Out <<std::setprecision(10);

    //print coordinates
    Out << "%  SYNC_PART_COORDS ";
		Out << getX()  <<" ";
		Out << getY()  <<" ";
		Out << getZ()  <<" ";
		Out <<" x, y, z positions in [m]";
    Out << std::endl;

    //print momentum
    Out << "%  SYNC_PART_MOMENTUM ";
		Out << getPX()  <<" ";
		Out << getPY()  <<" ";
		Out << getPZ()  <<" ";
		Out <<" px, py, pz momentum component in GeV/c";
    Out << std::endl;

    //print x-axis ort
    Out << "%  SYNC_PART_X_AXIS ";
		Out << getNormalXX()  <<" ";
		Out << getNormalXY()  <<" ";
		Out << getNormalXZ()  <<" ";
		Out <<" nxx, nxy, pxz - x-axis ort coordinates";
    Out << std::endl;

    //print energy
    Out << "%  info only: energy of the synchronous particle [GeV] = ";
		Out << getEnergy()  <<" ";
    Out << std::endl;

    //print momentum
    Out << "%  info only: momentum of the synchronous particle [GeV/c] = ";
		Out << getMomentum()  <<" ";
    Out << std::endl;

    //print beta
    Out << "%  info only: beta=v/c of the synchronous particle = ";
		Out << getBeta()  <<" ";
    Out << std::endl;

    //print gamma
    Out << "%  info only: gamma=1/sqrt(1-(v/c)**2) of the synchronous particle = ";
		Out << getGamma()  <<" ";
    Out << std::endl;

    //print time
    Out << "%  SYNC_PART_TIME ";
		Out << getTime()  <<" ";
		Out <<" time in [sec]";
    Out << std::endl;
  }

  if(rank_MPI == 0){Out.flush();}
  // ===== MPI end =====
}

