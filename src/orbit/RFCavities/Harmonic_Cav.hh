#ifndef HARMONIC_CAV_H
#define HARMONIC_CAV_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

using namespace std;

class Harmonic_Cav: public OrbitUtils::CppPyWrapper
{
  public:
    Harmonic_Cav(double ZtoPhi, double dESync   ,
                 double RFHNum, double RFVoltage,
                 double RFPhase);
    virtual ~Harmonic_Cav();
    void   setZtoPhi(double ZtoPhi);
    double getZtoPhi();
    void   setdESync(double dESync);
    double getdESync();
    void   setRFHNum(double RFHNum);
    double getRFHNum();
    void   setRFVoltage(double RFVoltage);
    double getRFVoltage();
    void   setRFPhase(double RFPhase);
    double getRFPhase();
    void   trackBunch(Bunch* bunch);

  private:
  double _ZtoPhi;
  double _dESync;
  double _RFHNum;
  double _RFVoltage;
  double _RFPhase;

  protected:

};

#endif
