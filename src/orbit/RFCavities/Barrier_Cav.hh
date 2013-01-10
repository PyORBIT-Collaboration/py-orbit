#ifndef BARRIER_CAV_H
#define BARRIER_CAV_H

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

class Barrier_Cav: public OrbitUtils::CppPyWrapper
{
  public:
    Barrier_Cav(double ZtoPhi,    double RFVoltage,
                double RFPhasep,  double RFPhasem,
                double dRFPhasep, double dRFPhasem);
    virtual ~Barrier_Cav();
    void   setZtoPhi(double ZtoPhi);
    double getZtoPhi();
    void   setRFVoltage(double RFVoltage);
    double getRFVoltage();
    void   setRFPhasep(double RFPhasep);
    double getRFPhasep();
    void   setRFPhasem(double RFPhasem);
    double getRFPhasem();
    void   setdRFPhasep(double dRFPhasep);
    double getdRFPhasep();
    void   setdRFPhasem(double dRFPhasem);
    double getdRFPhasem();
    void   trackBunch(Bunch* bunch);

  private:
  double _ZtoPhi;
  double _RFVoltage;
  double _RFPhasep;
  double _RFPhasem;
  double _dRFPhasep;
  double _dRFPhasem;

  protected:

};

#endif
