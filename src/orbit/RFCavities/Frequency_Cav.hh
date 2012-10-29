#ifndef FREQUENCY_CAV_H
#define FREQUENCY_CAV_H

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

class Frequency_Cav: public OrbitUtils::CppPyWrapper
{
  public:
    Frequency_Cav(double RFFreq, double RFE0TL, double RFPhase);
    virtual ~Frequency_Cav();
    void   setRFFreq(double RFFreq);
    double getRFFreq();
    void   setRFE0TL(double RFE0TL);
    double getRFE0TL();
    void   setRFPhase(double RFPhase);
    double getRFPhase();
    void   trackBunch(Bunch* bunch);

  private:
    double _RFFreq;
    double _RFE0TL;
    double _RFPhase;

  protected:

};

#endif
