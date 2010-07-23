//

#ifndef SC_SPACECHARGE_CALC_2P5D_H
#define SC_SPACECHARGE_CALC_2P5D_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Grid1D.hh"
#include "Grid2D.hh"
//#include "Boundary2D.hh"

using namespace std;

class SpaceChargeCalc2p5D: public OrbitUtils::CppPyWrapper
{
public:
	
	/** Constructor */
	SpaceChargeCalc2p5D(double lkick);

	/** Destructor */
	virtual ~SpaceChargeCalc2p5D();
	
	void calcPotential(Bunch* bunch, Grid2D* phiGrid, Grid1D* zGrid);

	
private:
	//memory allocation and step calculation for dx_ and dy_ 
	void init();
	
protected:
	double lkick;
};
//end of SC_SPACECHARGE_CALC_2P5D_H
#endif
