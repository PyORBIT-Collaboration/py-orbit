#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_aperture.hh"

#include "wrap_TAperture.hh"
#include "wrap_PhaseAperture.hh"
#include "wrap_EnergyAperture.hh"

#include "wrap_bunch.hh"

#include <iostream>

namespace wrap_aperture{

#ifdef __cplusplus
extern "C" {
#endif
	
  static PyMethodDef ApertureModuleMethods[] = { {NULL,NULL} };	
  
	//--------------------------------------------------
	//Initialization aperture module
	//--------------------------------------------------

	void initaperture(){
		//create new module
		PyObject* module = Py_InitModule("aperture",ApertureModuleMethods);
		wrap_aperture::initTAperture(module);
		wrap_phase_aperture::initPhaseAperture(module);
		wrap_energy_aperture::initEnergyAperture(module);
	}

#ifdef __cplusplus
}
#endif


}
