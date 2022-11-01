#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_aperture.hh"

#include "wrap_TAperture.hh"
#include "wrap_PhaseAperture.hh"
#include "wrap_EnergyAperture.hh"
#include "wrap_BaseAperture.hh"
#include "wrap_PyBaseApertureShape.hh"
#include "wrap_PrimitiveApertureShape.hh"
#include "wrap_CompositeApertureShape.hh"
#include "wrap_ConvexApertureShape.hh"

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
		wrap_base_aperture::initBaseAperture(module);
		wrap_py_base_aperture_shape::initPyBaseApertureShape(module);
		wrap_primitive_aperture_shape::initPrimitiveApertureShape(module);
		wrap_py_composite_aperture_shape::initCompositeApertureShape(module);
		wrap_convex_aperture_shape::initConvexApertureShape(module);
	}

#ifdef __cplusplus
}
#endif


}
