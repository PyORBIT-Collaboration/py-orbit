#include "orbit_mpi.hh"

#include "wrap_grid1D.hh"
#include "wrap_grid2D.hh"
#include "wrap_grid3D.hh"
#include "wrap_poissonsolverfft2d.hh"
#include "wrap_poissonsolverfft3d.hh"
#include "wrap_boundary2d.hh"
#include "wrap_spacecharge.hh"
#include "wrap_spacechargecalc2p5d.hh"
#include "wrap_spacechargecalc2p5d_rb.hh"
#include "wrap_spacechargecalc3d.hh"
#include "wrap_lspacechargecalc.hh"
#include "wrap_uniform_ellipsoid_field_calculator.hh"
#include "wrap_spacechargecalc_uniform_ellipse.hh"

static PyMethodDef spacechargeMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initspacecharge(){
    //create new module
    PyObject* module = Py_InitModule("spacecharge",spacechargeMethods);
		//add the other classes init
		wrap_spacecharge::initGrid1D(module);
		wrap_spacecharge::initGrid2D(module);
		wrap_spacecharge::initGrid3D(module);		
		wrap_spacecharge::initPoissonSolverFFT2D(module);
		wrap_spacecharge::initPoissonSolverFFT3D(module);
		wrap_spacecharge::initBoundary2D(module);
		wrap_spacecharge::initSpaceChargeCalc2p5D(module);
		wrap_spacecharge::initSpaceChargeCalc2p5Drb(module);
		wrap_spacecharge::initSpaceChargeCalc3D(module);
		wrap_spacecharge::initUniformEllipsoidFieldCalculator(	module);
		wrap_spacecharge::initSpaceChargeCalcUniformEllipse(	module);
		wrap_lspacechargecalc::initLSpaceChargeCalc(module);
  }
	
	PyObject* getSpaceChargeType(char* name){
		PyObject* mod = PyImport_ImportModule("spacecharge");
		PyObject* pyType = PyObject_GetAttrString(mod,name);
		Py_DECREF(mod);
		Py_DECREF(pyType);
		return pyType;
	}
		
#ifdef __cplusplus
}
#endif
