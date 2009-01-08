#include "orbit_mpi.hh"

#include "wrap_grid2D.hh"
#include "wrap_poissonsolverfft2d.hh"

static PyMethodDef spacechargeMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initspacecharge(){
    //create new module
    PyObject* module = Py_InitModule("spacecharge",spacechargeMethods);
		//add the other classes init
		wrap_spacecharge::initGrid2D(module);
		wrap_spacecharge::initPoissonSolverFFT2D(module);
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
