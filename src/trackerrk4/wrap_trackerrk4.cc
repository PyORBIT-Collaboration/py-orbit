#include "orbit_mpi.hh"

#include "wrap_trackerrk4.hh"
#include "wrap_runge_kutta_tracker.hh"
#include "wrap_py_external_effects.hh"

static PyMethodDef trackerrk4Methods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void inittrackerrk4(){
    //create new module
    PyObject* module = Py_InitModule("trackerrk4",trackerrk4Methods);
		//add the other classes init
		wrap_trackerrk4::initRungeKuttaTracker(module);
		wrap_trackerrk4_py_external_effects::initPyExternalEffects(module);
  }
	
#ifdef __cplusplus
}
#endif

