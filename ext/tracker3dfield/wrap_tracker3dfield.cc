#include "orbit_mpi.hh"

#include "wrap_tracker3dfield.hh"
#include "wrap_runge_kutta_tracker.hh"
#include "wrap_py_external_effects.hh"
#include "wrap_las_strip_external_effects.hh"
#include "wrap_cpp_base_field_source.hh"

static PyMethodDef tracker3dfieldMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void inittracker3dfield(){
    //create new module
    PyObject* module = Py_InitModule("tracker3dfield",tracker3dfieldMethods);
		//add the other classes init
		wrap_tracker3dfield::initRungeKuttaTracker(module);
		wrap_tracker3dfield_py_external_effects::initPyExternalEffects(module);
		wrap_tracker3dfield_las_strip_external_effects::initLasStripExternalEffects(module);
		wrap_tracker3dfield_cpp_base_field_source::initCppBaseFieldSource(module);

  }
	
#ifdef __cplusplus
}
#endif

