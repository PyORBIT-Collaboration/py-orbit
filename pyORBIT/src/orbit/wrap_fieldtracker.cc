//////////////////////////////// -*- C++ -*- //////////////////////////////
//
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "wrap_fieldtracker.hh"
#include "wrap_multexpansion3D.hh"
#include "wrap_magfieldtracker3D.hh"

namespace wrap_orbit_fieldtracker{

static PyMethodDef FieldTrackerMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initfieldtracker(){
    //create new module
    PyObject* module = Py_InitModule("fieldtracker",FieldTrackerMethods);

		//add the MultipoleExpansion python class
		wrap_orbit_multexpansion3D::initmultexpansion3D(module);
      
      //add the MagneticFieldTracker3D class
      wrap_orbit_magfieldtracker3D::initmagfieldtracker3D(module);
  }


#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_fieldtracker
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
