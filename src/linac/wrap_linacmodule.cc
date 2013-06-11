#include "orbit_mpi.hh"

#include "wrap_linacmodule.hh"
#include "wrap_BaseRfGap.hh"
#include "wrap_MatrixRfGap.hh"
#include "wrap_RfGapTTF.hh"
#include "wrap_SuperFishFieldSource.hh"
#include "wrap_RfGapThreePointTTF.hh"

static PyMethodDef linacmoduleMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif
  namespace wrap_linac{
   void initlinac(){
		 //create new module
		 PyObject* module = Py_InitModule("linac",linacmoduleMethods);
		 //add the other classes init
		 wrap_linac::initBaseRfGap(module);
		 wrap_linac::initMatrixRfGap(module);
		 wrap_linac::initRfGapTTF(module);	
		 wrap_linac::initSuperFishFieldSource(module);
		 wrap_linac::initRfGapThreePointTTF(module);
	 }
	 
	 PyObject* getLinacType(char* name){
		 PyObject* mod = PyImport_ImportModule("linac");
		 PyObject* pyType = PyObject_GetAttrString(mod,name);
		 Py_DECREF(mod);
		 Py_DECREF(pyType);
		 return pyType;
	 }
	}
		
#ifdef __cplusplus
}
#endif
