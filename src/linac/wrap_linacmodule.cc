#include "orbit_mpi.hh"

#include "wrap_linacmodule.hh"
#include "wrap_BaseRfGap.hh"

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
