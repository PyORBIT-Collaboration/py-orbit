#include "Python.h"
#include "orbit_mpi.hh"

#include "pyORBIT_Object.hh"

#include "linac_tracking.hh"

#include "wrap_linac_tracking.hh"

namespace wrap_linac_tracking
{
    void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C"
{
#endif

  //---------------------------------------------------------
  //linac_tracking method wrappers
  //---------------------------------------------------------

	//Tracking a bunch through a linac drift element wrapper
	static PyObject* wrap_linac_drift(PyObject *self, PyObject *args)
	{
		PyObject* pyBunch;
		double length;
		if(!PyArg_ParseTuple(	args, "Od:linac_drift",
			&pyBunch, &length))
		{
			error("linac tracking - linac_drift - cannot parse arguments!");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
		linac_tracking::linac_drift(cpp_bunch, length);
		Py_INCREF(Py_None);
		return Py_None;
	}
	
	//Tracking a bunch through a linac quad1 element (linear part of quad tracking) wrapper
	static PyObject* wrap_linac_quad1(PyObject *self, PyObject *args)
	{
		PyObject* pyBunch;
		double length;
		double kq;
		int useCharge = 1;
		if(!PyArg_ParseTuple(	args, "Odd:linac_quad1",
			&pyBunch, &length,&kq))
		{
			error("linac tracking - linac_quad1(pyBunch,length,kq) - cannot parse arguments!");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
		linac_tracking::linac_quad1(cpp_bunch, length, kq, useCharge);
		Py_INCREF(Py_None);
		return Py_None;
	}   		

	//Tracking a bunch through a linac quad2 element (non linear part of quad tracking) wrapper
	static PyObject* wrap_linac_quad2(PyObject *self, PyObject *args)
	{
		PyObject* pyBunch;
		double length;
		if(!PyArg_ParseTuple(	args, "Od:linac_quad2",
			&pyBunch, &length))
		{
			error("linac tracking - linac_quad2 - cannot parse arguments!");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
		linac_tracking::linac_quad2(cpp_bunch, length);
		Py_INCREF(Py_None);
		return Py_None;
	}   
	
	//Tracking a bunch through a linac quad3 element (non-linear part of quad tracking) wrapper
	static PyObject* wrap_linac_quad3(PyObject *self, PyObject *args)
	{
		PyObject* pyBunch;
		double length;
		double dB_dz = 0.;
		int useCharge = 1;
		if(!PyArg_ParseTuple(	args, "Odd:linac_quad3",
			&pyBunch, &length,&dB_dz))
		{
			error("linac tracking - linac_quad3(pyBunch,length,dB_dz) - cannot parse arguments!");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
		linac_tracking::linac_quad3(cpp_bunch, length, dB_dz);
		Py_INCREF(Py_None);
		return Py_None;
	} 	
	
	//Tracking a bunch through a linac kicker element with different kick for different energies
	static PyObject* wrap_linac_kick(PyObject *self, PyObject *args)
	{
		PyObject* pyBunch;
		double kx, ky, kE;
		int useCharge = 1;
		if(!PyArg_ParseTuple(	args, "Oddd|i:kick",
			&pyBunch, &kx, &ky, &kE, &useCharge))
		{
			error("linac tracking - linac_kick - cannot parse arguments!");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
		linac_tracking::kick(cpp_bunch, kx, ky, kE, useCharge);
		Py_INCREF(Py_None);
		return Py_None;
	} 	
	
	static PyMethodDef linactrackingMethods[] =
	{
		{"drift",            wrap_linac_drift,          METH_VARARGS, "Tracking a bunch through a linac drift "},
		{"quad1",            wrap_linac_quad1,          METH_VARARGS, "Tracking a bunch through a linear part of a linac quad"},
		{"quad2",            wrap_linac_quad2,          METH_VARARGS, "Tracking a bunch through a non-linear part of a linac quad"},
		{"quad3",            wrap_linac_quad3,          METH_VARARGS, "Tracking a bunch through a non-linear part of a linac quad: long.field"},
		{"kick",             wrap_linac_kick,           METH_VARARGS, "Tracking a bunch through a kicker"},		
		{ NULL, NULL }
	};
	
	void initlinactracking()
	{
		PyObject *m;
		m = Py_InitModule((char*)"linac_tracking", linactrackingMethods);
	}
	
	PyObject* getLinacTrackingType(char* name)
	{
		PyObject* mod = PyImport_ImportModule("linac_tracking");
		PyObject* pyType = PyObject_GetAttrString(mod, name);
		Py_DECREF(mod);
		Py_DECREF(pyType);
		return pyType;
	}
  
#ifdef __cplusplus
}
#endif

} //end of namespace
