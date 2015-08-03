#include "Python.h"
#include "orbit_mpi.hh"

#include "pyORBIT_Object.hh"

#include "errorbase.hh"

#include "wrap_errorbase.hh"

namespace wrap_errorbase
{
  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C"
{
#endif

  //---------------------------------------------------------
  // errorbase method wrappers
  //---------------------------------------------------------

  // Displace the coordinates of a bunch
  static PyObject* wrap_CoordDisplacement(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double dx, dxp, dy, dyp, dz, dE;
    if(!PyArg_ParseTuple(args, "Odddddd:CoordDisplacement",
                         &pyBunch, &dx, &dxp, &dy, &dyp, &dz, &dE))
    {
      error("errorbase - CoordDisplacement - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::CoordDisplacement(cpp_bunch, dx, dxp, dy, dyp, dz, dE);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Longitudinally displace a bunch
  static PyObject* wrap_LongDisplacement(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double ds;
    if(!PyArg_ParseTuple(args, "Od:LongDisplacement",
                         &pyBunch, &ds))
    {
      error("errorbase - LongDisplacement - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::LongDisplacement(cpp_bunch, ds);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // XY rotate a bunch
  static PyObject* wrap_StraightRotationXY(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double anglexy;
    if(!PyArg_ParseTuple(args, "Od:StraightRotationXY",
                         &pyBunch, &anglexy))
    {
      error("errorbase - StraightRotationXY - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::StraightRotationXY(cpp_bunch, anglexy);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // XS rotate a bunch entering element
  static PyObject* wrap_StraightRotationXSI(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double anglexsi, length;
    if(!PyArg_ParseTuple(args, "Odd:StraightRotationXSI",
                         &pyBunch, &anglexsi, &length))
    {
      error("errorbase - StraightRotationXSI - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::StraightRotationXSI(cpp_bunch, anglexsi, length);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // XS rotate a bunch leaving element
  static PyObject* wrap_StraightRotationXSF(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double anglexsf, length;
    if(!PyArg_ParseTuple(args, "Odd:StraightRotationXSF",
                         &pyBunch, &anglexsf, &length))
    {
      error("errorbase - StraightRotationXSF - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::StraightRotationXSF(cpp_bunch, anglexsf, length);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // YS rotate a bunch entering element
  static PyObject* wrap_StraightRotationYSI(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double angleysi, length;
    if(!PyArg_ParseTuple(args, "Odd:StraightRotationYSI",
                         &pyBunch, &angleysi, &length))
    {
      error("errorbase - StraightRotationYSI - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::StraightRotationYSI(cpp_bunch, angleysi, length);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // YS rotate a bunch leaving element
  static PyObject* wrap_StraightRotationYSF(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double angleysf, length;
    if(!PyArg_ParseTuple(args, "Odd:StraightRotationYSF",
                         &pyBunch, &angleysf, &length))
    {
      error("errorbase - StraightRotationYSF - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::StraightRotationYSF(cpp_bunch, angleysf, length);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Bend field strength error to a bunch entering element
  static PyObject* wrap_BendFieldI(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double drho;
    if(!PyArg_ParseTuple(args, "Od:BendFieldI",
                         &pyBunch, &drho))
    {
      error("errorbase - BendFieldI - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::BendFieldI(cpp_bunch, drho);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Bend field strength error to a bunch leaving element
  static PyObject* wrap_BendFieldF(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double drho;
    if(!PyArg_ParseTuple(args, "Od:BendFieldF",
                         &pyBunch, &drho))
    {
      error("errorbase - BendFieldF - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::BendFieldF(cpp_bunch, drho);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // X displacement error to a bunch entering bend
  static PyObject* wrap_BendDisplacementXI(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double anglexi, disp;
    if(!PyArg_ParseTuple(args, "Odd:BendDisplacementXI",
                         &pyBunch, &anglexi, &disp))
    {
      error("errorbase - BendDisplacementXI - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::BendDisplacementXI(cpp_bunch, anglexi, disp);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // X displacement error to a bunch leaving bend
  static PyObject* wrap_BendDisplacementXF(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double anglexf, disp;
    if(!PyArg_ParseTuple(args, "Odd:BendDisplacementXF",
                         &pyBunch, &anglexf, &disp))
    {
      error("errorbase - BendDisplacementXF - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::BendDisplacementXF(cpp_bunch, anglexf, disp);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Y displacement error to a bunch entering bend
  static PyObject* wrap_BendDisplacementYI(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double disp;
    if(!PyArg_ParseTuple(args, "Od:BendDisplacementYI",
                         &pyBunch, &disp))
    {
      error("errorbase - BendDisplacementYI - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::BendDisplacementYI(cpp_bunch, disp);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Y displacement error to a bunch leaving bend
  static PyObject* wrap_BendDisplacementYF(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double disp;
    if(!PyArg_ParseTuple(args, "Od:BendDisplacementYF",
                         &pyBunch, &disp))
    {
      error("errorbase - BendDisplacementYF - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::BendDisplacementYF(cpp_bunch, disp);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // L displacement error to a bunch entering bend
  static PyObject* wrap_BendDisplacementLI(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double angleli, disp;
    if(!PyArg_ParseTuple(args, "Odd:BendDisplacementLI",
                         &pyBunch, &angleli, &disp))
    {
      error("errorbase - BendDisplacementLI - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::BendDisplacementLI(cpp_bunch, angleli, disp);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // L displacement error to a bunch leaving bend
  static PyObject* wrap_BendDisplacementLF(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double anglelf, disp;
    if(!PyArg_ParseTuple(args, "Odd:BendDisplacementLF",
                         &pyBunch, &anglelf, &disp))
    {
      error("errorbase - BendDisplacementLF - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::BendDisplacementLF(cpp_bunch, anglelf, disp);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // General rotation error to a bunch entering element
  static PyObject* wrap_RotationI(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double anglei, rhoi, theta, length;
    const char* et   = NULL;
    const char* type = NULL;
    if(!PyArg_ParseTuple(args, "Oddddss:RotationI",
                         &pyBunch, &anglei, &rhoi, &theta,
                         &length, &et, &type))
    {
      error("errorbase - RotationI - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::RotationI(cpp_bunch, anglei, rhoi, theta,
                          length, et, type);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // General rotation error to a bunch leaving element
  static PyObject* wrap_RotationF(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double anglef, rhoi, theta, length;
    const char* et   = NULL;
    const char* type = NULL;
    if(!PyArg_ParseTuple(args, "Odddd:RotationF",
                         &pyBunch, &anglef, &rhoi, &theta,
                         &length, &et, &type))
    {
      error("errorbase - RotationF - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::RotationF(cpp_bunch, anglef, rhoi, theta,
                          length, et, type);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Oscillating dipole kick a bunch
  static PyObject* wrap_DipoleKickerOsc(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double k, phaselength, phase;
    if(!PyArg_ParseTuple(args, "Oddd:DipoleKickerOsc",
                         &pyBunch, &k, &phaselength, &phase))
    {
      error("errorbase - DipoleKickerOsc - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::DipoleKickerOsc(cpp_bunch, k, phaselength, phase);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Quadrupole kick a bunch
  static PyObject* wrap_QuadKicker(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double k;
    if(!PyArg_ParseTuple(args, "Od:QuadKicker",
                         &pyBunch, &k))
    {
      error("errorbase - QuadKicker - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::QuadKicker(cpp_bunch, k);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Oscillating quadrupole kick a bunch
  static PyObject* wrap_QuadKickerOsc(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    double k, phaselength, phase;
    if(!PyArg_ParseTuple(args, "Oddd:QuadKickerOsc",
                         &pyBunch, &k, &phaselength, &phase))
    {
      error("errorbase - QuadKickerOsc - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::QuadKickerOsc(cpp_bunch, k, phaselength, phase);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Drift a particle
  static PyObject* wrap_drifti(PyObject *self, PyObject *args)
  {
    PyObject* pyBunch;
    int i;
    double length;
    if(!PyArg_ParseTuple(args, "Oid:drifti",
                         &pyBunch, &i, &length))
    {
      error("errorbase - drifti - cannot parse arguments!");
    }
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
    error_base::drifti(cpp_bunch, i, length);
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Error function
  static PyObject* wrap_derf(PyObject *self, PyObject *args)
  {
    double x;
    if(!PyArg_ParseTuple(args, "d:derf", &x))
    {
      error("errorbase - derf - cannot parse arguments!");
    }
    double errf = error_base::derf(x);
    return Py_BuildValue("d", errf);
  }

  // Helps find Gaussian distribution
  static PyObject* wrap_root_normal(PyObject *self, PyObject *args)
  {
    double errtest, ymin, ymax, tol;
    if(!PyArg_ParseTuple(args, "dddd:root_normal",
                         &errtest, &ymin, &ymax, &tol))
    {
      error("errorbase - root_normal - cannot parse arguments!");
    }
    double root = error_base::root_normal(errtest, ymin, ymax, tol);
    return Py_BuildValue("d", root);
  }

  // Returns Gaussian distribution
  static PyObject* wrap_getGauss(PyObject *self, PyObject *args)
  {
    double mean, sigma, cutoff;
    if(!PyArg_ParseTuple(args, "ddd:getGauss",
                         &mean, &sigma, &cutoff))
    {
      error("errorbase - getGauss - cannot parse arguments!");
    }
    double sample = error_base::getGauss(mean, sigma, cutoff);
    return Py_BuildValue("d", sample);
  }

  static PyMethodDef errorbaseMethods[] =
  {
    {"CoordDisplacement",   wrap_CoordDisplacement,   METH_VARARGS, "Displace the coordinates of a bunch"},
    {"LongDisplacement",    wrap_LongDisplacement,    METH_VARARGS, "Longitudinally displace a bunch"},
    {"StraightRotationXY",  wrap_StraightRotationXY,  METH_VARARGS, "XY rotate a bunch"},
    {"StraightRotationXSI", wrap_StraightRotationXSI, METH_VARARGS, "XS rotate a bunch entering element"},
    {"StraightRotationXSF", wrap_StraightRotationXSF, METH_VARARGS, "XS rotate a bunch leaving element"},
    {"StraightRotationYSI", wrap_StraightRotationYSI, METH_VARARGS, "YS rotate a bunch entering element"},
    {"StraightRotationYSF", wrap_StraightRotationYSF, METH_VARARGS, "YS rotate a bunch leaving element"},
    {"BendFieldI",          wrap_BendFieldI,          METH_VARARGS, "Bend field strength error to a bunch entering element"},
    {"BendFieldF",          wrap_BendFieldF,          METH_VARARGS, "Bend field strength error to a bunch leaving element"},
    {"BendDisplacementXI",  wrap_BendDisplacementXI,  METH_VARARGS, "X displacement error to a bunch entering bend"},
    {"BendDisplacementXF",  wrap_BendDisplacementXF,  METH_VARARGS, "X displacement error to a bunch leaving bend"},
    {"BendDisplacementYI",  wrap_BendDisplacementYI,  METH_VARARGS, "Y displacement error to a bunch entering bend"},
    {"BendDisplacementYF",  wrap_BendDisplacementYF,  METH_VARARGS, "Y displacement error to a bunch leaving bend"},
    {"BendDisplacementLI",  wrap_BendDisplacementLI,  METH_VARARGS, "L displacement error to a bunch entering bend"},
    {"BendDisplacementLF",  wrap_BendDisplacementLF,  METH_VARARGS, "L displacement error to a bunch leaving bend"},
    {"RotationI",           wrap_RotationI,           METH_VARARGS, "General rotation error to a bunch entering element"},
    {"RotationF",           wrap_RotationF,           METH_VARARGS, "General rotation error to a bunch leaving element"},
    {"DipoleKickerOsc",     wrap_DipoleKickerOsc,     METH_VARARGS, "Oscillating dipole kick a bunch"},
    {"QuadKicker",          wrap_QuadKicker,          METH_VARARGS, "Quadrupole kick a bunch"},
    {"QuadKickerOsc",       wrap_QuadKickerOsc,       METH_VARARGS, "Oscillating quadrupole kick a bunch"},
    {"drifti",              wrap_drifti,              METH_VARARGS, "Drifts a macroparticle"},
    {"derf",                wrap_derf,                METH_VARARGS, "Error function"},
    {"root_normal",         wrap_root_normal,         METH_VARARGS, "Helps find Gaussian distribution"},
    {"getGauss",            wrap_getGauss ,           METH_VARARGS, "Returns Gaussian distribution"},
    { NULL, NULL }
  };

  void initerrorbase(void)
  {
    PyObject *m, *d;
    m = Py_InitModule((char*)"error_base", errorbaseMethods);
    d = PyModule_GetDict(m);
  }

  PyObject* getBaseERRORType(char* name)
  {
    PyObject* mod = PyImport_ImportModule("error_base");
    PyObject* pyType = PyObject_GetAttrString(mod, name);
    Py_DECREF(mod);
    Py_DECREF(pyType);
    return pyType;
  }
	
#ifdef __cplusplus
}
#endif

} //end of namespace wrap_errorbase
