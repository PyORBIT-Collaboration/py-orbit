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






    //Integration through a very simple ring type RF cavity
    static PyObject* wrap_RingRF(PyObject *self, PyObject *args)
    {
        PyObject* pyBunch;
        double voltage, phase_s, ring_length;
        int harmonics_numb;
        int useCharge = 1;
        if(!PyArg_ParseTuple(	args, "Odidd|i:RingRF",
                             &pyBunch, &ring_length, &harmonics_numb,
                             &voltage, &phase_s, &useCharge))
        {
            error("errorbase - RingRF - cannot parse arguments!");
        }
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
        error_base::RingRF(cpp_bunch, ring_length, harmonics_numb,
                            voltage, phase_s, useCharge);
        Py_INCREF(Py_None);
        return Py_None;
    }

    static PyMethodDef errorbaseMethods[] =
    {
			{"rotatexy",         wrap_rotatexy,       METH_VARARGS, "Rotates bunch around z axis "},
			{"drift",            wrap_drift,          METH_VARARGS, "Tracking a bunch through a drift "},
			{"wrapbunch",		 wrap_wrapbunch,		  METH_VARARGS, "Tracking a bunch through a wrapbunch routine"},
			{"multp",            wrap_multp,          METH_VARARGS, "Tracking a bunch through a multipole "},
			{"multpfringeIN",    wrap_multpfringeIN,  METH_VARARGS, "Tracking a bunch through an IN edge of a multipole "},
			{"multpfringeOUT",   wrap_multpfringeOUT, METH_VARARGS, "Tracking a bunch through an OUT edge of a multipole"},
			{"kick",             wrap_kick,           METH_VARARGS, "Kicker element: chnges in x-prime, y-prime and dE"},
			{"quad1",            wrap_quad1,          METH_VARARGS, "Quadrupole element one: linear transport matrix "},
			{"quad2",            wrap_quad2,          METH_VARARGS, "Quadrupole element two: drift in quadrupole "},
			{"quadfringeIN",     wrap_quadfringeIN,   METH_VARARGS, "Quadrupole element IN edge"},
			{"quadfringeOUT",    wrap_quadfringeOUT,  METH_VARARGS, "Quadrupole element OUT edge"},
			{"wedgerotate",      wrap_wedgerotate,    METH_VARARGS, "Rotates coordinates by e for fringe fields at non-SBEND "},
			{"wedgedrift",       wrap_wedgedrift,     METH_VARARGS, "Drifts particles through wedge for non-SBEND "},
			{"wedgebend",        wrap_wedgebend,      METH_VARARGS, "Straight bends particles through wedge for non-SBEND "},
			{"bend1",            wrap_bend1,          METH_VARARGS, "Linear bend transport "},
			{"bend2",            wrap_bend2,          METH_VARARGS, "Kinetic bend transport (same as nonlinear quad transport - quad2) "},
			{"bend3",            wrap_bend3,          METH_VARARGS, "Nonlinear curvature bend transport depending on py and dE in Hamiltonian "},
			{"bend4",            wrap_bend4,          METH_VARARGS, "Nonlinear curvature bend transport depending on px in Hamiltonian "},
			{"bendfringeIN",     wrap_bendfringeIN,   METH_VARARGS, "Hard edge fringe field for a bend IN"},
			{"bendfringeOUT",    wrap_bendfringeOUT,  METH_VARARGS, "Hard edge fringe field for a bend OUT"},
			{"soln",             wrap_soln,           METH_VARARGS, "Integration through a solenoid "},
			{"wedgebendCF",      wrap_wedgebendCF,    METH_VARARGS, "Straight bends particles through wedge for Combined Function non-SBEND "},
			{"RingRF",           wrap_RingRF,         METH_VARARGS, "Tracking particles through a simple ring RF cavity."},
			{ NULL, NULL }
    };

    void initerrorbase(void)
    {
        PyObject *m, *d;
        m = Py_InitModule((char*)"error_base", errorbaseMethods);
        d = PyModule_GetDict(m);
        error_base::init_factorial();
        wrap_errorbase_matrix_generator::initMatrixGenerator(m);
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
