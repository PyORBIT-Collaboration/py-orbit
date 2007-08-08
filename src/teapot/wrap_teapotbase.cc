#include "Python.h"
#include "orbit_mpi.hh"

#include "teapotbase.hh"

#include "wrap_teapotbase.hh"

namespace wrap_teapotbase{

     void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

		//---------------------------------------------------------
		//teapotbase methods wrappers
		//---------------------------------------------------------

		//Rotate bunch around z axis
    static PyObject* wrap_rotatexy(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double anglexy;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:rotatexy",&pyBunch,&anglexy)){
				error("teapotbase - rotatexy - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::rotatexy(cpp_bunch,anglexy);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Tracking a bunch through a drift
    static PyObject* wrap_drift(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double length;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:drift",&pyBunch,&length)){
				error("teapotbase - drift - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::drift(cpp_bunch,length);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}


		//Tracking a bunch through a multipole
    static PyObject* wrap_multp(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
      int pole, skew;
			double kl;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Oidi:multp",&pyBunch,&pole,&kl,&skew)){
				error("teapotbase - multp - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::multp(cpp_bunch, pole, kl, skew);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Tracking a bunch through an IN edge of a multipole
    static PyObject* wrap_multpfringeIN(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
      int pole, skew;
			double kl;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Oidi:multpfringeIN",&pyBunch,&pole,&kl,&skew)){
				error("teapotbase - multpfringeIN - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::multpfringeIN(cpp_bunch, pole, kl, skew);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Tracking a bunch through an OUT edge of a multipole
    static PyObject* wrap_multpfringeOUT(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
      int pole, skew;
			double kl;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Oidi:multpfringeOUT",&pyBunch,&pole,&kl,&skew)){
				error("teapotbase - multpfringeOUT - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::multpfringeOUT(cpp_bunch, pole, kl, skew);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Kicker element
    static PyObject* wrap_kick(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double kx, ky, kE;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Oddd:kick",&pyBunch,&kx,&ky,&kE)){
				error("teapotbase - kick - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::kick(cpp_bunch,kx,ky,kE);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Quadrupole element one: linear transport matrix
    static PyObject* wrap_quad1(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double length, kq;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Odd:quad1",&pyBunch,&length,&kq)){
				error("teapotbase - quad1 - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::quad1(cpp_bunch,length,kq);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Quadrupole element two: drift in quadrupole
    static PyObject* wrap_quad2(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double length;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:quad2",&pyBunch,&length)){
				error("teapotbase - quad2 - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::quad2(cpp_bunch,length);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Quadrupole element IN edge
    static PyObject* wrap_quadfringeIN(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double kq;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:quadfringeIN",&pyBunch,&kq)){
				error("teapotbase - quadfringeIN - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::quadfringeIN(cpp_bunch,kq);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Quadrupole element OUT edge
    static PyObject* wrap_quadfringeOUT(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double kq;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:quadfringeOUT",&pyBunch,&kq)){
				error("teapotbase - quadfringeOUT - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::quadfringeOUT(cpp_bunch,kq);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Rotates coordinates by e for fringe fields at non-SBEND
    static PyObject* wrap_wedgerotate(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double e;
			int frinout;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Odi:wedgerotate",&pyBunch,&e,&frinout)){
				error("teapotbase - wedgerotate - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::wedgerotate(cpp_bunch,e,frinout);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Drifts particles through wedge for non-SBEND
    static PyObject* wrap_wedgedrift(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double e;
			int inout;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Odi:wedgedrift",&pyBunch,&e,&inout)){
				error("teapotbase - wedgedrift - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::wedgedrift(cpp_bunch,e,inout);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Straight bends particles through wedge for non-SBEND
    static PyObject* wrap_wedgebend(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double e, rho;
			int inout, nsteps;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Odidi:wedgebend",&pyBunch,&e,&inout,&rho,&nsteps)){
				error("teapotbase - wedgebend - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::wedgebend(cpp_bunch, e, inout, rho, nsteps);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Linear bend transport
    static PyObject* wrap_bend1(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double length;
			double th;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Odd:bend1",&pyBunch,&length,&th)){
				error("teapotbase - bend1 - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::bend1(cpp_bunch,length,th);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Kinetic bend transport (same as nonlinear quad transport - quad2)
    static PyObject* wrap_bend2(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double length;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:bend2",&pyBunch,&length)){
				error("teapotbase - bend2 - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::bend2(cpp_bunch,length);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

	  //Nonlinear curvature bend transport depending on py and dE in Hamiltonian
    static PyObject* wrap_bend3(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double th;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:bend3",&pyBunch,&th)){
				error("teapotbase - bend3 - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::bend3(cpp_bunch,th);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Nonlinear curvature bend transport depending on px in Hamiltonian
    static PyObject* wrap_bend4(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double th;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:bend4",&pyBunch,&th)){
				error("teapotbase - bend4 - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::bend4(cpp_bunch,th);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Hard edge fringe field for a bend IN
    static PyObject* wrap_bendfringeIN(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double rho;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:bendfringeIN",&pyBunch,&rho)){
				error("teapotbase - bendfringeIN - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::bendfringeIN(cpp_bunch,rho);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Hard edge fringe field for a bend OUT
    static PyObject* wrap_bendfringeOUT(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double rho;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:bendfringeOUT",&pyBunch,&rho)){
				error("teapotbase - bendfringeOUT - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::bendfringeOUT(cpp_bunch,rho);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Integration through a solenoid
    static PyObject* wrap_soln(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double length, B;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Odd:soln",&pyBunch,&length,&B)){
				error("teapotbase - soln - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::soln(cpp_bunch,length,B);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Hard edge fringe field for a solenoid IN
    static PyObject* wrap_solnfringeIN(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double B;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:solnfringeIN",&pyBunch,&B)){
				error("teapotbase - solnfringeIN - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::solnfringeIN(cpp_bunch,B);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Hard edge fringe field for a solenoid OUT
    static PyObject* wrap_solnfringeOUT(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double B;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Od:solnfringeOUT",&pyBunch,&B)){
				error("teapotbase - solnfringeOUT - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::solnfringeOUT(cpp_bunch,B);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

		//Straight bends particles through wedge for Combined Function non-SBEND
    static PyObject* wrap_wedgebendCF(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double e, rho;
			int inout, vecnum, nsteps;

			std::vector<int> poleV;
			std::vector<double> klV;
			std::vector<int> skewV;

			PyObject* polePySeq;
			PyObject* klPySeq;
			PyObject* skewPySeq;

			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"OdidiOOOi:wedgebendCF",&pyBunch,
					&e,&inout,&rho,&vecnum,&polePySeq,&klPySeq,&skewPySeq,&nsteps)){
				error("teapotbase - wedgebendCF - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			//unpack all sequences
			if(!PySequence_Check(polePySeq) ||
				 !PySequence_Check(klPySeq) ||
				 !PySequence_Check(skewPySeq)){
			    error("teapotbase - wedgebendCF - sequences with poles,kls or skews are wrong!");
				}

				if( (PySequence_Size(polePySeq) != vecnum) ||
					(PySequence_Size(klPySeq)   != vecnum) ||
					(PySequence_Size(skewPySeq) != vecnum)){
					error("teapotbase - wedgebendCF - size of sequences with poles,kls or skews are wrong!");
				}

				for(int i = 0; i < vecnum; i++){
					PyObject* polePy = PySequence_GetItem(polePySeq,i);
					PyObject* klPy = PySequence_GetItem(klPySeq,i);
					PyObject* skewPy = PySequence_GetItem(skewPySeq,i);

					PyObject* tuplePy = Py_BuildValue("(OOO)",polePy,klPy,skewPy);

					double kl;
					int pole, skew;

					if(!PyArg_ParseTuple(tuplePy,"idi",&pole,&kl,&skew)){
						error("teapotbase - wedgebendCF - values in poles,kls or skews are wrong!");
					}

					poleV.push_back(pole);
					klV.push_back(kl);
					skewV.push_back(skew);

					Py_DECREF(tuplePy);

					Py_DECREF(polePy);
					Py_DECREF(klPy);
					Py_DECREF(skewPy);
				}

				//the teapotbase function call
			  teapot_base::wedgebendCF(cpp_bunch,e,inout,rho,vecnum,poleV,klV,skewV,nsteps);

			  //clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			  Py_DECREF(py_bunch_ref);

			  Py_INCREF(Py_None);
			  return Py_None;
		  }

		//Integration through a very simple ring type RF cavity
    static PyObject* wrap_ringRF(PyObject *self, PyObject *args) {
			PyObject* pyBunch;
			double voltage, phase_s, ring_length;
			int harmonics_numb;
			//NO NEW OBJECT CREATED BY PyArg_ParseTuple! - NO NEED OF Py_DECREF()
			if(!PyArg_ParseTuple(	args,"Odidd:ringRF",&pyBunch,&ring_length,&harmonics_numb,&voltage,&phase_s)){
				error("teapotbase - ringRF - cannot parse arguments!");
			}

			PyObject* py_bunch_ref = PyObject_GetAttrString( pyBunch ,"cpp_ptr");
			Bunch* cpp_bunch = (Bunch*) PyCObject_AsVoidPtr(py_bunch_ref);

			teapot_base::ringRF(cpp_bunch,ring_length,harmonics_numb,voltage,phase_s);

			//clear the reference from PyObject_GetAttrString( pyBunch ,"cpp_ptr")
			Py_DECREF(py_bunch_ref);

			Py_INCREF(Py_None);
			return Py_None;
		}

  static PyMethodDef teapotbaseMethods[] = {
    {"rotatexy",         wrap_rotatexy,       METH_VARARGS, "Rotates bunch around z axis "},
    {"drift",            wrap_drift,          METH_VARARGS, "Tracking a bunch through a drift "},
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
    {"solnfringeIN",     wrap_solnfringeIN,   METH_VARARGS, "Hard edge fringe field for a solenoid IN "},
    {"solnfringeOUT",    wrap_solnfringeOUT,  METH_VARARGS, "Hard edge fringe field for a solenoid OUT "},
    {"wedgebendCF",      wrap_wedgebendCF,    METH_VARARGS, "Straight bends particles through wedge for Combined Function non-SBEND "},
    {"ringRF",           wrap_ringRF,         METH_VARARGS, "Tracking particles through a simple ring RF cavity."},
    { NULL, NULL }
  };


  void initteapotbase(void) {
    PyObject *m, *d;
    m = Py_InitModule((char*)"teapot_base",teapotbaseMethods);
    d = PyModule_GetDict(m);

		teapot_base::init_factorial();
  }


#ifdef __cplusplus
}
#endif

} //end of namespace
