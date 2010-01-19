// MEAD.j -- interface for MEAD        -*- C++ -*-
// ID: $Id: MEAD.j,v 1.18 2001/07/13 20:39:05 bashford Exp $
// SOURCE: $Source: /cvs-repository/bashford/cvsroot/mead/swig/MEAD.j,v $

%module MEAD

// ------ Exception handling ---------------

%include exception.i

%{
#ifdef __cplusplus
#define BEGIN_CPLUSPLUS_SECTION	}
#else
#define BEGIN_CPLUSPLUS_SECTION
#endif

#ifdef __cplusplus
#define END_CPLUSPLUS_SECTION extern "C" {
#else
#define END_CPLUSPLUS_SECTION
#endif

#ifdef USE_EXCEPTIONS
#include "MEAD/MEADexcept.h"

#define MEAD_try_catch(func)					\
  try {								\
    func;							\
  }								\
  catch(MEADexcept e) {						\
    string error = e.get_error1();				\
    if (error == "") {						\
      error = "The MEAD library has thrown a MEADexcept";	\
    }								\
    else {							\
      error += (e.get_error2() + e.get_error3());		\
    }								\
    SWIG_exception(SWIG_RuntimeError, (char *) error.c_str());	\
    return NULL;						\
  }
#else
#define MEAD_try_catch(func) func
#endif

%}

%except(python) {
  MEAD_try_catch($function)
}

%pragma(python) include="MEADexcept.py"

// ------ Functions to emulate Python objects ---------------

%{
#include "./Python_funcs.h"
%}

// ------ Typemaps ---------------

%include ./meadtypes.i

//----------------------------------------------------------------------------
// Shadow Classes

%include Coord.j
%include ChargeDist.j
%include DielectricEnvironment.j
%include ElectrolyteEnvironment.j
%include ElstatPot.j
%include Atom.j
%include AtomID.j
%include AtomSet.j
%include PointCharge.j
%include AtomChargeSet.j
%include ManyPointCharge.j
%include OnePointCharge.j
%include UniformDielectric.j
%include DielectricSphere.j
%include DielectricSlab.j
%include DielByAtoms.j
%include DielMembAtoms.j
%include UniformElectrolyte.j
%include ElySphere.j
%include ElectrolyteByAtoms.j
%include SolvAccVol.j
%include SphericalHarmonic.j
%include CubeLatSpec.j
%include PhysCond.j
%include FinDiffMethod.j
%include AnalySphere.j
%include AnalySlab.j
%include Debye.j
%include FinDiffElstatPot.j
%include ElstatPotCombination.j
%include MomentAnalysis.j
