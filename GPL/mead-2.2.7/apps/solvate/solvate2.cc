/* Calculate the electrostatic solvation energy of a molecule.

    Copyright (c) 1993--1995 by Donald Bashford

    This source code file is part of the MEAD (Macroscopic
    Electrostatics with Atomic Detail) package of objects and
    programs.  MEAD is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 1, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; see the file COPYING.  If not, write to
    the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
    02139, USA.

    Donald Bashford can be contacted by electronic mail by the address,
    bashford@scripps.edu, or by paper mail at Department of Molecular
    Biology, The Scripps Research Institute, 10666 North Torrey Pines
    Road, La Jolla, California 92037.

$Id: solvate2.cc,v 2.1 2001/06/26 18:32:19 bergsma Exp $
*/
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <strstream>
using std::istrstream;
#include <fstream>
using std::ifstream;

#include <string>
using std::string;

#include "MEAD/PhysCond.h"
#include "MEAD/AtomSet.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/DielByAtoms.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ElectrolyteByAtoms.h"
#include "MEAD/UniformElectrolyte.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/FinDiffElstatPot.h"
#include "MEAD/FDGridLevel.h"

// For DEBUGGING
extern void big_rigid_stats();
extern void big_rigid_delete_all();
extern "C" {
  void malloc_stats();
}

main(int argc, char* argv[])
{
  string molname;
  float epsvac = 1.0;
  float doReactionField = 0;
  float epsin = 0; // This must be specified on the command line.  0 is error.
  // Parse the command line
  for (int iarg=1; iarg<argc; ++iarg) {
    if (iarg+1==argc) {  // This is the last one, must be molname
      molname = argv[iarg];
    }
    else if ((string) argv[iarg] == (string) "-epsin") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -epsin");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epsin option");
      epsin = f;
    }
    else if ((string) argv[iarg] == (string) "-epssol") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -epssol");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epssol option");
      PhysCond::set_epsext(f);
    }
    else if ((string) argv[iarg] == (string) "-epsvac") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -epsvac");
      istrstream ist(argv[++iarg]);
      ist >> epsvac;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epsvac option");
    }
    else if ((string) argv[iarg] == (string) "-solrad") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -solrad");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -solrad option");
      PhysCond::set_solrad(f);
    }
    else if ((string) argv[iarg] == (string) "-sterln") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -sterln");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -sterln option");
      PhysCond::set_sterln(f);
    }
    else if ((string) argv[iarg] == (string) "-ionicstr") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -ionicstr");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -ionicstr option");
      PhysCond::set_ionicstr(f);
    }
    else if ((string) argv[iarg] == (string) "-T") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -T");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -T option");
      PhysCond::set_T(f);
    }
    else if ((string) argv[iarg] == (string) "-kBolt") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -kBolt");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -kBolt option");
      PhysCond::set_kBolt(f);
    }
    else if ((string) argv[iarg] == (string) "-conconv") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -conconv");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -conconv option");
      PhysCond::set_conconv(f);
    }
    else if ((string) argv[iarg] == (string) "-econv") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -econv");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -econv option");
      PhysCond::set_econv(f);
    }
    else if ((string) argv[iarg] == (string) "-bohr_radius") {
      if (iarg+1>=argc)
	error ("main: ERROR, not enough args after -bohr_radius");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -bohr_radius option");
      PhysCond::set_bohr_radius(f);
    }
    else if ((string) argv[iarg] == (string) "-proton_charge") {
      if (iarg+1>=argc)
	error ("main: ERROR, not enough args after -proton_charge");
      istrstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -proton_charge opt.");
      PhysCond::set_proton_charge(f);
    }
    else if ((string) argv[iarg] == (string) "-ReactionField")
      doReactionField = 1;
    // Process -blab flags
    else if ((string) argv[iarg] == (string) "-epsave_oldway")
      FDGridLevel::set_epsave_oldway(true);
    else if ((string) argv[iarg] == (string) "-converge_oldway")
      FDGridLevel::set_use_fixed_maxrmsdiff(true);
    // Process -blab flags
    else if ((string) argv[iarg] == (string) "-blab1") {
      blab1pt = &cout;
    }
    else if ((string) argv[iarg] == (string) "-blab2") {
      blab1pt = &cout;
      blab2pt = &cout;
    }
    else if ((string) argv[iarg] == (string) "-blab3") {
      blab1pt = &cout;
      blab2pt = &cout;
      blab3pt = &cout;
    }
    else {
      cerr << "main: ERROR command option, " << argv[iarg]
	<< ", not recognized" << endl;
      error();
    }
  }

  if (molname == (string) "")
    ::error("main: ERROR, molname not given on command line");
  if (epsin == 0.0)
    ::error("main: ERROR, -epsin value not specified.",
	    "   This flag is now mandatory.");

  // Done with comand line processing
  // Print summary of input info from the command line:
  cout << "Starting " << argv[0] << " for molecule named " << molname << "\n";
  cout << "with interior dielectric constant = " << epsin << endl;
  cout << "using the following physical conditions:\n";
  PhysCond::print();
  cout << "and vacuum dielectric constant = " << epsvac << endl;
  if (blab3pt == &cout) cout << "Blab level set to " << 3 << endl;
  else if (blab2pt == &cout) cout << "Blab level set to " << 2 << endl;
  else if (blab1pt == &cout) cout << "Blab level set to " << 1 << endl;
  else cout << "No blab level set (so no blabbing)" << endl;

//  cout << "Malloc stats before reading atoms" << endl;
//  cout.flush();
//  cnull.flush();
//  malloc_stats();
//  big_rigid_stats();
//  {

  AtomSet a(molname);

/* For Debugging only..
  Coord crd(0.0, 0.0, 0.0);
  Atom at;
  at.atname = "CL";
  at.resname = "CL";
  at.resnum = 1;
  at.coord = crd;
  at.charge = 1.0;
  at.rad = 2.0;
  AtomID key(at.resnum, at.atname);
  a.insert(at);
  //a[key] = at;
*/
  a.read();


  AtomChargeSet rho(a);

  TwoValueDielectricByAtoms eps(a, epsin);

  ElectrolyteByAtoms ely(a);

  FinDiffMethod fdm;

  string ogm_filename = molname + ".ogm";
  fdm.read(ogm_filename);

/* Instead of reading in
  fdm.add_level(41, 1.0, crd);	// the coarsest level
  fdm.add_level(41, 0.25, crd);	// a finer one in a smaller region.
*/

  Coord interesting (0,0,0);
  fdm.resolve(a.geom_cent(), interesting);

  cout << "Using finite difference method with lattice levels:" << endl;
  cout << fdm;

  float prod_sol;

  FinDiffElstatPot phi(fdm, &eps, &rho, &ely);
  phi.solve();
  prod_sol = phi * rho;
  cout << "prod_sol = " << prod_sol << endl;



  float safe_eps_sol = PhysCond::get_epsext(); // What for?
  PhysCond::set_epsext(epsvac);
  PhysCond::set_ionicstr(0.0);

  UniformElectrolyte elyvac(0.0);  // No electrolyte is the default
  TwoValueDielectricByAtoms vac_eps(a, epsin);
  FinDiffElstatPot vac_phi(fdm, &vac_eps, &rho, &elyvac);
  vac_phi.solve();
  float prod_vac = vac_phi * rho;
  cout << "prod_vac = " << prod_vac << endl;
  float solvation_energy = (prod_sol - prod_vac) / 2 * PhysCond::get_econv();
  cout << "\n\nSOLVATION ENERGY = " << solvation_energy
    << "\n(probably in kcal/mole)" << endl;



  if (doReactionField) {
    string fieldpoint_filename = molname + ".fpt" ;
    ifstream fpt (fieldpoint_filename.c_str());
    if (!fpt.good())
      ::error(argv[0], "main: failed to open for reading, ",
	      fieldpoint_filename.c_str());
    string reacfield_filename = molname + ".rf";
    ofstream rf(reacfield_filename.c_str());
    if (!rf.good())
      ::error(argv[0], "main: failed to open for writing, ",
	      reacfield_filename.c_str());
    while (fpt.good()) {
      Coord c;
      fpt >> c;
      if (fpt.eof()) break;
      if (!fpt.good())
	::error(argv[0], "main: input failure while reading ",
		fieldpoint_filename.c_str());
      rf << phi.value(c) - vac_phi.value(c) << "\n";
      if (!rf.good())
	::error(argv[0], "main: output failure while writing ",
		reacfield_filename.c_str());
    }
  }


//  }
//  big_rigid_delete_all();
//  cout << "Malloc stats at end of Main" << endl;
//  cout.flush();
//  cnull.flush();
//  malloc_stats();
//  big_rigid_stats();
//  return 0;

}

// solvate.cc ends here
