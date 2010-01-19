/* Program to calculate potential due to molecule and write to file.

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

$Id: potential.cc,v 2.17 2007/11/26 22:38:22 bashford Exp $
*/
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>
using std::istringstream;
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
#include "MEAD/FinDiffMethod.h"
#include "MEAD/FinDiffElstatPot.h"
#include "MEAD/FDGridLevel.h"

main(int argc, char* argv[])
{
  string molname;
  string initfield;
  string outfield;
  float epsin = 0; // This must be specified on the command line.  0 is error.
  // Parse the command line
  for (int iarg=1; iarg<argc; ++iarg) {
    if (iarg+1==argc) {  // This is the last one, must be molname
      molname = argv[iarg];
    }
    else if ((string) argv[iarg] == (string) "-CoarseFieldInit") {
      if (iarg+1>=argc)
	error("main: ERROR, not enough args after -CoarseFieldInit");
      initfield = argv[++iarg];
    }
    else if ((string) argv[iarg] == (string) "-CoarseFieldOut") {
      if (iarg+1>=argc)
	error("main: ERROR, not enough args after -CoarseFieldOut");
      outfield = argv[++iarg];
    }
    else if ((string) argv[iarg] == (string) "-AvsScaleFactor") {
      if (iarg+1>=argc)
	error("main: ERROR, not enough args after -CoarseFieldOut");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epsin option");
      FDGridLevel::set_avsScaleFactor(f);
    }
    else if ((string) argv[iarg] == (string) "-epsin") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -epsin");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epsin option");
      epsin = f;
    }
    else if ((string) argv[iarg] == (string) "-epsext") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -epsext");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epsext option");
      PhysCond::set_epsext(f);
    }
    else if ((string) argv[iarg] == (string) "-solrad") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -solrad");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -solrad option");
      PhysCond::set_solrad(f);
    }
    else if ((string) argv[iarg] == (string) "-sterln") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -sterln");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -sterln option");
      PhysCond::set_sterln(f);
    }
    else if ((string) argv[iarg] == (string) "-ionicstr") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -ionicstr");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -ionicstr option");
      PhysCond::set_ionicstr(f);
    }
    else if ((string) argv[iarg] == (string) "-T") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -T");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -T option");
      PhysCond::set_T(f);
    }
    else if ((string) argv[iarg] == (string) "-kBolt") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -kBolt");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -kBolt option");
      PhysCond::set_kBolt(f);
    }
    else if ((string) argv[iarg] == (string) "-conconv") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -conconv");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -conconv option");
      PhysCond::set_conconv(f);
    }
    else if ((string) argv[iarg] == (string) "-econv") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -econv");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -econv option");
      PhysCond::set_econv(f);
    }
    else if ((string) argv[iarg] == (string) "-bohr_radius") {
      if (iarg+1>=argc)
	error ("main: ERROR, not enough args after -bohr_radius");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -bohr_radius option");
      PhysCond::set_bohr_radius(f);
    }
    else if ((string) argv[iarg] == (string) "-proton_charge") {
      if (iarg+1>=argc)
	error ("main: ERROR, not enough args after -proton_charge");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -proton_charge opt.");
      PhysCond::set_proton_charge(f);
    }
    else if ((string) argv[iarg] == (string) "-epsave_oldway")
      FDGridLevel::set_epsave_oldway(true);
    else if ((string) argv[iarg] == (string) "-converge_oldway")
      FDGridLevel::set_use_fixed_maxrmsdiff(true);
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
  if (blab3pt == &cout) cout << "Blab level set to " << 3 << endl;
  else if (blab2pt == &cout) cout << "Blab level set to " << 2 << endl;
  else if (blab1pt == &cout) cout << "Blab level set to " << 1 << endl;
  else cout << "No blab level set (so no blabbing)" << endl;

  AtomSet a(molname);
  a.read();
  ChargeDist_lett* prho = new AtomChargeSet(a);
  DielectricEnvironment_lett* peps = new TwoValueDielectricByAtoms(a, epsin);
  ElectrolyteEnvironment_lett* pely = new ElectrolyteByAtoms(a);


  FinDiffMethod fdm;
  string ogm_filename = molname + ".ogm";
  fdm.read(ogm_filename);
  Coord interesting (0,0,0);
  fdm.resolve(a.geom_cent(), interesting);
  cout << "Using finite difference method with lattice levels:" << endl;
  cout << fdm;

  FinDiffElstatPot phi(fdm, peps, prho, pely);
  if (initfield.length()) {
    string ffn = initfield;
    ffn += ".fld"; 
    const char *ffnc = ffn.c_str();
    // open it just as an existence test
    ifstream initfield_file_existence(ffnc);
    if (initfield_file_existence.good()) {
      initfield_file_existence.close();
      phi.solve_using_coarse_init(initfield);
    }
    else {
      cerr << "ERROR from " << argv[0] << " main program:\n"
	   << "could not open CourseFieldInit file, " << ffn << endl;
      return 1;
    }
  }
  else
    phi.solve();

  if (outfield.length() != 0)
    phi.write_coarse_field(outfield);

  string fieldpoint_filename_string = molname + ".fpt" ;
  const char *fieldpoint_filename = fieldpoint_filename_string.c_str();
  ifstream fpt (fieldpoint_filename, std::ios::in);
  if (!fpt.good() && outfield.length() == 0)
    cerr << "WARNING from " << argv[0] << " main program:\n"
	 << "Could not open field point file, " << fieldpoint_filename
	 << ", for reading,\n and no CoarseFieldOut was specified\n"
	 << "Exiting without giving any potentials." << endl;
  else {
    while(fpt.good()) {
      Coord c;
      fpt >> c;
      if (fpt.eof()) break;
      if (!fpt.good())
	::error(argv[0], "main: input failure while reading ",
		fieldpoint_filename);
      cout << phi.value(c) << "\n";
    }
  }
}

// potential.cc ends here
