/* Program pair_interactions.

   Calculate the potential on mol1 due the presence of mol2, and
   visa-versa. Two calculations are performed in which the charges
   are zeroed out in mol1 but present in mol2, then zeroed out in
   mol2 but present in mol1.

 */


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>
using std::istringstream;
#include <string>
#include <math.h>
#include <list>
#include <algorithm>

#include "MEAD/PhysCond.h"
#include "MEAD/AtomSet.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/UniformElectrolyte.h"
#include "MEAD/DielByAtoms.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/FinDiffElstatPot.h"

// FIXME do this with STL functions instead
void
merge_atoms(AtomSet& merged, const AtomSet& a1, const AtomSet& a2)
{
  list<Atom> merged_atlist;
  for (AtomSet::const_iterator i = a1.begin(); i != a1.end(); ++i)
    merged_atlist.push_back(i->second);
  for (AtomSet::const_iterator i = a2.begin(); i != a2.end(); ++i)
    merged_atlist.push_back(i->second);
  merged = AtomSet(merged_atlist);
}

void
zero_charges(AtomSet& ats)
{
  for (AtomSet::iterator i = ats.begin(); i != ats.end(); ++i) {
    i->second.charge = 0.0;
  }
}

main(int argc, char* argv[])
{
  string molname1;
  string molname2;
  string progname(argv[0]);
  float epsin = 0; // This must be specified on the command line.  0 is error.

  // Parse the command line
  for (int iarg=1; iarg<argc; ++iarg) {
    if (iarg+2==argc) {  // These are the last ones, must be molnames
      molname1 = argv[iarg];
      ++iarg;
      molname2 = argv[iarg];
    }
    else if ((string) argv[iarg] == (string) "-epsin") {
      if (iarg+1>=argc) error (progname.c_str(), ": ERROR, not enough args after -epsin");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error (progname.c_str(), ": FAILED trying to read value after -epsin option");
      epsin = f;
    }
    else if ((string) argv[iarg] == (string) "-epsext") {
      if (iarg+1>=argc) error (progname.c_str(), ": ERROR, not enough args after -epsext");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error (progname.c_str(), ": FAILED trying to read value after -epsext option");
      PhysCond::set_epsext(f);
    }
    else if ((string) argv[iarg] == (string) "-solrad") {
      if (iarg+1>=argc) error (progname.c_str(), ": ERROR, not enough args after -solrad");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error (progname.c_str(), ": FAILED trying to read value after -solrad option");
      PhysCond::set_solrad(f);
    }
    else if ((string) argv[iarg] == (string) "-sterln") {
      if (iarg+1>=argc) error (progname.c_str(), ": ERROR, not enough args after -sterln");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error (progname.c_str(), ": FAILED trying to read value after -sterln option");
      PhysCond::set_sterln(f);
    }
    else if ((string) argv[iarg] == (string) "-ionicstr") {
      if (iarg+1>=argc) error (progname.c_str(), ": ERROR, not enough args after -ionicstr");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error (progname.c_str(), ": FAILED trying to read value after -ionicstr option");
      PhysCond::set_ionicstr(f);
    }
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
      cerr << progname << ": ERROR command option, " << argv[iarg]
	<< ", not recognized" << endl;
      error();
    }
  }

  float econv = PhysCond::get_econv();

  if (molname1 == (string) "")
    ::error(progname.c_str(), ": ERROR, molname1 not given on command line");
  if (molname1 == (string) "")
    ::error(progname.c_str(), ": ERROR, molname1 not given on command line");
  if (epsin == 0.0)
    ::error(progname.c_str(), ": ERROR, -epsin value not specified.",
	    "   This flag is now mandatory.");

  string m1m2 = molname1 + "_" + molname2;

  string interaction_1_2_filename = m1m2
    + "_interaction_1_2.dat";
  string interaction_2_1_filename = m1m2
    + "_interaction_2_1.dat";

  ofstream interaction_1_2_ost(interaction_1_2_filename.c_str());
  if (!interaction_1_2_ost) {
    ::error(progname.c_str(), ": ERROR, Cannot open output file interaction_1_2_filename");
  }
  ofstream interaction_2_1_ost(interaction_2_1_filename.c_str());
  if (!interaction_2_1_ost) {
    ::error(progname.c_str(), ": ERROR, Cannot open output file interaction_2_1_filename");
  }

  // Done with comand line processing
  // Print summary of input info from the command line:
  cout << "Starting " << argv[0] << " for pair of molecules named "
       << molname1 << " and " << molname2 << "\n";
  cout << "with interior dielectric constant = " << epsin << endl;
  cout << "using the following physical conditions:\n";
  PhysCond::print();
  if (blab3pt == &cout) cout << "Blab level set to " << 3 << endl;
  else if (blab2pt == &cout) cout << "Blab level set to " << 2 << endl;
  else if (blab1pt == &cout) cout << "Blab level set to " << 1 << endl;
  else cout << "No blab level set (so no blabbing)" << endl;

  AtomSet a1(molname1);
  a1.read();
  AtomSet a2(molname2);
  a2.read();

  FinDiffMethod fdm;
  string ogm_filename = molname1 + ".ogm";
  fdm.read(ogm_filename);

  // Resolve the lattice and do some checks.
  Coord geom_cent1 = a1.geom_cent();
  Coord geom_cent2 = a2.geom_cent();
  Coord geom_cent = geom_cent1 + (geom_cent2 - geom_cent1) / 2.0;
  Coord center_of_interest = geom_cent;
  fdm.resolve(geom_cent, center_of_interest);
  cout << "Using finite difference method with lattice levels:" << endl;
  cout << fdm << endl;

  UniformElectrolyte ely(PhysCond::get_ionicstr()); // No electrolyte is default

  // Create an atomset with mol1 charges zeroed out
  AtomSet a1zero(a1);
  zero_charges(a1zero);

  // Create an atomset with mol2 charges zeroed out
  AtomSet a2zero(a2);
  zero_charges(a2zero);

  // Calculate potential with zero mol1 charges
  AtomSet ats;
  merge_atoms(ats, a1zero, a2);
  TwoValueDielectricByAtoms *eps = new TwoValueDielectricByAtoms(ats, epsin);
  AtomChargeSet *acs = new AtomChargeSet(ats);
  FinDiffElstatPot phi1(fdm, eps, acs, &ely);
  phi1.solve();
  cout << "Calculation with molecule 1 charges zeroed is done" << endl;

  // Calculate potential with zero mol2 charges
  ats.clear();
  merge_atoms(ats, a1, a2zero);
  delete eps;
  eps = new TwoValueDielectricByAtoms(ats, epsin);
  delete acs;
  acs = new AtomChargeSet(ats);
  FinDiffElstatPot phi2(fdm, eps, acs, &ely);
  phi2.solve();
  cout << "Calculation with molecule 2 charges zeroed is done" << endl;

  delete eps;
  delete acs;

  // Write out interaction energy between atoms in molecule 1 with molecule 2
  //                          and between atoms in molecule 2 with molecule 1

  for (AtomSet::const_iterator ind = a1.begin(); ind != a1.end(); ++ind) {
    interaction_1_2_ost << ind->second.atname << " " << ind->second.resname << " " << ind->second.resnum << " " << ind->second.charge * phi1.value(ind->second.coord) * econv << endl;
  }

  for (AtomSet::const_iterator ind = a2.begin(); ind != a2.end(); ++ind) {
    interaction_2_1_ost << ind->second.atname << " " << ind->second.resname << " " << ind->second.resnum << " " << ind->second.charge * phi2.value(ind->second.coord) * econv << endl;
  }

  interaction_1_2_ost.close();
  interaction_2_1_ost.close();
}

// pair_interactions.cc ends here
