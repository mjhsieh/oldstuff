/* Titration calculations for multi-site molecule, option for flexibility

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

$Id: multiflex.cc,v 2.18 2004/11/15 05:30:17 bashford Exp $
*/

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <sstream>
using std::istringstream;
#include <string>
using std::string;

#include "MultiFlexiSiteMol.h"
#include "MEAD/AtomSet.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/PhysCond.h"
#include "FlexiSite.h"
#include "MEAD/Bigmem.h"
#include "MEAD/DielByAtoms.h"
#include "MEAD/ElectrolyteByAtoms.h"
#include "MEAD/FDGridLevel.h"
#include "SiteUtil.h"
#include "MEAD/SolvAccVol.h"
#include "MEAD/UniformDielectric.h"
#include "MEAD/UniformElectrolyte.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/Debye.h"
#include "MEAD/PairGenBorn.h"

class PGBMaker : public AtomicElstatMaker {
public:
  PGBMaker(float eps_in) : epsin(eps_in) {}
  ElstatPot operator() (const AtomSet& ats, const AtomChargeSet& gen,
			const Coord& cen_intr);
private:
  float epsin;
};

ElstatPot PGBMaker::operator() (const AtomSet& ats, const AtomChargeSet& gen,
				const Coord& cen_intr)
{
  AtomChargeSet*  acs = new AtomChargeSet(gen);
  TwoValueDielectricByAtoms* tvdba = new TwoValueDielectricByAtoms(ats,epsin);
  UniformElectrolyte* ue = new  UniformElectrolyte;
  ElstatPot_lett *pgb = new PairGenBorn(tvdba, acs, ue);
  return ElstatPot(pgb);
}

class CoulMaker : public AtomicElstatMaker {
public:
  CoulMaker(float eps_in) : epsin(eps_in) {}
  ElstatPot operator() (const AtomSet& ats, const AtomChargeSet& gen,
			const Coord& cen_intr);
private:
  float epsin;
};

ElstatPot CoulMaker::operator() (const AtomSet& ats, const AtomChargeSet& gen,
				const Coord& cen_intr)
{
  ChargeDist_lett*  acs = new AtomChargeSet(gen);
  UniformDielectric* ud = new UniformDielectric(epsin);
  UniformElectrolyte* ue = new  UniformElectrolyte(0.0);
  return ElstatPot(new Debye(ud, acs, ue));
}

class FD2DielEMaker : public AtomicElstatMaker {
public:
  FD2DielEMaker(AtomSet& ats, float eps_in, const FinDiffMethod& fm)
    : atlist(ats), epsin(eps_in), fdm(fm), has_eps_ely(false) {}
  bool atsdiffer(const AtomSet& ats);
  ElstatPot operator() (const AtomSet& ats, const AtomChargeSet& gen,
			const Coord& cen_intr);
  virtual DielectricEnvironment_lett* make_diel(const AtomSet& ats)
    {return new TwoValueDielectricByAtoms(ats, epsin);}
  virtual ElectrolyteEnvironment_lett* make_ely(const AtomSet& ats)
    {return new ElectrolyteByAtoms(ats);}
protected:
  float epsin;
private:
  AtomSet atlist;
  DielectricEnvironment eps;
  ElectrolyteEnvironment ely;
  FinDiffMethod fdm;
  bool has_eps_ely;
};


bool FD2DielEMaker::atsdiffer(const AtomSet& ats)
{
  if (ats.size() != atlist.size()) return true;
  AtomSet::const_iterator i = ats.begin();
  AtomSet::const_iterator j = atlist.begin();
  while(true) {
    if (i == ats.end() && j == atlist.end()) break;
    if (i == ats.end() || j == atlist.end()) return true;
    if (i->first != j->first || i->second != j->second
	|| i->second.coord != j->second.coord
	|| i->second.rad != j->second.rad) return true;
    ++i;
    ++j;
  }
  return false;
}

ElstatPot
FD2DielEMaker::operator() (const AtomSet& ats, const AtomChargeSet& gen,
			   const Coord& cen_intr)
{
  // If called with a different atom set, use it this one time to make
  // electrolyte and dielectric environments, otherwise use the pre-stored
  // ones.
  ChargeDist rho(new AtomChargeSet(gen));
  fdm.resolve(ats.geom_cent(), cen_intr);
  if (atsdiffer(ats)) {
    return ElstatPot(fdm, DielectricEnvironment(make_diel(ats)),
		     rho, ElectrolyteEnvironment(make_ely(ats)));
  }
  else {
    // make the stored default eps and ely members if necessary
    if (! has_eps_ely) {
      eps = make_diel(atlist);
      ely = make_ely(atlist);
      has_eps_ely=true;
    }
    return ElstatPot(fdm, eps, rho, ely);
  }
}


class FDUnifEMaker : public AtomicElstatMaker {
public:
  FDUnifEMaker (float epsilon, const FinDiffMethod& fm)
    : epsin(epsilon), fdm(fm),
      eps(new UniformDielectric(epsilon)),
      ely(new UniformElectrolyte(0.0))
    {}
  ElstatPot operator() (const AtomSet& ats, const AtomChargeSet& gen,
			const Coord& cen_intr)
    {
      fdm.resolve(ats.geom_cent(), cen_intr);
      ChargeDist rho(new AtomChargeSet(gen));
      return ElstatPot(fdm, eps, rho, ely);
    }
private:
  float epsin;
  FinDiffMethod fdm;
  DielectricEnvironment eps;
  ElectrolyteEnvironment ely;
};

class SavEMaker : public FD2DielEMaker {
public:
  SavEMaker(AtomSet& ats, float eps_in, const FinDiffMethod& fm)
    : FD2DielEMaker(ats, eps_in, fm) {}
  virtual ElectrolyteEnvironment_lett* make_ely(const AtomSet& ats)
    {
      return new ElectrolyteByAtoms(ats, PhysCond::get_ionicstr(),
				    SolvAccVol(ats));
    }
};

#include "MEAD/DielMembAtoms.h"

class FDMembEMaker : public FD2DielEMaker {
public:
  FDMembEMaker(AtomSet& ats, float eps_in, const FinDiffMethod& fm,
	     float zl, float zu, Coord c, float holrad_arg)
    : FD2DielEMaker(ats, eps_in, fm),
      zbot(zl), ztop(zu), holecen(c), holerad(holrad_arg)
    {}
  virtual DielectricEnvironment_lett* make_diel(const AtomSet& ats)
    {
      return new TwoValueDielMembAtoms(ats, epsin, zbot, ztop,
				       holecen, holerad);
    }
private:
  float zbot, ztop, holerad;
  Coord holecen;
};


#ifdef USE_MPI_COARSE
#include <stdio.h>
#include <unistd.h>
#include <iomanip.h>
#include <mpi.h>

int myid, numprocs;
#endif

main(int argc, char* argv[])
{

#ifdef USE_MPI_COARSE
  // Standard MPI preliminaries
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif

  string molname;
  float epsin = 0; // Must be specified on the command line.  0 is error.
  int ismemb = 0;
  float ztop=0, zbot=0;
  float holerad = 0, holex = 0, holey = 0;
  int dosite = 0; // If zero do all sites, postive do just the dosite'th site
  bool ely_by_sav = false; // If true, SolvAccVol is electrolyte boundary.
  bool use_genborn = false; // If true, use pairwise generalized Born model
  // Parse the command line
  for (int iarg=1; iarg<argc; ++iarg) {
    if (iarg+1==argc) {  // This is the last one, must be molname
      molname = argv[iarg];
    }
    // The -nopotats inhibits writing (but not using) potats
    else if ((string) argv[iarg] == (string) "-nopotats")
      SiteUtil_no_potats = true;
    else if ((string) argv[iarg] == (string) "-SelfBack") {
      FlexiSite::do_self_of_back = 1;
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
    else if ((string) argv[iarg] == (string) "-site") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -site");
      istringstream ist(argv[++iarg]);
      int i;
      ist >> i;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -site option");
      dosite = i;
    }

    else if ((string) argv[iarg] == (string) "-epssol") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -epssol");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -epssol option");
      PhysCond::set_epsext(f);
    }
    else if ((string) argv[iarg] == (string) "-epsvac") {
      error ("The -epsvac string is not a valid option for ",
	     argv[0], "\n");
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
    else if ((string) argv[iarg] == string("-ElyBySAV")) {
      ely_by_sav = true;
    }
    else if ((string) argv[iarg] == string("-GenBorn")) {
      use_genborn = true;
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
    // process the form "-membz float float"
    else if ((string) argv[iarg] == (string) "-membz") {
      ismemb = 1;
      if (iarg+2>=argc) error ("main: ERROR, not enough args after -membz");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read first value after -membz opt.");
      zbot = f;
      istringstream ist2(argv[++iarg]);
      ist2 >> f;
      if (!(ist2.good() || ist2.eof()))
	error ("main: FAILED trying to read second value after -membz opt.");
      ztop = f;
    }
    // process the form "-membhole float [float float]"
    else if ((string) argv[iarg] == (string) "-membhole") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -membhole");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read first value after -memhole opt.");
      holerad = f;
      if (iarg+2 < argc) { // There _could_ be two more floats...
	istringstream istx(argv[iarg+1]);
	istx >> f;
	if (istx.good() || istx.eof()) { // Yes, expect this float and another.
	  holex = f;
	  iarg+=2;
	  istringstream isty(argv[iarg]);
	  isty >> f;
	  if (!(isty.good() || isty.eof()))
	    error("main: ERROR, only two floats seen after -membhole.\n",
		  "One or three expected.\n");
	  else
	    holey = f;
	}
      }
    }

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

#ifdef USE_MPI_COARSE
  // Here deal with the problem that parallel processes each produce
  // their own output to stdio which is initiall tied to the same file,
  // hence chaotic output.  Therefore slave processes change their stdout
  // to files labelled with process number and sometimes also site number
  char *myofnamebase = 0;
  char *myofname = 0;
  if (numprocs > 1 && myid == 0) {
    // This is the master process.  Others are slaves.
    // Generate the names for the tmp outfiles of slaves and send them.
    cout << "RUNNING IN PARALLEL MODE.  THIS IS THE MASTER SPEAKING.\n"
      << "STANDARD OUTPUT OF SLAVES WILL BE DIVERTED TO TMP FILES\n"
	<< "AND REINCORPORATED INTO MASTERS STDOUT LATER." << endl;
    for(int sp=1; sp<numprocs; ++sp) {
      ostrstream ost;
      ost << molname << "_p" << setw(3) << setfill('0') << sp << ends;
      char *namebase = ost.str(); // so names like molname_proc{001,002,...}
      blab1 << "sending name base " << namebase << " to process " << sp <<endl;
      MPI_Send(namebase, FILENAME_MAX, MPI_CHAR, sp, 0, MPI_COMM_WORLD);
      delete namebase;
    }
  }
  else if(numprocs>1) {
    myofnamebase = new char[FILENAME_MAX];
    MPI_Status stat;
    MPI_Recv(myofnamebase, FILENAME_MAX, MPI_INT, 0, MPI_ANY_TAG,
	     MPI_COMM_WORLD, &stat);
    ostrstream ost;
    ost << myofnamebase << "_prelim" << ends;
    myofname = ost.str();
    freopen(myofname, "w", stdout);
    cout << "Hello, I am process " << myid << endl;
  }
#endif //  USE_MPI_COARSE

  cout << "Starting " << argv[0] << " for molecule named " << molname << "\n";
  cout << "using the following physical conditions:\n";
  PhysCond::print();
  cout << "\nDielectric constant for molecular interior = " << epsin <<endl;
  if ((PhysCond::get_solrad() > PhysCond::get_sterln())
      && PhysCond::get_ionicstr() != 0 && !ely_by_sav) {
    cerr << "WARNING: The ion excusion radius, or \"Stern length\" has been\n"
	 << "set to a value less than the probe radius that determines the\n"
	 << "solvent accessible volume (SAV) which defines the \"interior\".\n"
	 << "This could cause small bits of \"electrolyte\" to appear in\n"
	 << "the \"interior\".  Is this really what you want?\n"
	 << "Perhaps you should use the flag \"-ElyBySAV\" to specifiy that\n"
	 << "the electrolyte boundary is the SAV, that is, the electrolyte\n"
	 << "and dielectric boundaries are the same.\n"
	 << "For now though, " << argv[0] << " will continue with the\n"
	 << "parameters given." << endl;
  }

  if (ismemb) {
    cout << "A membrane is in the region from z = " << zbot << " to " << ztop;
    cout << ",\nwith a hole of radius " << holerad << " at (x,y) = ("
      << holex << "," << holey << ")\n";
    if (PhysCond::get_ionicstr() != 0)
      error ("Ability to have membrane and nonzero ionic strength not implemented\n");
  }
  if (blab3pt == &cout) cout << "Blab level set to " << 3 << endl;
  else if (blab2pt == &cout) cout << "Blab level set to " << 2 << endl;
  else if (blab1pt == &cout) cout << "Blab level set to " << 1 << endl;
  else cout << "No blab level set (so no blabbing)" << endl;

  AtomSet a(molname);
  a.read();

  FinDiffMethod go;
  string ogm_filename = molname;
  ogm_filename += ".ogm";
  go.read(ogm_filename);
  FinDiffMethod gm;
  string mgm_filename = molname;
  mgm_filename += ".mgm";
  gm.read(mgm_filename);
  cout << "For molecule using finite difference method with lattice levels:\n"
    << go << "\nFor model compounds using finite difference method "
      << "with lattice levels:\n" << gm << endl;

#ifdef USE_MPI_COARSE
  {
    SolvAccVol sav(a);
    if( sav.read_top_in_ascii() == 0 ) {
      blab2 << "Unable to read in SolAccVol information"
	    << " from ascii file either." << endl;
      if(myid==0 && numprocs>1) {
	blab2  << "Calling SolvAccVol::anal_calc to calc. SAV from atom data."
	       << endl;
	sav.anal_calc();
	blab2 << "Returned to main from SolvAccVol::anal_calc" << endl;
	blab2 << "Writing SolvAccVol information to file" << endl;
	sav.write_top_in_ascii();
      }
      // Do a broadcast to signal slaves that the SAV file is ready.
      int dummy;
      MPI_Bcast(&dummy, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (myid!=0 && numprocs>1) {
	sleep(2);
	if (0==sav.read_top_in_ascii())
	  cerr << "ERROR process " << myid
	    << " failed to read its SAV" << endl;
      }
    }
  }
#endif


  AtomicElstatMaker *mac_maker=0, *mod_maker=0, *u_maker=0;
  if (use_genborn) {
    if (ely_by_sav) {
      cerr << "WARNING: The \"-ElyBySAV\" flag was specified, but but this\n"
	   << "run is using generalized Born. "
	   << "Ignoring the \"-ElyBySAV\" flag." << endl;
    }
    if (ismemb) {
      ::error("Membranes not implemented within the Generalized Born model");
    }
    mac_maker = new PGBMaker(epsin);
    mod_maker = new PGBMaker(epsin);
    u_maker = new CoulMaker(epsin);
  }
  else if (ely_by_sav) {
    if (PhysCond::get_ionicstr() != 0) {
      // Set mac_maker etc to something appropriate
      mac_maker = new SavEMaker(a, epsin, go);
      AtomSet empty;
      mod_maker = new SavEMaker(empty, epsin, gm);
      u_maker = new FDUnifEMaker(epsin, go);
    }
    else {
      cerr << "WARNING: The \"-ElyBySAV\" flag was specified, but the ionic\n"
	   << "strength is zero.  Ignoring the flag." << endl;
      mac_maker = new FD2DielEMaker(a, epsin, go);
      AtomSet empty;
      mod_maker = new FD2DielEMaker(empty, epsin, gm);
      u_maker = new FDUnifEMaker(epsin, go);
    }
  }
  else if (ismemb) {
    Coord c(holex, holey, 0.0);
    //  using zbot, ztop, holerad and c, set mac_maker,
    //   etc to something appropriate
    mac_maker = new FDMembEMaker(a, epsin, go, zbot, ztop, c, holerad);
    AtomSet empty;
    mod_maker = new FD2DielEMaker(empty, epsin, gm);
    u_maker = new FDUnifEMaker(epsin, go);
  }
  else { // The normal case
    mac_maker = new FD2DielEMaker(a, epsin, go);
    AtomSet empty;
    mod_maker = new FD2DielEMaker(empty, epsin, gm);
    u_maker = new FDUnifEMaker(epsin, go);
  }

  MultiFlexiSiteMol mol(molname, a, mac_maker, mod_maker, u_maker);
  mol.set_up_sites();

#ifdef USE_MPI_COARSE
  if (myid==0 && numprocs>1) {
    if (dosite>0)
      cerr << "WARNING: the \"-site n\" option is being ignored since\n"
	<< "it makes no sense in a sitewise parallel scheme" << endl;
    int sitessent=0;
    int procend = numprocs;
    if (mol.numsites() < numprocs-1) {
      cerr << "WARNING:  The sitewise parallel scheme can use at most\n"
	<< "   nsites+1 = " << (mol.numsites()+1) << " processors of the "
	  << numprocs << " allocated." << endl;
      procend = mol.numsites()+1;
    }
    for (int ip=1; ip<procend; ++ip) {
      int dummy;
      blab1 << "Master process sending site " << (ip-1) << " to process "
	<< ip << endl;
      MPI_Send(&dummy, 1, MPI_INT, ip, ip-1, MPI_COMM_WORLD);
      ++sitessent;
    }
    // Now wait for slaves to reply when finished and send them new sites.
    for (int irecvd=0; irecvd<mol.numsites(); ++irecvd) {
      MPI_Status stat;
      int dummy;
      MPI_Recv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &stat);
      int sender = stat.MPI_SOURCE;
      blab1 << "Master has recieved word that process " << sender
	<< " has finished with site " << stat.MPI_TAG << endl;
      if (sitessent < mol.numsites()) {
	int dummy;
	blab1 << "Master process sending site " << sitessent
	  << " to the newly available process, " << sender << endl;
	MPI_Send(&dummy, 1, MPI_INT, sender, sitessent, MPI_COMM_WORLD);
	++sitessent;
      }
      else {
	blab1 << "Master process telling process " << sender
	  << " to stop waiting for sites to do." << endl;
	int dummy;
	MPI_Send(&dummy, 1, MPI_INT, sender,
		 mol.numsites()+1, MPI_COMM_WORLD);
      }
    }
  }
  else if (numprocs>1) { // myid != master, so this is a slave process
    while(1) { // BUG FIXME! failure of MPI_Recv could get us stuck here.
      int master = 0;
      int dummy;
      MPI_Status stat;
      MPI_Recv(&dummy, 1, MPI_INT, master, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &stat);
      int isite = stat.MPI_TAG;  // The site to do.
      if (isite >= mol.numsites()) break; // means no more sites to do
      // Change the output file name to reflect site number
      // as per Valerie's suggestion.
      ostrstream ost;
      ost << myofnamebase << "_s" << setfill('0') << setw(3) << isite << ends;
      if (myofname) delete myofname;
      myofname = ost.str();
      freopen(myofname, "w", stdout);
      blab1 << "This is process " << myid
	<< " having recieved assingment to site " << isite << endl;
      mol.calc_site(isite);
      blab1 << "This is process " << myid
	<< " signalling the master that site " << isite << " is done" << endl;
      MPI_Send(&dummy, 1, MPI_INT, master, isite, MPI_COMM_WORLD);
    }
  }
  if (myofname) delete myofname;
  if (myofnamebase) delete myofnamebase;

  // At this point a parallel jobs should be done with "hard" calcs for all
  // sites, but it remains to have the master process re-do the calc_site
  // loop to quickly get all data from potat files and consolidate it.
#endif

#ifdef USE_MPI_COARSE
  if(myid==0) {
    // For a non-parallel job, this is the first run through the sites
    // otherwise, see above.
#endif

    if (dosite > 0 && dosite <= mol.numsites()) { // Do just one site
      mol.calc_site(dosite-1);
    }
    else { // Do them all
      for (int isite=0; isite<mol.numsites(); ++isite) {
	mol.calc_site(isite);
      }
      mol.symmetrize_interactions();
    }
    mol.write_titr_info();
    delete mac_maker;
    delete mod_maker;
    delete u_maker;
    big_rigid_delete_all();
#ifdef USE_MPI_COARSE
  }
  MPI_Finalize();
#endif
  return 0;
}

// multiflex.cc ends here
