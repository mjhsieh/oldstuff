/* A titrating site and a multi-site molecule

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


$Id: SiteInMulti.cc,v 1.5 2005/01/21 20:41:53 bashford Exp $
*/
#include "SiteInMulti.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/PhysCond.h"
#include "SiteUtil.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/SolvAccVol.h"

#include <iostream>
using std::istream;
using std::ios;
#include <sstream>
using std::ostringstream;
using std::istringstream;
using std::flush;
#include <fstream>
using std::ifstream;

#include <list>

typedef AtomSet::iterator alitr;
typedef AtomSet::const_iterator calitr;

SiteInMulti::SiteInMulti(string molname_arg,
			 const AtomSet* all_atp_arg,
			 AtomChargeSet* neut_arg,
			 AtomicElstatMaker* mac_emaker_arg,
			 AtomicElstatMaker* mod_emaker_arg)
: molname(molname_arg),
  sitename(""),
  sitetype(""),
  resnum(0),
  sitesfile_read_done(0),
  setup_done(0),
  calc_done(0),
  pKmod(0),
  charge_state1(),
  charge_state2(),
  refstatep(0),
  allstate(),
  all_atp(all_atp_arg),
  ref_atp(neut_arg),
  mac_emaker(mac_emaker_arg),
  mod_emaker(mod_emaker_arg),
  firstkey(),
  pKint(0),
  delta_pK_self(0),
  delta_pK_back(0),
  macself1(0),
  macself2(0),
  macback1(0),
  macback2(0),
  modself1(0),
  modself2(0),
  modback1(0),
  modback2(0),
  info_string("")
{
}


SiteInMulti::SiteInMulti(const SiteInMulti& sim)
{
  error("ERROR: Copy construction of a SiteInMulti object not implemented!\n");
}

SiteInMulti& SiteInMulti::operator=(const SiteInMulti& sim)
{
  error("ERROR: assignment of a SiteInMulti object not implemented!\n");
  return *this;
}

SiteInMulti::~SiteInMulti() {
}



istream&
SiteInMulti::read_sitesfile_item(istream& ist)
{
  string line;
  getline(ist, line);
  istringstream lnstr(line);
  lnstr >> resnum;
  lnstr >> sitetype;
  if (lnstr.fail()) {
    // Something is wrong, but caller should deal with it.
    // set state of ist so caller can know.
    ist.setstate(lnstr.rdstate());
    return ist;
  }
  lnstr >> chainid;
  ostringstream ost;
  ost << sitetype << "-" << resnum;
  if (chainid.size() > 0) ost << "-" << chainid;
  ost << flush;
  sitename = ost.str();
  sitesfile_read_done = 1;
  return ist;
}

static string                           // puke!
maybe_dirname (string s)
{
  string slash("/");
  int last_slash = s.rfind (slash);
  // if "/" is not found, npos should be returned, but don't trust
  // implementation to do this.   Instead ...
  bool has_slash = last_slash >=0 && last_slash < s.length();
  string result = has_slash ? s.substr (0, 1 + last_slash) : string ("");
  return result;
}

void
SiteInMulti::setup_site_adjusting_refstate()
{
  if (setup_done)
    ::error("ERROR SiteInMulti::setup_site_adjusting_refstate called twice ",
	    "for same object\n");
  if (!sitesfile_read_done)
    ::error("ERROR SiteInMulti::setup_site_adjusting_refstate called\n",
	    "before reading sitesfile info for this object\n");
  // Open the .st file for this site type and read stuff before atoms:
  string tfn = maybe_dirname (molname) + sitetype + ".st";
  ifstream typein(tfn.c_str(), ios::in);
  if (!typein)
    error ("SiteInMulti::setup_site_adjusting_refstate:  Could not open\n",
	   tfn.c_str(),
	   " for reading site type information.\n");

  char c;
  typein >> c;
  if (c == '/') {
    // FIXME:  This needs to be implemented.
    error("SiteInMulti::setup_site_adjusting_refstate:\n",
	  "new .st file format not yet implemented.\n");
  }
  else {
    typein.putback(c);
    typein >> pKmod;
  }

  // Read in atom-by-atom info.
  float state1_tot = 0, state2_tot = 0;
  list<Atom> atlist1, atlist2, atlist_all;
  bool donefirst = false;
  while (1) {
    string resname, atname;
    // NOTE FIXME resname not really used.  Eliminate it in new format.
    typein >> resname >> atname;
    if (!typein) break;
    AtomID key(resnum, atname, chainid);
    if (!donefirst) {
      firstkey = key;
      donefirst = true;
    }
    bool was_in_1 = false;
    typein >> c;
    typein.putback(c);
    if (c == 'n') {
      string s;
      typein >> s;
      if (s != "nil")
	error("SiteInMulti::setup_site_adjusting_refstate: ERROR reading\n",
	      tfn.c_str(),
	      "; expected \"nil\" or floating-point number.\n");
    }
    else {
      float chr;
      typein >> chr;
      Atom at = (*all_atp)[key];
      at.charge = chr;
      atlist1.push_back(at);
      state1_tot += chr;
      atlist_all.push_back(at);
      was_in_1 = true;
    }

    typein >> c;
    typein.putback(c);
    if (c == 'n') {
      string s;
      typein >> s;
      if (s != "nil")
	error("SiteInMulti::setup_site_adjusting_refstate: ERROR reading\n",
	      tfn.c_str(),
	      "; expected \"nil\" or floating-point number.\n");
    }
    else {
      float chr;
      typein >> chr;
      Atom at = (*all_atp)[key];
      at.charge = chr;
      atlist2.push_back(at);
      state2_tot += chr;
      if (!was_in_1) atlist_all.push_back(at);
    }
    if (typein.eof()) break;
    if (!typein.good())
      error("ERROR: SiteInMulti::setup_site_adjusting_refstate: error reading",
	    tfn.c_str(), "\n");
  }
  charge_state1 = AtomChargeSet(atlist1);
  charge_state2 = AtomChargeSet(atlist2);
  allstate = AtomChargeSet(atlist_all);

  // Now find which state is the nearest to neutral for a reference state
  float abs1 = state1_tot >= 0 ? state1_tot : - state1_tot;
  float abs2 = state2_tot >= 0 ? state2_tot : - state2_tot;
  if (abs1 < abs2) { // State 1 is the neutral state.
    refstatep = &charge_state1;
    blab1 << "Site " << sitename << ": State 1 is the reference state\n"
      << "with total charge = " << state1_tot
	<< ".  State 2 has total charge " << state2_tot << endl;
  }
  else {
    refstatep = &charge_state2;
    blab1 << "Site " << sitename << ": State 2 is the reference state\n"
      << "with total charge = " << state2_tot
	<< ".  State 1 has total charge " << state1_tot << endl;
  }
  float diff = state1_tot - state2_tot;
  if (diff > 1.01 || diff < 0.99) {
    cerr << "WARNING: SiteInMulti::setup_site_adjusting_refstate:\n"
      << "If this is a pH titration calculation it is generally expected\n"
      << "that state1 is the protonated state, and state 2 is deprotonated,\n"
      << "and that the charge_state1 - charge_state2 = 1.\n"
      << "But for this site (" << sitename << ") it is " << diff << "."
      << endl;
  }

  // Now reconcile this site's refstate with the overall refstate.
  DeltaAtomSet delta(atlist_all, *refstatep);
  ref_atp->adjust_atoms_absolute(delta);
  setup_done = 1;
}


#include "AtomRegex.h"

AtomChargeSet
SiteInMulti::make_model_compound(const AtomChargeSet& atlist)
{
  // The model compound is generated from atlist by matching each atom
  // against a regular expression.  Those that match become part of the
  // model compound.  The regular expression is generated using an
  // AtomRegex which is like a POSIX extended regex (as in POSIX egrep)
  // except as follows: "&" is replaced by the residue number
  // associated with the current site "&[+-]n" is replaced by the residue
  // number plus or minus n, where n is an integer.  See AtomRegex.h

  // FIXME.  For now, specify the protoreg here.  In future get it
  // from the .st files

  char *pattern = "^((.*:.*:&:@)|([CO]:.*:&-1:@)|((N|H|HN|CA):.*:&+1:@))$";
  AtomRegex model_regex(pattern, resnum, chainid);

  // Now actally make the model compound
  AtomChargeSet model_list;
  blab2 << "Model compound includes:" << endl;
  for (calitr iat = atlist.begin(); iat!=atlist.end(); ++iat) {
    const Atom& test_atom = iat->second;
    ostringstream test_ost;
    test_ost << test_atom.atname << ":" << test_atom.resname
	     << ":" << test_atom.resnum
	     << ":" << test_atom.chainid
	     << flush;
    string ts = test_ost.str();
    if (model_regex.matches(ts.c_str())) {
      model_list.insert(test_atom);
      blab3 << test_atom;
    }
  }
  blab2 << "Total of " << model_list.size() << " atoms" << endl;
  return AtomChargeSet(model_list);
}


void SiteInMulti::calc_elstat(const vector<SiteInMulti*>& sip,
			      vector<SiteSiteInter>* ssi_row)
{
  if (calc_done)
    ::error("ERROR SiteInMulti::calc_elstat called twice ",
	    "for same object\n");
  if (!setup_done)
    ::error("ERROR SiteInMulti::calc_elstat called\n",
	    "before doing setup for this object\n");

  blab1 << "STARTING CALC_ELSTAT FOR SITE: " << sitename << endl;
  string mol_site_name = molname + "." + sitename;
  Coord cen_intr = (*all_atp)[firstkey].coord;

  vector<SiteSiteInter> sstemp1;
  vector<SiteSiteInter> sstemp2;

  blab1 << "Solving potentials for the macromolecule in charge state 1"
    << endl;
  string mol_site_state_name = mol_site_name + ".p";
  DeltaAtomSet delta1(allstate, charge_state1);
  if (!charge_state1.has_charges()) {
    blab2 << "No charges for this state.  Zeroing terms." << endl;
    macself1 = macback1 = 0;
    sstemp1.reserve(sip.size());
    for (int i = 0; i < sip.size(); ++i) sstemp1.push_back(0);
  }
  else {
    Potat* state1_pot = 0;
    if (delta1.coords_differ()) {
      blab2 << "Atoms for this state differ from \"allstate\"; new eps, ely:"
	    << endl;
      AtomSet atlist(*all_atp);
      atlist.adjust_atoms_absolute(delta1);

      state1_pot = elstat_mol_charge_state(charge_state1, mol_site_state_name,
					   atlist, mac_emaker, cen_intr);
    }
    else {
      blab2 << "Atoms for this state same as \"allstate\"; use all_eps, all_ely:"
	    << endl;

      state1_pot = elstat_mol_charge_state(charge_state1, mol_site_state_name,
					   *all_atp, mac_emaker, cen_intr);
    }
    macself1 = (*state1_pot) * charge_state1;
    // macback is interaction with background charges (ref_atp)
    // EXCEPT those of this site (refstatep).
    macback1 = (*state1_pot) * (*ref_atp) - (*state1_pot) * (*refstatep);
    sstemp1.reserve(sip.size());
    for (int i = 0; i<sip.size(); ++i) {
      sstemp1.push_back((*state1_pot)*sip[i]->get_state1()
			- (*state1_pot)*sip[i]->get_state2());
    }
    delete state1_pot;
  }


  blab1 << "Solving potentials for the macromolecule in charge state 2"
    << endl;
  mol_site_state_name = mol_site_name + ".u";
  DeltaAtomSet delta2(allstate, charge_state2);
  if (!charge_state2.has_charges()) {
    macself2 = macback2 = 0;
    sstemp2.reserve(sip.size());
    for (int i = 0; i < sip.size(); ++i) sstemp2.push_back(0);
  }
  else {
    Potat* state2_pot = 0;
    if (delta2.coords_differ()) {
      blab2 << "Atoms for this state differ from \"allstate\"; new eps, ely:"
	    << endl;
      AtomSet atlist(*all_atp);
      atlist.adjust_atoms_absolute(delta2);
      state2_pot = elstat_mol_charge_state(charge_state2, mol_site_state_name,
					   atlist, mac_emaker, cen_intr);
    }
    else {
      blab2 << "Atoms for this state same as \"allstate\"; use all_eps, all_ely:"
	    << endl;
      state2_pot = elstat_mol_charge_state(charge_state2, mol_site_state_name,
					   *all_atp, mac_emaker, cen_intr);
    }
    macself2 = (*state2_pot) * charge_state2;
    macback2 = (*state2_pot) * (*ref_atp) - (*state2_pot) * (*refstatep);
    sstemp2.reserve(sip.size());
    for (int i = 0; i<sip.size(); ++i) {
      sstemp2.push_back((*state2_pot)*sip[i]->get_state1()
			- (*state2_pot)*sip[i]->get_state2());
    }
    delete state2_pot;
  }
  int i;
  ssi_row->clear();
  ssi_row->reserve(sip.size());
  for (i = 0; i<sip.size(); ++i)
    ssi_row->push_back(sstemp1[i] - sstemp2[i]);

  blab1 << "Creating the model compound" << endl;
  // The atoms to be used in calculating dielectric boundaries, etc.
  AtomChargeSet model_compound = make_model_compound(*all_atp);
  // The charges to be used in calculating the background term.
  // Use *ref_atp so that if an of the model compound atom happens
  // to be part of some other site (as might happen at chain termini)
  // the charges already reflect the neutral state.
  AtomChargeSet model_back_chrg = make_model_compound(*ref_atp);
  // Any atoms that are "titrating" in *this* site should be zero in
  // the background set.  This requires adjustment...
  for (alitr b = model_back_chrg.begin(); b!=model_back_chrg.end() ; ++b) {
    if (refstatep->contains(b->first)) (b->second).charge = 0;
  }
  vector<SiteInMulti*> ssempty(size_t(0));
  blab1 << "Solving potentials for the model compound in charge state 1"
    << endl;
  mol_site_state_name = mol_site_name + ".mp";
  if (!charge_state1.has_charges()) {
    modself1 = modback1 = 0;
  }
  else {
    Potat* mod1_pot = 0;
    if (delta1.coords_differ()) {
      blab2 << "Atoms for this state differ from \"allstate\"; new eps, ely:"
	    << endl;
      AtomSet atlist(model_compound);
      atlist.adjust_atoms_absolute(delta1);
      mod1_pot = elstat_mol_charge_state(charge_state1, mol_site_state_name,
					 atlist, mod_emaker, cen_intr);
    }
    else {
      blab2 << "Atoms for this state same as \"allstate\"; use mod_eps, mod_ely:"
	    << endl;
      mod1_pot = elstat_mol_charge_state(charge_state1, mol_site_state_name,
					 model_compound, mod_emaker, cen_intr);
    }
    modself1 = (*mod1_pot) * charge_state1;
    modback1 = (*mod1_pot) * model_back_chrg;
    delete mod1_pot;
  }

  blab1 << "Solving potentials for the model compound in charge state 2"
    << endl;
  mol_site_state_name = mol_site_name + ".mu";
  if (!charge_state2.has_charges()) {
    modself2 = modback2 = 0;
  }
  else {
    Potat* mod2_pot = 0;
    if (delta2.coords_differ()) {
      blab2 << "Atoms for this state differ from \"allstate\"; new eps, ely:"
	    << endl;
      AtomSet atlist(model_compound);
      atlist.adjust_atoms_absolute(delta2);
      mod2_pot = elstat_mol_charge_state(charge_state2, mol_site_state_name,
					 atlist, mod_emaker, cen_intr);
    }
    else {
      blab2 << "Atoms for this state same as \"allstate\"; use mod_eps, mod_ely:"
	    << endl;
      mod2_pot = elstat_mol_charge_state(charge_state2, mol_site_state_name,
					 model_compound, mod_emaker, cen_intr);
    }
    modself2 = (*mod2_pot) * charge_state2;
    modback2 = (*mod2_pot) * model_back_chrg;
    delete mod2_pot;
  }


  delta_pK_self = -(macself1-macself2-modself1+modself2)/2.0
    / PhysCond::get_ln10kT();
  delta_pK_back = -(macback1-macback2-modback1+modback2)
    / PhysCond::get_ln10kT();


  blab1 << "macself1-macself2 = " << macself1-macself2
    << "\nmacback1-deprotpack = " << macback1-macback2
      << "\nmodself1-modself2 = " << modself1-modself2
	<< "\nmodback1-modback2 = "<< modback1-modback2 <<"\n";
  blab1 << "delta self = " << delta_pK_self * PhysCond::get_ln10kT()
    << " or " << delta_pK_self << " pK units\n";
  blab1 << "delta back = " << delta_pK_back * PhysCond::get_ln10kT()
    <<  " or " << delta_pK_back << " pK units\n";
  blab1 << "Site-site interactions:" << endl;
  for (i=0; i<ssi_row->size(); ++i)
    blab1 << i << "   " << (*ssi_row)[i] << endl;

  pKint = pKmod + (delta_pK_self + delta_pK_back);
  calc_done = 1;

}



// FIXME this pkint file format should be changed, especiall the C or A
// thing.

ostream&
SiteInMulti::write_pkint_entry(ostream& ost) const
{
  if (!calc_done) return ost;
  ost.setf(ios::scientific); // What is floatfield?
  ost << pKint << " ";
  if (refstatep == &charge_state1)
    ost << 'A';  // An anionic site
  else if (refstatep == &charge_state2)
    ost << 'C';  // An cationic site;
  else
    cerr << "ERROR SiteInMulti::write_pkint_line: Could not determine "
      << "C or A." << endl;
  ost << " " << sitename << endl;
  return ost;
}

ostream&
SiteInMulti::write_summ_header(ostream& ost)
{
  ost << "\
   site name           pKmod      delta self    delta back      pkint"
//34567890123456789012345678901234567890123456789012345678901234567890123456789
//       1         2         3         4         5         6         7
//           |             |             |             |             |
    << endl;
  return ost;
}

#include <iomanip>
using std::setw;

ostream&
SiteInMulti::write_summ_entry(ostream& ost) const
{
  if (!calc_done) return ost;
  ost << " " << setw(12) << sitename << " " << setw(13) << pKmod << " "
    << setw(13) << delta_pK_self << " " << setw(13) << delta_pK_back
      << " " << setw(13) << pKint << endl;
  return ost;
}

// SiteInMulti.cc ends here
