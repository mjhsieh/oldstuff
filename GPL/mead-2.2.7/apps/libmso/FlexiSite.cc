/* A titrating site with several conformational states.

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

$Id: FlexiSite.cc,v 1.7 2007/02/21 22:02:50 bashford Exp $
*/
#include "FlexiSite.h"
#include "MEAD/AtomSet.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/UniformDielectric.h"
#include "MEAD/PhysCond.h"
#include "SiteUtil.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/SolvAccVol.h"

#include <math.h>
#include <iomanip>
using std::ios;
using std::setw;

int FlexiSite::do_self_of_back = 0;

typedef AtomSet::const_iterator calitr;

FlexiSite::FlexiSite(string molname_arg,
		     const AtomSet* all_atomp,
		     AtomChargeSet* neut_atomp,
		     AtomicElstatMaker* mac_emaker,
		     AtomicElstatMaker* mod_emaker,
		     AtomicElstatMaker* u_emaker_arg)
: SiteInMulti(molname_arg, all_atomp, neut_atomp, mac_emaker, mod_emaker),
u_emaker(u_emaker_arg),
nconfs(0),
do_models(1),          // Now mandatory.  Later this will be an option.
state1_free_energy(0),
mod_state1_free_energy(0),
state2_free_energy(0),
mod_state2_free_energy(0),
confname(0),
mac_non_elstat(0),
mod_non_elstat(0),
self1_arr(0),
self2_arr(0),
back1_arr(0),
back2_arr(0),
self_of_back(0),
mod_self1_arr(0),
mod_self2_arr(0),
mod_back1_arr(0),
mod_back2_arr(0),
mod_self_of_back(0),
prob_state1_conf(0),
prob_state2_conf(0),
prob_mod_state1_conf(0),
prob_mod_state2_conf(0)
{
}

FlexiSite::FlexiSite(const FlexiSite& fs)
 : SiteInMulti(fs)
{
  error("ERROR: Copy construction of a FlexiSite object not implemented!\n");
}

FlexiSite& FlexiSite::operator=(const FlexiSite& fs)
{
// SiteInMulti::operator=(fs);
  error("ERROR: assignment of a FlexiSite object not implemented!\n");
  return *this;
}

FlexiSite::~FlexiSite()
{
  delete [] confname;
  delete [] mac_non_elstat;
  delete [] mod_non_elstat;
  delete [] self1_arr;
  delete [] self2_arr;
  delete [] back1_arr;
  delete [] back2_arr;
  delete [] self_of_back;
  delete [] mod_self1_arr;
  delete [] mod_self2_arr;
  delete [] mod_back1_arr;
  delete [] mod_back2_arr;
  delete [] mod_self_of_back;
  delete [] prob_state1_conf;
  delete [] prob_state2_conf;
  delete [] prob_mod_state1_conf;
  delete [] prob_mod_state2_conf;

}


int FlexiSite::read_conf_file()
{
// Read the list of protonated conformer names.
// First pass to count number of conformers.
  blab2 << "Looking for a conf file" << endl;
  string conf_fname = molname + "." + sitename + ".confs";
  ifstream fconf1(conf_fname.c_str(), ios::in);
  if (!fconf1) {
    blab1 << "FlexiSite::read_conf_file: Could not open .confs file, \n"
      << conf_fname << endl;
    return 0;
  }
  while (1) {
    string n;
    float e;
    fconf1 >> n >> e;
    if (do_models) fconf1 >> e;
    if (!fconf1.good()) break;
    ++nconfs;
  }
  blab1 << nconfs << " conformer entries read from " << conf_fname
    << " file." << endl;
  // rewind the file
  fconf1.clear();
  fconf1.seekg(0, ios::beg); 
  if (!fconf1) ::error("FlexiSite::read_conf_file: rewind of conf file failed");
// Allocate storage for per-conformer terms and make main pass through confs
  confname = new string [nconfs];
  mac_non_elstat = new float [nconfs];

  self1_arr = new float [nconfs];
  back1_arr = new float [nconfs];
  prob_state1_conf = new float[nconfs];

  self2_arr = new float [nconfs];
  back2_arr = new float [nconfs];
  prob_state2_conf = new float[nconfs];
  if (do_self_of_back)
    self_of_back = new float [nconfs];

  if (do_models) {
    mod_non_elstat = new float [nconfs];
    mod_self1_arr = new float [nconfs];
    mod_back1_arr = new float [nconfs];
    prob_mod_state1_conf = new float[nconfs];
    mod_self2_arr = new float [nconfs];
    mod_back2_arr = new float [nconfs];
    prob_mod_state2_conf = new float[nconfs];
    if (do_self_of_back)
      mod_self_of_back = new float [nconfs];
  }
  int iconf;
  for (iconf = 0; iconf<nconfs; ++iconf) {
    mac_non_elstat[iconf] = 0.0;
    self1_arr[iconf] = 0.0;
    back1_arr[iconf] = 0.0;
    self2_arr[iconf] = 0.0;
    back2_arr[iconf] = 0.0;
    if (do_self_of_back)
      self_of_back[iconf] = 0;
    if (do_models) {
      mod_non_elstat[iconf] = 0.0;
      mod_self1_arr[iconf] = 0.0;
      mod_back1_arr[iconf] = 0.0;
      mod_self2_arr[iconf] = 0.0;
      mod_back2_arr[iconf] = 0.0;
      if (do_self_of_back)
	mod_self_of_back[iconf] = 0;
    }
  }
  for (iconf=0; iconf < nconfs; ++iconf) {
    fconf1 >> confname[iconf] >> mac_non_elstat[iconf];
    mac_non_elstat[iconf] /= PhysCond::get_econv();
    if (do_models) {
      fconf1 >> mod_non_elstat[iconf];
      // FIXME!  Here the energies in the confs file are assumed to
      // be in kcal/mole, but it would be more consistent to have them
      // in pK units!
      mod_non_elstat[iconf] /= PhysCond::get_econv();
      if (!fconf1) ::error("FlexiSite::read_conf_file: ",
			   "attempt to read from conf file "
			   "fails or gives unexpected end-of-file.");
    }
  }
  return 1;
}

void set_coords_of_a_to_b(AtomSet& a, const AtomSet& b)
{
  for (AtomSet::iterator i = a.begin(); i!=a.end(); ++i) {
    AtomSet::const_iterator p = b.find(i->first);
    if (p == b.end()) {
      cerr << "ERROR set_flexibly to: Guide atom set is missing the atom, "
	<<  i->first << endl;
      ::error();
    }
    else {
      i->second.coord = p->second.coord;
    }
  }
}


void FlexiSite::calc_elstat (const vector<SiteInMulti*>& sip,
			     vector<SiteSiteInter>* ssi_row)
{
  if (read_conf_file()) {
    blab1 << "FlexiSite::calc_elstat: doing electrostatic calculations for\n"
      << "multiple conformers of site, " << sitename << "." << endl;
  }
  else {
    blab1 << "FlexiSite::calc_elstat:  Doing electrostatic calculations for\n"
      << "base conformer only for site, " << sitename << "." << endl;
    SiteInMulti::calc_elstat(sip, ssi_row);
    return;
  }

  int nsites = sip.size();
  float *ave_int1 = 0;
  float *ave_int2 = 0;
  if (nsites) {
    ave_int1 = new float[nsites];
    ave_int2 = new float[nsites];
  }

  vector< vector<SiteSiteInter> > intr1_arr(nconfs);
  vector< vector<SiteSiteInter> > intr2_arr(nconfs);
  for (int i=0; i<nconfs; ++i) {
    intr1_arr[i].reserve(nsites);
    intr2_arr[i].reserve(nsites);
  }

  SiteSiteInter *ssnil = 0;

  AtomChargeSet back_chrg(*ref_atp);
  for (calitr b = refstatep->begin(); b!=refstatep->end(); ++b) {
    back_chrg[b->first].charge = 0;
  }
  for (int iconf=0; iconf < nconfs; ++iconf) {
    blab2 << "Calculations for conformer, " << confname[iconf] << endl;
    string molconfname = molname + "." + sitename + "." + confname[iconf];
    AtomChargeSet atlist(molconfname);
    atlist.read();
    Coord cen_intr = atlist[firstkey].coord;

// FIXME!!  Unlike the newest version of SiteInMulti, there is nothing
// here to handle the case that protonated atoms might be different
// in either identity or coordinates from unprotonated ones.
// This should be OK in the case of Tony's lysozyme calculations
// where prot and deprot are the same.

    // Set up and do the macromolecule calculation.
    set_coords_of_a_to_b(back_chrg, atlist);
    set_coords_of_a_to_b(charge_state1, atlist);

    if (do_self_of_back) {

      blab2 << "Doing the self energy of the macromolecule background charges"
	<< "\n for conformer " << confname[iconf]
	  << ".\nFirst the calculation in solvent..." << endl;
      string mol_site_state_name = molconfname + ".b";
      Potat* back_pot = elstat_mol_charge_state(back_chrg, mol_site_state_name,
						atlist, mac_emaker, cen_intr);
      float backself = (*back_pot) * back_chrg;
      delete back_pot;
      mol_site_state_name = molconfname + ".vb";
      Potat* unif_back_pot
	= elstat_mol_charge_state(back_chrg, mol_site_state_name,
			      atlist, u_emaker, cen_intr);
      float unifback = (*unif_back_pot) * back_chrg;
      delete unif_back_pot;

      self_of_back[iconf] = (backself - unifback) / 2.0;
      blab2 << "For state 1 of site " << sitename << " conformer "
	<< confname[iconf] << " of the protein,\n(backself - unifback)/2 = "
	  << "(" << backself << " - " << unifback << ")/2 = "
	    << self_of_back[iconf] << endl;

    }
    blab2 << "Doing the macromolecule calculation for state 1 of conformer "
      << confname[iconf] << endl;
    string mol_site_state_name = molconfname + ".p";
    Potat* mac1_pot =
      elstat_mol_charge_state(charge_state1, mol_site_state_name,
			      atlist, mac_emaker, cen_intr);
    float macself = (*mac1_pot) * charge_state1;
    back1_arr[iconf] = (*mac1_pot) * back_chrg;
    for (int i = 0; i<sip.size(); ++i) {
      intr1_arr[iconf].push_back((*mac1_pot)*sip[i]->get_state1()
				 - (*mac1_pot)*sip[i]->get_state2());
    }
    delete mac1_pot;
    // Set up and do the uniform dielectric calculation.
    // FIXME!  A faster way would be to use only the finest grid.
    blab2 << "Doing uniform dielectric calculation for state 1 of conformer "
      << confname[iconf] << endl;
    mol_site_state_name = molconfname + ".vp";

    Potat* unif1_pot =
      elstat_mol_charge_state(charge_state1, mol_site_state_name,
			      charge_state1, u_emaker, cen_intr);
    float unifself1 = (*unif1_pot) * charge_state1;
    self1_arr[iconf] = (macself - unifself1) / 2.0;
    delete unif1_pot;
    blab2 << "For state 1 of site " << sitename << " conformer "
      << confname[iconf] << " of the protein,\n(macself - unifself)/2 = "
	<< "(" << macself << " - " << unifself1 << ")/2 = "
	  << self1_arr[iconf] << endl;

    blab2 << "Doing the macromolecule calculation for state 2 of conformer "
      << confname[iconf] << endl;
    set_coords_of_a_to_b(charge_state2, atlist);
    mol_site_state_name = molconfname + ".u";
    macself = 0;
    Potat* mac2_pot =
      elstat_mol_charge_state(charge_state2, mol_site_state_name,
			      atlist, mac_emaker, cen_intr);
    macself = (*mac2_pot) * charge_state2;
    back2_arr[iconf] = (*mac2_pot) * back_chrg;
    for (int i = 0; i<sip.size(); ++i) {
      intr2_arr[iconf].push_back((*mac2_pot)*sip[i]->get_state1()
				 - (*mac2_pot)*sip[i]->get_state2());
    }
    delete mac2_pot;

    // Set up and do the uniform dielectric calculation.
    // FIXME!  A faster way would be to use only the finest grid.
    blab2 << "Doing uniform dielectric calculation for state 2 of conformer "
      << confname[iconf] << endl;
    mol_site_state_name = molconfname + ".vu";
    Potat* unif2_pot =
      elstat_mol_charge_state(charge_state2, mol_site_state_name,
			    charge_state2, u_emaker, cen_intr);
    float unifself2 = (*unif2_pot) * charge_state2;
    delete unif2_pot;
    self2_arr[iconf] = (macself - unifself2) / 2.0;

    blab2 << "For state 2 of site " << sitename << " conformer "
      << confname[iconf] << " of the protein,\n(macself - unifself)/2 = "
	<< "(" << macself << " - " << unifself2 << ")/2 = "
	  << self2_arr[iconf] << endl;

    // Do the model compound calculation
    // FIXME! This is not needed if we have predetermined elstat energy
    if (do_models) {
      blab2 << "Doing model compound calculation for state 1 of conformer "
	<< confname[iconf] << endl;
      AtomChargeSet model_compound = make_model_compound(atlist);
      AtomChargeSet model_back_chrg = make_model_compound(back_chrg);
      if (do_self_of_back) {
	blab2 << "Do self energy of the model compound background charges"
	  << "\n for conformer " << confname[iconf]
	    << ".\nFirst the calculation in solvent ... " << endl;
	string mol_site_state_name = molconfname + ".bm";

	cout << "DOES model_back_chrg HAVE CHARGES?.....   ";
	if (model_back_chrg.has_charges())
	  cout << "YES!" << endl;
	else
	  cout << "NO!" << endl;

	Potat* back_pot =
	  elstat_mol_charge_state(model_back_chrg, mol_site_state_name,
				model_compound, mod_emaker, cen_intr);
        float backself = (*back_pot) * model_back_chrg;
	delete back_pot;
	mol_site_state_name = molconfname + ".vbm";
	Potat* unif_back_pot =
	  elstat_mol_charge_state(model_back_chrg, mol_site_state_name,
				  model_compound, u_emaker, cen_intr);
        float unifback = (*unif_back_pot) * model_back_chrg;
	delete unif_back_pot;
	mod_self_of_back[iconf] = (backself - unifback) / 2.0;
	blab2 << "For site " << sitename << " conformer " << confname[iconf]
	  << " of the model compound,\n(backself - unifback)/2 = "
	    << "(" << backself << " - " << unifback << ")/2 = "
	      << mod_self_of_back[iconf] << endl;

      }

      mol_site_state_name = molconfname + ".mp";
      Potat* mod1_pot =
	elstat_mol_charge_state(charge_state1, mol_site_state_name,
				model_compound, mod_emaker, cen_intr);
      float modself = (*mod1_pot) * charge_state1;
      mod_back1_arr[iconf] = (*mod1_pot) * model_back_chrg;
      mod_self1_arr[iconf] = (modself - unifself1) / 2.0;
      delete mod1_pot;
      blab2 << "For state 1 of site " << sitename << " conformer "
	    <<confname[iconf]
	    << " of the model compound,\n(modself - unifself1)/2 = "
	    << "(" << modself << " - " << unifself1 << ")/2 = "
	    << mod_self1_arr[iconf] << endl;


      blab2 << "Doing model compound calculation for state 2 of conformer "
	    << confname[iconf] << endl;
      mol_site_state_name = molconfname + ".mu";
      Potat* mod2_pot =
	elstat_mol_charge_state(charge_state2, mol_site_state_name,
				model_compound, mod_emaker, cen_intr);
      modself = (*mod2_pot) * charge_state2;
      mod_back2_arr[iconf] = (*mod2_pot) * model_back_chrg;
      mod_self2_arr[iconf] = (modself - unifself2) / 2.0;
      delete mod2_pot;

      blab2 << "For state 2 of site " << sitename << " conformer "
	    <<confname[iconf]
	    << " of the model compound,\n(modself - unifself2)/2 = "
	    << "(" << modself << " - " << unifself2 << ")/2 = "
	    << mod_self2_arr[iconf] << endl;
    }
  }

  blab2 << "Finished with loop through conformations.  Now average them."
    << endl;
  float kT = PhysCond::get_kBolt() * PhysCond::get_T();

  float min_energy1 = mac_non_elstat[0] + self1_arr[0] + back1_arr[0];
  if (do_self_of_back)
    min_energy1 += self_of_back[0];
  int i;
  for (i=0; i < nconfs; ++i) {
    float energy = mac_non_elstat[i] + self1_arr[i] + back1_arr[i];
    if (do_self_of_back)
      energy += self_of_back[i];
    min_energy1 = energy < min_energy1 ? energy : min_energy1;
  }

  float canon1_sum = 0;
  for (i=0; i < nconfs; ++i) {
    float energy = mac_non_elstat[i] + self1_arr[i] + back1_arr[i];
    if (do_self_of_back)
      energy += self_of_back[i];
    energy -= min_energy1;
    canon1_sum += (prob_state1_conf[i] = exp(-energy/kT));
  }
  state1_free_energy = -kT * log(canon1_sum) + min_energy1;
  blab1 << "The free energy of charge state 1 in the macromolecule is "
    << state1_free_energy*PhysCond::get_econv() << " kcal/mol\n"
      << "from the following states (energies in kcal/mol):\n"
        << "   non-elstat    self    back   ";
  if (do_self_of_back)
    blab1 << "self of back    ";
  blab1 << "total   gibbs term    probability" << endl;
  for (i=0; i < nconfs; ++i) {
    float energy = mac_non_elstat[i] + self1_arr[i] + back1_arr[i];
    if (do_self_of_back)
      energy += self_of_back[i];
    blab1 << mac_non_elstat[i]*PhysCond::get_econv() << "  "
      << self1_arr[i]*PhysCond::get_econv() << "  "
        << back1_arr[i]*PhysCond::get_econv() << "  ";
    if (do_self_of_back)
      blab1 << self_of_back[i]*PhysCond::get_econv() << "  ";
    blab1 << energy*PhysCond::get_econv() << "  "
      << prob_state1_conf[i] << "*exp(" << min_energy1/kT << ")    ";
    prob_state1_conf[i] /= canon1_sum;
    blab1 << prob_state1_conf[i] << endl;
  }
  blab1 << "Average interaction of charge state 1 with other sites"
    << " (kcal/mole):" << endl;
  int isite;
  for(isite=0; isite<nsites; ++isite) {
    ave_int1[isite] = 0;
    for (int j=0; j < nconfs; ++j)
      ave_int1[isite] += intr1_arr[j][isite] * prob_state1_conf[j];
    blab1 << isite << "   " << ave_int1[isite]*PhysCond::get_econv() << endl;
  }


  float min_energy2 = mac_non_elstat[0] + self2_arr[0] + back2_arr[0];
  if (do_self_of_back)
    min_energy2 += self_of_back[0];
  for (i=0; i < nconfs; ++i) {
    float energy = mac_non_elstat[i] + self2_arr[i] + back2_arr[i];
    if (do_self_of_back)
      energy += self_of_back[i];
    min_energy2 = energy < min_energy2 ? energy : min_energy2;
  }

  float canon2_sum = 0;
  for (i=0; i < nconfs; ++i) {
    float energy = mac_non_elstat[i] + self2_arr[i] + back2_arr[i];
    energy -= min_energy2;
    if (do_self_of_back)
      energy += self_of_back[i];
    canon2_sum += (prob_state2_conf[i] = exp(-energy/kT));
  }
  state2_free_energy = -kT * log(canon2_sum) + min_energy2;
  blab1 << "The free energy of charge state 2 in the macromolecule is "
    << state2_free_energy*PhysCond::get_econv() << " kcal/mol\n"
      << "from the following states (energies in kcal/mol):\n"
        << "   non-elstat    self    back   ";
  if (do_self_of_back)
    blab1 << "self of back    ";
  blab2 << "total   gibbs term    probability" << endl;
  for (i=0; i < nconfs; ++i) {
    float energy = mac_non_elstat[i] + self2_arr[i] + back2_arr[i];
    if (do_self_of_back)
      energy += self_of_back[i];
    blab1 << mac_non_elstat[i]*PhysCond::get_econv() << "  "
      << self2_arr[i]*PhysCond::get_econv() << "  "
        << back2_arr[i]*PhysCond::get_econv() << "  ";
    if (do_self_of_back)
      blab1 << self_of_back[i]*PhysCond::get_econv() << "  ";
    blab1 << energy*PhysCond::get_econv() << "  "
      << prob_state2_conf[i] << "*exp(" << min_energy2/kT << ")    ";
    prob_state2_conf[i] /= canon2_sum;
    blab1 << prob_state2_conf[i] << endl;
  }
  blab1 << "Average interaction of charge state 2 with other sites"
    << " (kcal/mole):" << endl;
  for(isite=0; isite<nsites; ++isite) {
    ave_int2[isite] = 0;
    for (int j=0; j < nconfs; ++j)
      ave_int2[isite] += intr2_arr[j][isite] * prob_state2_conf[j];
    blab1 << isite << "   " << ave_int2[isite]*PhysCond::get_econv() << endl;
  }


  // Now do averages for the model compound
  if(do_models) {
    float min_energy1 = mod_non_elstat[0] + mod_self1_arr[0]
      + mod_back1_arr[0];
    if (do_self_of_back)
      min_energy1 += mod_self_of_back[0];
    int i;
    for (i=0; i < nconfs; ++i) {
      float energy = mod_non_elstat[i] + mod_self1_arr[i] + mod_back1_arr[i];
      if (do_self_of_back)
	energy += mod_self_of_back[i];
      min_energy1 = energy < min_energy1 ? energy : min_energy1;
    }

    float canon1_sum = 0;
    for (i=0; i < nconfs; ++i) {
      float energy = mod_non_elstat[i] + mod_self1_arr[i] + mod_back1_arr[i];
      if (do_self_of_back)
	energy += mod_self_of_back[i];
      energy -= min_energy1;
      canon1_sum += (prob_mod_state1_conf[i] = exp(-energy/kT));
    }
    mod_state1_free_energy = -kT * log(canon1_sum) + min_energy1;
    blab1 << "The free energy of charge state 1 in the model compound is "
      << mod_state1_free_energy*PhysCond::get_econv() << " kcal/mol\n"
	<< "from the following states (energies in kcal/mol):\n"
	  << "   non-elstat    self    back   ";
    if (do_self_of_back)
      blab1 << "self of back    ";
    blab1 << "total   gibbs term    probability" << endl;
    for (i=0; i < nconfs; ++i) {
      float energy = mod_non_elstat[i] + mod_self1_arr[i] + mod_back1_arr[i];
      if (do_self_of_back)
	energy += mod_self_of_back[i];
      blab1 << mod_non_elstat[i]*PhysCond::get_econv() << "  "
	<< mod_self1_arr[i]*PhysCond::get_econv() << "  "
	  << mod_back1_arr[i]*PhysCond::get_econv() << "  ";
      if (do_self_of_back)
	blab1 << mod_self_of_back[i]*PhysCond::get_econv() << "  ";
      blab1 << energy*PhysCond::get_econv() << "  "
	<< prob_mod_state1_conf[i] << "*exp(" << min_energy1/kT << ")    ";
      prob_mod_state1_conf[i] /= canon1_sum;
      blab1 << prob_mod_state1_conf[i] << endl;
    }

    float min_energy2 = mod_non_elstat[0] + mod_self2_arr[0]
      + mod_back2_arr[0];
    if (do_self_of_back)
      min_energy2 += mod_self_of_back[0];
    for (i=0; i < nconfs; ++i) {
      float energy = mod_non_elstat[i] + mod_self2_arr[i] + mod_back2_arr[i];
      if (do_self_of_back)
	energy += mod_self_of_back[i];
      min_energy2 = energy < min_energy2 ? energy : min_energy2;
    }

    float canon2_sum = 0;
    for (i=0; i < nconfs; ++i) {
      float energy = mod_non_elstat[i] + mod_self2_arr[i] + mod_back2_arr[i];
      if (do_self_of_back)
	energy += mod_self_of_back[i];
      energy -= min_energy2;
      canon2_sum += (prob_mod_state2_conf[i] = exp(-energy/kT));
    }
    mod_state2_free_energy = -kT * log(canon2_sum) + min_energy2;
    blab1 << "The free energy of charge state 2 in the model compound is "
      << mod_state2_free_energy*PhysCond::get_econv() << " kcal/mol\n"
	<< "from the following states (energies in kcal/mol):\n"
	  << "   non-elstat    self    back   ";
    if (do_self_of_back)
      blab1 << "self of back    ";
    blab1 << "total   gibbs term    probability" << endl;
    for (i=0; i < nconfs; ++i) {
      float energy = mod_non_elstat[i] + mod_self2_arr[i] + mod_back2_arr[i];
      if (do_self_of_back)
	energy += mod_self_of_back[i];
      blab1 << mod_non_elstat[i]*PhysCond::get_econv() << "  "
	<< mod_self2_arr[i]*PhysCond::get_econv() << "  "
	  << mod_back2_arr[i]*PhysCond::get_econv() << "  ";
      if (do_self_of_back)
	blab1 << mod_self_of_back[i]*PhysCond::get_econv() << "  ";
      blab1 << energy*PhysCond::get_econv() << "  "
	<< prob_mod_state2_conf[i] << "*exp(" << min_energy2/kT << ")    ";
      prob_mod_state2_conf[i] /= canon2_sum;
      blab1 << prob_mod_state2_conf[i] << endl;
    }
  }

  pKint = pKmod
    + (   - state1_free_energy +     state2_free_energy
       +mod_state1_free_energy - mod_state2_free_energy)
                                                 /PhysCond::get_ln10kT();

  blab1 << "For site " << sitename << ", pK_mod of " << pKmod
    << ", shifted to intrinsic pK = " << pKint
      << "\nSite-site interactions:" << endl;
  ssi_row->clear();
  ssi_row->reserve(nsites);
  for(isite=0; isite<nsites; ++isite) {
    ssi_row->push_back(ave_int1[isite] - ave_int2[isite]);
    blab1 << isite+1 << "   " << ssi_row->back() << endl;
  }
  calc_done = 1;
  delete [] ave_int1;
  delete [] ave_int2;
  if(do_self_of_back)
    write_summ_file_with_self_of_back();
  else
    write_summ_file();
}


void FlexiSite::write_summ_file() const
{
  if (!calc_done) return;  // Or should this be an error?  FIXME?
  string summ_fname = molname + "." + sitename + ".summ";
  ofstream ost(summ_fname.c_str());
  if (!ost) {
    cerr << "WARNING: FlexiSite::write_summ_file could not open "
      << summ_fname << " for writing." << endl;
    return;
  }
  // For this file output use scientific notation for floats
  ost.setf(ios::scientific|ios::left);
  ost << "Site " << sitename << " has " << nconfs << " conformers." << endl;

  ost << "Charge state 1 (all energies in pK units):" << endl;
// A little header...
  ost << "\
 Conf. name    non-elstat     self        background     probablilty"
//34567890123456789012345678901234567890123456789012345678901234567890123456789
//       1         2         3         4         5         6         7
//           |             |             |             |             |
    << endl;
  int iconf;
  for (iconf = 0; iconf<nconfs; ++iconf) {
    ost << " " << setw(12) << confname[iconf] << " "
      << setw(13) << mac_non_elstat[iconf]/PhysCond::get_ln10kT() << " "
	<< setw(13) << self1_arr[iconf]/PhysCond::get_ln10kT() << " "
	  << setw(13) << back1_arr[iconf]/PhysCond::get_ln10kT() << " "
	    << setw(13) << prob_state1_conf[iconf]
	      << endl;
  }
  ost << "Free energy of charge state 1 = "
    << state1_free_energy/PhysCond::get_ln10kT() << " (pK units)" << endl;


  ost << "Charge state 2 (all energies in pK units):" << endl;
// A little header...
  ost << "\
 Conf. name    non-elstat     self        background     probablilty"
//34567890123456789012345678901234567890123456789012345678901234567890123456789
//       1         2         3         4         5         6         7
//           |             |             |             |             |
    << endl;
  for (iconf = 0; iconf<nconfs; ++iconf) {
    ost << " " << setw(12) << confname[iconf] << " "
      << setw(13) << mac_non_elstat[iconf]/PhysCond::get_ln10kT() << " "
	<< setw(13) << self2_arr[iconf]/PhysCond::get_ln10kT() << " "
	  << setw(13) << back2_arr[iconf]/PhysCond::get_ln10kT() << " "
	    << setw(13) << prob_state2_conf[iconf]
	      << endl;
  }
  ost << "Free energy of charge state 2 = "
    << state2_free_energy/PhysCond::get_ln10kT() << " (pK units)" << endl;


  if (do_models) {
    ost << "Charge state 1 of model compound (energies in pK units):" << endl;
    // A little header...
    ost << "\
 Conf. name    non-elstat     self        background     probablilty"
//34567890123456789012345678901234567890123456789012345678901234567890123456789
//       1         2         3         4         5         6         7
//           |             |             |             |             |
      << endl;
    int iconf;
    for (iconf = 0; iconf<nconfs; ++iconf) {
      ost << " " << setw(12) << confname[iconf] << " "
	<< setw(13) << mod_non_elstat[iconf]/PhysCond::get_ln10kT() << " "
	  << setw(13) << mod_self1_arr[iconf]/PhysCond::get_ln10kT() << " "
	    << setw(13) << mod_back1_arr[iconf]/PhysCond::get_ln10kT() << " "
	      << setw(13) << prob_mod_state1_conf[iconf]
		<< endl;
    }
    ost << "Free energy of model compound charge state 1 = "
      << mod_state1_free_energy/PhysCond::get_ln10kT()
	<< " (pK units)" << endl;


    ost << "Charge state 2 of model compound (energies in pK units):" << endl;
    // A little header...
    ost << "\
 Conf. name    non-elstat     self        background     probablilty"
//34567890123456789012345678901234567890123456789012345678901234567890123456789
//       1         2         3         4         5         6         7
//           |             |             |             |             |
      << endl;
    for (iconf = 0; iconf<nconfs; ++iconf) {
      ost << " " << setw(12) << confname[iconf] << " "
	<< setw(13) << mod_non_elstat[iconf]/PhysCond::get_ln10kT() << " "
	  << setw(13) << mod_self2_arr[iconf]/PhysCond::get_ln10kT() << " "
	    << setw(13) << mod_back2_arr[iconf]/PhysCond::get_ln10kT() << " "
	      << setw(13) << prob_mod_state2_conf[iconf]
		<< endl;
    }
    ost << "Free energy of model compound charge state 2 = "
      << mod_state2_free_energy/PhysCond::get_ln10kT()
	<< " (pK units)" << endl;
  }
  ost << "Site " << sitename << " pKint = " << pKint << endl;
}

void FlexiSite::write_summ_file_with_self_of_back() const
{
  if (!do_self_of_back) {
    cerr << "ERROR: FlexiSite::write_summ_file_with_self_of_back:\n"
      << "called for object for which self of back calculations are not "
	<< "being done." << endl;
    return;
  }
  if (!calc_done) return;  // Or should this be an error?  FIXME?
  string summ_fname = molname + "." + sitename + ".summ";
  ofstream ost(summ_fname.c_str());
  if (!ost) {
    cerr << "WARNING: FlexiSite::write_summ_file could not open "
      << summ_fname << " for writing." << endl;
    return;
  }
  // For this file output use scientific notation for floats
  ost.setf(ios::scientific|ios::left);
  ost << "Site " << sitename << " has " << nconfs << " conformers." << endl;

  ost << "Charge state 1 (all energies in pK units):" << endl;
// A little header...
  ost << "\
 Conf. name    non-elstat     self        background    self-of-back   probablilty"
//345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
//       1         2         3         4         5         6         7         8
//           |             |             |             |             |             |
    << endl;
  int iconf;
  for (iconf = 0; iconf<nconfs; ++iconf) {
    ost << " " << setw(12) << confname[iconf] << " "
      << setw(13) << mac_non_elstat[iconf]/PhysCond::get_ln10kT() << " "
	<< setw(13) << self1_arr[iconf]/PhysCond::get_ln10kT() << " "
	  << setw(13) << back1_arr[iconf]/PhysCond::get_ln10kT() << " "
	    << setw(13) << self_of_back[iconf]/PhysCond::get_ln10kT() << " "
	    << setw(13) << prob_state1_conf[iconf]
	      << endl;
  }
  ost << "Free energy of charge state 1 = "
    << state1_free_energy/PhysCond::get_ln10kT() << " (pK units)" << endl;


  ost << "Charge state 2 (all energies in pK units):" << endl;
// A little header...
  ost << "\
 Conf. name    non-elstat     self        background    self-of-back   probablilty"
//345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
//       1         2         3         4         5         6         7         8
//           |             |             |             |             |             |
    << endl;
  for (iconf = 0; iconf<nconfs; ++iconf) {
    ost << " " << setw(12) << confname[iconf] << " "
      << setw(13) << mac_non_elstat[iconf]/PhysCond::get_ln10kT() << " "
	<< setw(13) << self2_arr[iconf]/PhysCond::get_ln10kT() << " "
	  << setw(13) << back2_arr[iconf]/PhysCond::get_ln10kT() << " "
	    << setw(13) << self_of_back[iconf]/PhysCond::get_ln10kT() << " "
	    << setw(13) << prob_state2_conf[iconf]
	      << endl;
  }
  ost << "Free energy of charge state 2 = "
    << state2_free_energy/PhysCond::get_ln10kT() << " (pK units)" << endl;

  if (do_models) {
    ost << "Charge state 1 of model compound (energies in pK units):" << endl;
    // A little header...
    ost << "\
 Conf. name    non-elstat     self        background    self-of-back   probablilty"
//345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
//       1         2         3         4         5         6         7         8
//           |             |             |             |             |             |
      << endl;
    int iconf;
    for (iconf = 0; iconf<nconfs; ++iconf) {
      ost << " " << setw(12) << confname[iconf] << " "
	<< setw(13) << mod_non_elstat[iconf]/PhysCond::get_ln10kT() << " "
	  << setw(13) << mod_self1_arr[iconf]/PhysCond::get_ln10kT() << " "
	    << setw(13) << mod_back1_arr[iconf]/PhysCond::get_ln10kT() << " "
	      << setw(13) << mod_self_of_back[iconf]/PhysCond::get_ln10kT() << " "
	      << setw(13) << prob_mod_state1_conf[iconf]
		<< endl;
    }
    ost << "Free energy of model compound charge state 1 = "
      << mod_state1_free_energy/PhysCond::get_ln10kT()
	<< " (pK units)" << endl;


    ost << "Charge state 2 of model compound (energies in pK units):" << endl;
    // A little header...
    ost << "\
 Conf. name    non-elstat     self        background    self-of-back   probablilty"
//345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
//       1         2         3         4         5         6         7         8
//           |             |             |             |             |             |
      << endl;
    for (iconf = 0; iconf<nconfs; ++iconf) {
      ost << " " << setw(12) << confname[iconf] << " "
	<< setw(13) << mod_non_elstat[iconf]/PhysCond::get_ln10kT() << " "
	  << setw(13) << mod_self2_arr[iconf]/PhysCond::get_ln10kT() << " "
	    << setw(13) << mod_back2_arr[iconf]/PhysCond::get_ln10kT() << " "
	      << setw(13) << mod_self_of_back[iconf]/PhysCond::get_ln10kT() << " "
	      << setw(13) << prob_mod_state2_conf[iconf]
		<< endl;
    }
    ost << "Free energy of model compound charge state 2 = "
      << mod_state2_free_energy/PhysCond::get_ln10kT()
	<< " (pK units)" << endl;
  }
  ost << "Site " << sitename << " pKint = " << pKint << endl;
}

ostream&
FlexiSite::write_summ_entry(ostream& ost) const
{
  if (nconfs == 0) return SiteInMulti::write_summ_entry(ost);
  else {
    if (!calc_done) return ost;
    ost << "A " << nconfs << " conformer site with pKint = " << pKint << endl;
    return ost;
  }
}

// FlexiSite.cc ends here
