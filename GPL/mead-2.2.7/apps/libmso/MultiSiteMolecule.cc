/* Molecule with multiple titrating sites.

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

$Id: MultiSiteMolecule.cc,v 1.3 2004/12/06 17:57:46 bashford Exp $
*/
#include "MultiSiteMolecule.h"
#include "SiteInMulti.h"
#include "MEAD/DielByAtoms.h"
#include "MEAD/ElectrolyteByAtoms.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/PhysCond.h"
#include "MEAD/Potat.h"
#include "MEAD/globals.h"
#include "MEAD/SolvAccVol.h"

#include <math.h>
#include <fstream>
using std::ofstream;
using std::ifstream;
using std::ios;

MultiSiteMolecule::MultiSiteMolecule (string nm,
				      const AtomSet & AtomList,
				      AtomicElstatMaker* mac_maker_arg,
				      AtomicElstatMaker* mod_maker_arg)
: name(nm), atlist(AtomList), neut_ats(AtomList),
  mac_maker(mac_maker_arg), mod_maker(mod_maker_arg)
{}

MultiSiteMolecule::~MultiSiteMolecule ()
{
  for (int i=0; i<silist.size(); ++i) {
    delete silist[i];
  }
}

void MultiSiteMolecule::calc_site (int ind)
{
  if (ind < 0 || ind >= silist.size()) {
    cerr << "ERROR: MultiSiteMolecule::calc_site: "
      << ".\nCannot do calculation for site number " << ind
	<< " (numbering from zero.\n"
	  << "Either this MultiSiteMolecule not set up right or index "
	    << "out of range" << endl;
    error();
  }
  silist[ind]->calc_elstat(silist, &ssi[ind]);
  ssi[ind][ind] = 0; // No self-interaction through the ssi matrix allowed
}

void MultiSiteMolecule::symmetrize_interactions()
{
  int nsites = silist.size();
  for (int isite=0; isite<nsites; ++isite)
    if (ssi[isite].size()!=nsites) {
      cerr << "WARNING MultiSiteMolecule::symmetrize_interactions called\n"
	<< "before all interactions are calculated."
	  << "Returning with no symmetrization." << endl;

      return;
    }
  cout << "Site-site interaction symmetry check:" << endl;
  cout.setf(ios::scientific);
  for (int i=0; i<nsites-1; ++i) {
    for (int j=i+1; j<nsites; ++j) {
      float ave = (ssi[i][j] + ssi[j][i]) / 2.0;
      float dev;
      if (ave == 0) {
	cout << "WARNING: ave = 0\n";
	dev = 0;
      }
      else {
	dev = (ssi[i][j] - ssi[j][i])/ave;
	dev = dev > 0 ? dev : -dev;
      }
      cout << i+1 << " " << j+1 << "   "
	<< ssi[i][j] << " " << ssi[j][i] << "   dev = " << dev << endl;
      ssi[i][j] = ssi[j][i] = ave;
    }
  }
  cout.unsetf(ios::scientific);
  cout << "Site site interactions have now been symmetrized." << endl;
}

SiteInMulti* MultiSiteMolecule::newsite()
{
  return new SiteInMulti(name, &atlist, &neut_ats, mac_maker, mod_maker);
}


void MultiSiteMolecule::set_up_sites ()
{
  if (silist.size() != 0)
    error ("ERROR: MultiSiteMolecule::set_up_sites:\n",
	   "Sites appear to be already defined\n");
  string site_filename = name + DOT + "sites";
  ifstream sitesin(site_filename.c_str(), ios::in);
  if (!sitesin)
    error ("ERROR: MultiSiteMolecule::set_up_sites: Unable to open file,\n",
	   site_filename.c_str(),
	   ", for reading list of sites.\n");
  ssi.clear();
  while(true) {
    SiteInMulti* sp = newsite();
    if ( ! sp->read_sitesfile_item(sitesin) ) break;
    silist.push_back(sp);
    ssi.push_back( vector<SiteSiteInter>() );
    silist.back()->setup_site_adjusting_refstate();
    // Make sure there are no duplicates.
    string iname = silist.back()->get_sitename();
    for (int j=0; j<silist.size()-1; ++j)
      if (silist[j]->get_sitename() == iname)
	cerr << "WARNING MultiSiteMolecule::set_up_sites: Duplicate site, "
	  << iname << endl;
  }
  blab2 << "MultiSiteMolecule::set_up_sites: " << silist.size()
    << " sites read from " << site_filename << endl;
}


void MultiSiteMolecule::write_titr_info()
{
  int nsites = silist.size();
  if (nsites<1) return;
  string pkint_filename = name + ".pkint";
  ofstream fpkint (pkint_filename.c_str());
  if (!fpkint)
    cerr << "WARNING: MultiSiteMolecule::write_titr_info: could not open "
      << pkint_filename << " for writing." << endl;
  else {
    for (int i=0; i<nsites; ++i) {
      silist[i]->write_pkint_entry(fpkint);
    }
  }

  string pksumm_filename = name + ".summ";
  ofstream fsumm (pksumm_filename.c_str());
  if (!fsumm)
    cerr << "WARNING: MultiSiteMolecule::write_titr_info: could not open "
      << pksumm_filename << " for writing." << endl;
  else {
    SiteInMulti::write_summ_header(fsumm);
    for (int i=0; i<nsites; ++i) {
      silist[i]->write_summ_entry(fsumm);
    }
  }

  string g_filename = name + DOT + "g";
  ofstream fg (g_filename.c_str());
  if (!fg)
    cerr << "WARNING: MultiSiteMolecule::write_titr_info: could not open "
      << g_filename << " for writing." << endl;
  else {
    fg.setf(ios::scientific);
    for (int i = 0; i<nsites; ++i)
      if (ssi[i].size())
	for (int j=0; j<nsites; ++j)
	  fg << i+1 << " " << j+1 << "   " << ssi[i][j] << endl;
  }
}

// MultiSiteMolecule.cc ends here
