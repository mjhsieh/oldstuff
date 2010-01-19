/* A titrating site and a multi-site molecule  -*- C++ -*-

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

$Id: SiteInMulti.h,v 1.4 2005/01/21 20:42:23 bashford Exp $
*/
#ifndef _SiteInMulti_h
#define _SiteInMulti_h 1

#include <string>
#include <iostream>
using std::istream;
#include <vector>
#include "MEAD/ChargeDist.h"
#include "MEAD/Coord.h"
#include "MEAD/AtomSet.h"
#include "MEAD/AtomChargeSet.h"

class DielectricEnvironment;
class ElectrolyteEnvironment;
class DielectricEnvironment_lett;
class ElectrolyteEnvironment_lett;
class FinDiffMethod;
class SolvAccVol;

typedef float SiteSiteInter;

class SiteInMulti {
public:
  SiteInMulti(string molname,
	      const AtomSet* all_atp,
	      AtomChargeSet* neut_atp,
	      class AtomicElstatMaker* mac_emaker,
	      class AtomicElstatMaker* mod_emaker);
  virtual ~SiteInMulti();
  istream& read_sitesfile_item(istream &);
  virtual void setup_site_adjusting_refstate();
  virtual void calc_elstat(const vector<SiteInMulti*>& sip,
			   vector<SiteSiteInter>* ssi_row);

// FIXME should restore this  void print();
  string get_info() const {return info_string;}
  string get_sitename() const {return sitename;}
  const AtomChargeSet& get_state1() const {return charge_state1;}
  const AtomChargeSet& get_state2() const {return charge_state2;}
  AtomChargeSet make_model_compound() {return make_model_compound(*all_atp);}
  ostream& write_pkint_entry(ostream&) const ;
  virtual ostream& write_summ_entry(ostream&) const ;
  static ostream& write_summ_header(ostream& ost);
protected:
  // Copy construction and assignment are not allowed!
  SiteInMulti(const SiteInMulti&);
  SiteInMulti& operator=(const SiteInMulti&);
  AtomChargeSet make_model_compound(const AtomChargeSet&);

  // Identification and calculation state data
  string molname;
  string sitename;
  string sitetype;
  int resnum;
  string chainid;
  int sitesfile_read_done;
  int setup_done;
  int calc_done;

  // Inputs pertaining only to local properties of the site
  float pKmod;
  AtomChargeSet charge_state1;
  AtomChargeSet charge_state2;
  AtomChargeSet *refstatep; // Will point to the "neutral" state
  AtomChargeSet allstate; // Union of all atoms of all possible states.
                          // (i.e. contains all protons);

  // Inputs pertaining to the whole macromolecule
  const AtomSet *all_atp; // Like allstate, but for whole macromolecule
                          // including all sites and "background" atoms.
  AtomChargeSet *ref_atp; // Charge state for whole molecule in which
                          // all sites are in their reference state
                          // (i.e. basis for "intrinsic pK")

  // Methods for preparing Dielectric and Electrolyte environments from atoms.
  AtomicElstatMaker *mac_emaker, *mod_emaker;
  AtomID firstkey; // Defines center of interest

  // Outputs
  float pKint;
  float delta_pK_self;
  float delta_pK_back;
  float macself1;
  float macself2;
  float macback1;
  float macback2;
  float modself1;
  float modself2;
  float modback1;
  float modback2;
  string info_string;
};

class AtomicElstatMaker {
public:
  virtual ElstatPot operator() (const AtomSet&, const AtomChargeSet&,
				const Coord& cen_itr) = 0;
};


#endif

// SiteInMulti.h ends here
