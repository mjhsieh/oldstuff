/* A titrating site with several conformational states.  -*- C++ -*-

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


$Id: FlexiSite.h,v 1.2 2004/12/06 17:57:43 bashford Exp $
*/
#ifndef _FlexiSite_h
#define _FlexiSite_h 1

#include <string>
#include "SiteInMulti.h"
#include "MEAD/AtomSet.h"

class FlexiSite : public SiteInMulti {
public:
  FlexiSite(string molname,
	    const AtomSet* all_atp,
	    AtomChargeSet* neut_atp,
	    class AtomicElstatMaker* mac_emaker,
	    class AtomicElstatMaker* mod_emaker,
	    class AtomicElstatMaker* u_emaker);
  virtual void calc_elstat (const vector<SiteInMulti*>& sip,
			    vector<SiteSiteInter>* ssi_row);
  virtual ostream& write_summ_entry(ostream&) const ;
  virtual ~FlexiSite();
  static int do_self_of_back;
private:
  FlexiSite(const FlexiSite&);
  FlexiSite& operator=(const FlexiSite&);
  int read_conf_file();
  void write_summ_file() const ;
  void write_summ_file_with_self_of_back() const ;


  int nconfs;
  int do_models;
  AtomicElstatMaker* u_emaker;
  float state1_free_energy;
  float mod_state1_free_energy;
  float state2_free_energy;
  float mod_state2_free_energy;

  // The following are to be arrays of size nconf for per conformation
  // quantities.  "1" and "2" refer to charge states; "mod" to the model
  // compound; "mac" or no suffix to the macromolecule.

  string *confname;
  float *mac_non_elstat;   // Energies from CHARMM or whatever
  float *mod_non_elstat;
  float *self1_arr;
  float *self2_arr;
  float *back1_arr;
  float *back2_arr;
  float *self_of_back;
  float *mod_self1_arr;
  float *mod_self2_arr;
  float *mod_back1_arr;
  float *mod_back2_arr;
  float *mod_self_of_back;
  float *prob_state1_conf; // Probability of finding state1 macmol in conf
  float *prob_state2_conf;
  float *prob_mod_state1_conf;
  float *prob_mod_state2_conf;

};

#endif

// FlexiSite.h ends here
