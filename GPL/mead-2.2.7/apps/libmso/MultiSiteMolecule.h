/* Molecule with multiple titrating sites.  -*- C++ -*-

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

$Id: MultiSiteMolecule.h,v 1.2 2004/12/06 17:57:47 bashford Exp $
*/
#ifndef _MultiSiteMolecule_h
#define _MultiSiteMolecule_h 1

#include <string>
#include <vector>
#include "MEAD/AtomSet.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "SiteInMulti.h"

class MultiSiteMolecule {
public:
  MultiSiteMolecule (string name,
		     const AtomSet & AtomList,
		     AtomicElstatMaker* mac_maker,
		     AtomicElstatMaker* mod_maker
		     );
  virtual ~MultiSiteMolecule ();
  void set_up_sites ();
  virtual void calc_site (int i);
  int numsites () const {return silist.size();}
  void symmetrize_interactions();
  void write_titr_info();
protected:
  virtual SiteInMulti* newsite();
  string name;
  AtomSet atlist;
  AtomChargeSet neut_ats;
  AtomicElstatMaker* mac_maker, *mod_maker;
  vector<SiteInMulti*> silist;
  vector< vector<SiteSiteInter> > ssi;
};


#endif

// MultiSiteMolecule.h ends here
