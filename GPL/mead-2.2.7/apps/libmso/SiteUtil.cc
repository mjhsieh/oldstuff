/* Utility used for titrating site calculations

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

$Id: SiteUtil.cc,v 1.1 2001/01/03 23:06:32 bashford Exp $
*/
#include "SiteUtil.h"
#include "MEAD/AtomSet.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/Potat.h"

bool SiteUtil_no_potats = false;  // If true, don't make potat files.

Potat* elstat_mol_charge_state(const AtomChargeSet& gen,
			       string mol_site_state_name,
			       const AtomSet& atlist,
			       AtomicElstatMaker* emaker,
			       const Coord& cen_intr)
{
  // See if there are no charges and we have it easy.
  if (!gen.has_charges()) {
    blab2 << "No charges for this state.  Zeroing terms." << endl;
    Potat* p = new Potat(atlist);
    p->zero();
    return p;
  }

  InPotat* potat = new InPotat(atlist);
  string potat_fname = mol_site_state_name + ".potat";
  if (potat->read(potat_fname)) {
    return potat;
  }
  ElstatPot phi = (*emaker)(atlist, gen, cen_intr);
  phi.solve();
  OutPotat* potat_out = new OutPotat(atlist, phi);
  if (!SiteUtil_no_potats) {
    potat_out->write(potat_fname);
  }
  return potat_out;
}

// SiteUtil.cc ends here
