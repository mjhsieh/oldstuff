// This is -*- C++ -*-
/*
$Id: SiteUtil.h,v 1.1 2001/01/03 23:06:53 bashford Exp $
*/
#ifndef _SiteUtil_h
#define _SiteUtil_h 1

#include <string>
#include "SiteInMulti.h"
#include "MEAD/Potat.h"

extern bool SiteUtil_no_potats;  // If true, don't make potat files.

Potat* elstat_mol_charge_state(const AtomChargeSet& gen,
			       string mol_site_state_name,
			       const AtomSet& atlist,
			       AtomicElstatMaker* emaker,
			       const Coord& cen_intr);

#endif

// SiteUtil.h ends here
