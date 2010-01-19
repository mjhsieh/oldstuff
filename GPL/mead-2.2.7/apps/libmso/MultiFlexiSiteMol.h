// This is -*- C++ -*-
/*
$Id: MultiFlexiSiteMol.h,v 1.2 2004/12/06 17:57:45 bashford Exp $
*/
#ifndef _MultiFlexiSiteMol_h
#define _MultiFlexiSiteMol_h 1

#include <string>
#include "MultiSiteMolecule.h"

class MultiFlexiSiteMol : public MultiSiteMolecule {
public: 
  MultiFlexiSiteMol (string name,
		     const AtomSet & AtomList,
		     AtomicElstatMaker* mac_maker_arg,
		     AtomicElstatMaker* mod_maker_arg,
		     AtomicElstatMaker* u_maker_arg);
  virtual ~MultiFlexiSiteMol();
protected:
  virtual SiteInMulti* newsite();
  AtomicElstatMaker* u_maker;
};

#endif
