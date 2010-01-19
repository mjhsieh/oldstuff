/*
$Id: MultiFlexiSiteMol.cc,v 1.2 2004/12/06 17:57:44 bashford Exp $
*/
#include "MultiFlexiSiteMol.h"
#include "FlexiSite.h"

MultiFlexiSiteMol::MultiFlexiSiteMol (string name_arg,
				      const AtomSet & AtomList,
				      AtomicElstatMaker* mac_maker_arg,
				      AtomicElstatMaker* mod_maker_arg,
				      AtomicElstatMaker* u_maker_arg)
: MultiSiteMolecule(name_arg, AtomList, mac_maker_arg, mod_maker_arg),
  u_maker(u_maker_arg)
{}


MultiFlexiSiteMol::~MultiFlexiSiteMol () {}


SiteInMulti* MultiFlexiSiteMol::newsite()
{
  return new FlexiSite(name, &atlist, &neut_ats, mac_maker, mod_maker,
		       u_maker);
}
