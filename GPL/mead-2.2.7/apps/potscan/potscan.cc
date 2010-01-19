#include "MEAD/Potat.h"
#include "MEAD/AtomSet.h"
#include "MEAD/AtomChargeSet.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>
using std::istringstream;
#include <fstream>
using std::ifstream;

main(int argc, char* argv[])
{
  cerr << "WARNING: THIS PROGRAM WILL NOT WORK ON POTAT FILES WRITTEN\n"
       << "WITH VERY OLD VERSIONS OF MEAD THAT USED libg++ CONTAINER Atom.Map"
       << endl;
  blab1pt = &cout;
  blab2pt = &cout;
  blab3pt = &cout;

  string molname = argv[1];
  string sitetype = argv[2];
  string resnumstring = argv[3];
  istringstream istres(argv[3]);
  int resnum;
  istres >> resnum;
  if (!istres) cerr << "WARNING: failed to get resnum from argv[3]" << endl;
  AtomSet a(molname);
  a.read();
  AtomSet atset(a);  // A second copy.  Maybe different internal order?
//  AtomSet atset(molname);
//  atset.read();
  string stfilename = sitetype;
  stfilename += ".st";
  ifstream stfile(stfilename.c_str());
  float pKmod; // This won't be used!
  stfile >> pKmod;
  list<Atom> prot_atlist;
  list<Atom> unprot_atlist;
  while (stfile) {
    string resname, atname;
    stfile >> resname >> atname; // Read atom specifier
    if (!stfile) break;
    AtomID key(resnum, atname);
    float pcharge, ucharge;      // Read the charges
    stfile >> pcharge >> ucharge;
    if (atset.contains(key)) {
      Atom at = atset[key];
      at.charge = pcharge;
      prot_atlist.push_back(at);
      at.charge = ucharge;
      unprot_atlist.push_back(at);
    }
    else {
      cerr << "ERROR: The atom, " << key << ", not found." << endl;
    }
  }
  AtomChargeSet prot(prot_atlist);
  AtomChargeSet unprot(unprot_atlist);
  cout << prot.size() << " atom lines read from " << stfilename << endl;

  int len = atset.size();
  float *protint = new float [len];
  float *unprotint = new float [len];
  float *diffint = new float [len];

  InPotat protpotat(atset);
  InPotat unprotpotat(atset);
  string basepotatfilename = molname;
  basepotatfilename += "." + sitetype + "-" + resnumstring;
  string protpotatfilename = basepotatfilename;
  protpotatfilename += ".p.potat";
  string unprotpotatfilename = basepotatfilename;
  unprotpotatfilename += ".u.potat";
  if (prot.has_charges()) {
    if (protpotat.read(protpotatfilename))
      cout << "protpotat.read returns success" << endl;
    else
      cerr << "WARNING: protpotat.read returns failure" << endl;
  }
  if (unprot.has_charges()) {
    if (unprotpotat.read(unprotpotatfilename))
      cout << "unprotpotat.read returns success" << endl;
    else
      cerr << "WARNING: unprotpotat.read returns failure" << endl;
  }

  AtomSet::const_iterator b = atset.begin();
  for (int i=0; i<len; ++i, ++b) {
    if (prot.contains(b->first)) { // Skip self-interactions
      protint[i] = 0;
      unprotint[i] = 0;
      diffint[i] = 0;
    }
    else {
      AtomChargeSet dummy;
      dummy.insert(b->second);
      protint[i] = prot.has_charges() ? protpotat*dummy : 0;
      unprotint[i] = unprot.has_charges() ? unprotpotat*dummy : 0;
      diffint[i] = protint[i] - unprotint[i];
    }
    cout << i << "  " << b->first
      << "  " << protint[i]
	<< "  " << unprotint[i]
	  << "  " << diffint[i] << endl;
  }

}

// potscan.cc ends here
