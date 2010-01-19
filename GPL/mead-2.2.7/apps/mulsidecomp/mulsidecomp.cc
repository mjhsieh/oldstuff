#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <sstream>
using std::istringstream;
#include <fstream>
using std::ifstream;
#include <string>
using std::string;
#include <map>
#include <vector>

#include "MEAD/AtomSet.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/Potat.h"
#include "SiteInMulti.h"
#include "MEAD/SolvAccVol.h"
#include "MEAD/globals.h"
#include "MEAD/PhysCond.h"
#include "MEAD/AtomSet.h"
// Function template for turning atomwise interactions in the map,
// atitner, for the AtomSet, ats, into some kind of groupwise
// interactions.  Grouping must be a class with a one-argument
// constructor taking an Atom argument, and a operator< that can be
// used to order the groups.  On return, result contains summed
// interactions by grouping, so that if g is a Grouping object
// result[g] will be the sum of interactions of all atoms that fall in
// that grouping.  Return value is the grand sum of all the group
// interactions (useful for sanity checks.  See ByResnum below for
// example of a Grouping

using std::multimap;

template<class Grouping>
float sum_by_group(const AtomSet& ats,
		  const map<AtomID, float>& atinter,
		  map<Grouping, float>* result)
{
  multimap<Grouping,float> minter;
  for (map<AtomID, float>::const_iterator i = atinter.begin();
       i != atinter.end(); ++i) {
    minter.insert( pair<const Grouping, float> (Grouping(ats[i->first]), i->second));
  }
  result->clear();
  typedef typename multimap<Grouping,float>::const_iterator g_citer;
  float mainsum = 0.0;
  for (g_citer j = minter.begin(); j != minter.end(); ) {
    float sum = 0;
    g_citer nextj = minter.upper_bound(j->first);
    for (g_citer i = j; i != nextj; ++i)
      sum += i->second;
    result->insert( pair<Grouping,float> (j->first, sum) );
    mainsum += sum;
    j = nextj;
  }
  return mainsum;
}

// Grouping by residue number.  If ByResnum a(atom1), b(atom2), then
// objects a and b are "equal" (that is "a < b" and "b < a" are both
// false) if atom1 and atom2 have the same residue number.

class ByResnum {
private:
  int _resnum;
public:
  ByResnum() : _resnum(0) {}
  ByResnum(const Atom& a) : _resnum(a.resnum) {}
  ByResnum(int resnum) : _resnum(resnum) {}
  int resnum() const {return _resnum;}
  bool operator<(const ByResnum r) const {return _resnum < r._resnum;}
  virtual ostream& print(ostream& ost) const {return ost << _resnum;}
};
inline ostream& operator<<(ostream& ost, const ByResnum& r)
{return r.print(ost);}

// Grouping by backbone and sidechain of each residue:
// First a utility function ...
inline bool is_BB(const string& atname)
{
  return atname == "N"
    || atname == "H"
    || atname == "HN"
    || atname == "CA"
    || atname == "HA"
    || atname == "C"
    || atname == "O"
    || atname == "HN1"
    || atname == "HN2"
    || atname == "HN3"
    || atname == "O1"
    || atname == "O2"
    || atname == "O3"
    ;
}
// ... then the class.
class BBSC :  public ByResnum {
public:
  typedef enum {backbone, sidechain} bbsc_type;
  BBSC(const Atom& a)
    : _bbsc(is_BB(a.atname) ? backbone : sidechain), ByResnum(a) {}
  BBSC(int resnum, bbsc_type type) : ByResnum(resnum), _bbsc(type) {}
  bool operator<(const BBSC& r) const {
    if (ByResnum::operator<(r)) return true;
    else if (r.ByResnum::operator<(*this)) return false;
    else return _bbsc < r._bbsc;
  }
  bbsc_type type() const {return _bbsc;}
  string bbsc_string() const
    {return _bbsc == backbone? "backbone" : "sidechain";}
  ostream& print(ostream& ost) const
    {ByResnum::print(ost); return ost << " " << bbsc_string();}
private:
  bbsc_type _bbsc;
};


// Fill *atinter with each atom's contribution to background term for
// the site, *siteptr, from a one-conformer multiflex pK calculation.
// Backchrg is the full set of background charges and may also include
// the "generating" charges of *sitptr.  Molname is the molname given
// to multiflex.  By default, the results are expressed as free energy
// of protonation in the macromolecule relative to the model compound
// in units of charge^2/length.  conv_factor can be used to change to
// other units.  For example, to get pK units:
//    conv_factor = -1.0/PhysCond::get_ln10kT()
// If any of the necessary potat files can't be found or read, return
// false, otherwise, succeed and return true.

bool get_atomwise_backterms(SiteInMulti* siteptr,
			    const AtomChargeSet& backchrg,
			    const string& molname,
			    map<AtomID,float> * atinter,
			    float conv_factor = 1.0)
{
  // This function works as follows:  Zero atinter.   For the protonated
  // state of the macromole, read the potat file and add the product
  // of the potential at each atom times the charge at each atom to
  // that atom's element in *atinter.  For the deprotonated
  // macromolecule, do the same except subtract.  Create the model
  // compound, read the potat for its protonated state and subtract its
  // atomwise interactions from *atinter.  For the deprotonated model
  // compound, to the same except subtract.  Finish up by zeroing any
  // elements the correspond to generating atoms (i.e. remove self
  // energy terms), and multiplying all values in *atinter by conv_factor.

  typedef map<AtomID, float>::iterator ati_iter;

  const AtomChargeSet& state1_chrg = siteptr->get_state1();
  const AtomChargeSet& state2_chrg = siteptr->get_state2();

  string mol_site_name = molname + "." + siteptr->get_sitename();

  for (ati_iter k = atinter->begin(); k != atinter->end(); ++k)
    k->second = 0.0;

  if (state1_chrg.has_charges()) { // Do protonated macromol
    InPotat ppot(backchrg);
    if (ppot.read(mol_site_name + ".p.potat"))
      for (AtomChargeSet::const_iterator cp = backchrg.begin();
	   cp != backchrg.end(); ++cp)
	(*atinter)[cp->first] = cp->second.charge * ppot[cp->first];
    else return false;
  }

  if (state2_chrg.has_charges()) { // Do deprotonated macromol
    InPotat ppot(backchrg);
    if (ppot.read(mol_site_name + ".u.potat"))
      for (AtomChargeSet::const_iterator cp = backchrg.begin();
	   cp != backchrg.end(); ++cp)
	(*atinter)[cp->first] -= cp->second.charge * ppot[cp->first];
    else return false;
  }

  AtomChargeSet mod_back = siteptr->make_model_compound();

  if (state1_chrg.has_charges()) { // Do protonated mod comp
    InPotat ppot(mod_back);
    if (ppot.read(mol_site_name + ".mp.potat"))
      for (AtomChargeSet::const_iterator cp = mod_back.begin();
	   cp != mod_back.end(); ++cp)
	(*atinter)[cp->first] -= cp->second.charge * ppot[cp->first];
    else return false;
  }

  if (state2_chrg.has_charges()) { // Do protonated mod comp
    InPotat ppot(mod_back);
    if (ppot.read(mol_site_name + ".mu.potat"))
      for (AtomChargeSet::const_iterator cp = mod_back.begin();
	   cp != mod_back.end(); ++cp)
	(*atinter)[cp->first] += cp->second.charge * ppot[cp->first];
    else return false;
  }

  // Zero generating charge's self interactions
  for (AtomChargeSet::const_iterator cp = state1_chrg.begin();
       cp != state1_chrg.end(); ++cp)
    (*atinter)[cp->first] = 0.0;
  for (AtomChargeSet::const_iterator cp = state2_chrg.begin();
       cp != state2_chrg.end(); ++cp)
    (*atinter)[cp->first] = 0.0;
  // Scale to desired untes
  for (ati_iter k = atinter->begin(); k != atinter->end(); ++k)
    k->second *= conv_factor;

  return true;
}

// Workaround:
// The SiteInMulti constructor needs args with these type signatures...
#include "MEAD/ElstatPot.h"
class BogoMaker : public AtomicElstatMaker {
public:
  ElstatPot operator() (const AtomSet&, const AtomChargeSet&,
				const Coord& cen_itr) {return ElstatPot();}
};


main(int argc, char* argv[])
{
  string molname;
  bool onesite = false;
  string site_to_do;
  for (int iarg=1; iarg<argc; ++iarg) {
    if (iarg+1==argc) {  // This is the last one, must be molname
      molname = argv[iarg];
    }
    else if ((string) argv[iarg] == (string) "-site") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -site");
      onesite = true;
      site_to_do = string(argv[++iarg]);
    }
    else if ((string) argv[iarg] == (string) "-T") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -T");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -T option");
      PhysCond::set_T(f);
    }
    else {
      cerr << "main: ERROR command option, " << argv[iarg]
	<< ", not recognized" << endl;
      error();
    }
  }
  AtomSet ats(molname);
  ats.read();
  AtomChargeSet backchrg(ats);

  // Read the sites file, create and set up a vector of pointers to sites.
  ifstream sitesin(string(molname + ".sites").c_str());
  if (!sitesin)
    ::error("ERROR: the sites file not found.\n");
  vector<SiteInMulti*> site;
  site.reserve(100);
  while (sitesin) {
    DielectricEnvironment* bogus_epsp = 0;
    ElectrolyteEnvironment* bogus_elyp = 0;
    AtomicElstatMaker* bogo_maker = new BogoMaker;
    SiteInMulti* sp = new SiteInMulti(molname, &ats, &backchrg,
				      bogo_maker, bogo_maker);
    if (sp->read_sitesfile_item(sitesin))
      site.push_back(sp);
  }

  // Complete the setup of sites, also read in .summ, .pkint, and .g info.
  map<int,int> resnum_to_sitenum; // for translation of residue # to site #
  vector< vector<float> > site_site_int; // Interaction matrix from .g file
  site_site_int.reserve(site.size());
  vector<float> pKback; // Background term from .summ file
  vector<char> cora; // Is site cation or anion (from .pkint file)

  bool got_summ = false;
  bool got_pKint = false;
  bool got_pKback = false;
  char linebuf[500];
  ifstream summfile(string(molname + ".summ").c_str());
  if (summfile) {
    if (summfile.getline(linebuf, 500)) // Pop off header
      got_summ = true;
  }
  ifstream gfile(string(molname + ".g").c_str());
  ifstream pkintfile(string(molname + ".pkint").c_str());
  for (int i = 0; i < site.size(); ++i) {
    site[i]->setup_site_adjusting_refstate();
    int rnum = -500;
    // A hack to get the residue number of this site
    const AtomChargeSet& st1 = site[i]->get_state1();
    if (st1.has_charges()) rnum = st1.begin()->second.resnum;
    else {
      const AtomChargeSet& st2 = site[i]->get_state2();
      if (st2.has_charges()) rnum = st2.begin()->second.resnum;
    }
    resnum_to_sitenum.insert( pair<int,int> (rnum, i) );
    if (summfile && summfile.getline(linebuf, 500)) {
      string s1, s2, s3; // Junk absorbers
      float v;
      istringstream ist(linebuf);
      ist >> s1 >> s2 >> s3 >> v;
      //      if (ist) pKback.push_back(v);
      //      else summfile.setstate(ios::failbit);
      // FIXME?  setstate, as above, is valid C++, but doesn't work on DEC C++
      // so below, use v unconditioanlly
      pKback.push_back(v);
    }
    if (pkintfile && pkintfile.getline(linebuf, 500)) {
      string s1;
      char c;
      istringstream ist(linebuf);
      ist >> s1 >> c;
      if (ist && ( c=='C' || c=='A' )) cora.push_back(c);
      // else pkintfile.setstate(ios::failbit); // FIXME? (see above)
    }
    if (gfile) {
      site_site_int.push_back( vector<float>() );
      site_site_int.back().reserve(site.size());
      for (int j=0; j < site.size(); ++j) {
	if(gfile.getline(linebuf, 500)) {
	  istringstream ist(linebuf);
	  int j1, j2;  // junk absorbers
	  float v;
	  if (ist >> j1 >> j2 >> v) site_site_int.back().push_back(v);
	  //	  else gfile.setstate(ios::failbit); FIXME? (see above)
	}
      }
    }

  }


  // Set up an associative array of atoms with interactions and zero it.
  map<AtomID, float> atinter;
  for (AtomSet::const_iterator i = ats.begin(); i != ats.end(); ++i)
    atinter.insert( pair<AtomID, float> (i->first, 0.0));

  // Main loop for scanning and analyzing .potat file contents.
  for (int i = 0; i < site.size(); ++i) {
    if (onesite && site[i]->get_sitename() != site_to_do) continue;

    if (!get_atomwise_backterms(site[i], backchrg, molname, &atinter,
				-1.0/PhysCond::get_ln10kT())) {
      cerr << "Skipping site " << site[i]->get_sitename()
	   << " because potat files missing or not read properly" << endl;
      continue;
    }
    cout << "SITE: " << site[i]->get_sitename() << endl;

    float this_pKback = 0.0;
    if (pKback.size() > i)
      this_pKback = pKback[i];
    else
      cerr << "WARNING a pKback value from a summ file was not available\n"
	   << "for this site.  The subsequent sum comparison will be bogus."
	   << endl;

    map<ByResnum,float> resinter;
    float res_test_sum = sum_by_group(ats, atinter, &resinter);
    float diff = res_test_sum - this_pKback;
    float sqdiff = diff*diff;
    if (sqdiff > 0.0001) {
      cerr << "WARNING: sum of the by-residue terms differs from pKback: "
	   << res_test_sum << " vs. " << this_pKback << endl;
    }


    map<BBSC,float> bbscinter;
    float bbsc_test_sum = sum_by_group(ats, atinter, &bbscinter);
    diff = bbsc_test_sum - this_pKback;
    sqdiff = diff*diff;
    if (sqdiff > 0.0001) {
      cerr << "WARNING: sum of the by-residue, backbond or sidechain terms\n"
	   << "differs from pKback: "
	   << bbsc_test_sum << " vs. " << this_pKback << endl;
    }

    bool do_inter = (site_site_int.size() > i
		     && site_site_int[i].size() == site.size()
		     && cora.size() == site.size());
    if (!do_inter)
      cerr << "WARNING: site--site interaction data are not available\n"
	   << "will skip printing that column" << endl;
    for (map<ByResnum,float>::const_iterator k = resinter.begin();
	 k != resinter.end(); ++k) {
      int rnum = k->first.resnum();
      cout << rnum << "\t";
      AtomID cak(rnum,"CA");
      string rname = ats.contains(cak) ? ats[cak].resname : "UNK";
      cout << rname << "\t" << k->second << "\t";
      BBSC isb(rnum, BBSC::backbone);
      if (bbscinter.find(isb) != bbscinter.end())
	cout << bbscinter[isb] << '\t';
      else
	cout << "<nobackbone>" << '\t';
      BBSC iss(rnum, BBSC::sidechain);
      if (bbscinter.find(iss) != bbscinter.end())
	cout << bbscinter[iss] << '\t';
      else
	cout << "<nosidechain>" << '\t';
      if (do_inter) {
	map<int,int>::const_iterator mi = resnum_to_sitenum.find(rnum);
	float v = 0;
	if ( mi != resnum_to_sitenum.end() ) {
	  v = site_site_int[i][mi->second];
	  switch(cora[mi->second]) {
	  case 'A':
	    v *= -1;
	    break;
	  case 'C':
	    break;
	  default:
	    cerr << "ERROR: bogus value in cora vector (see source)" << endl;
	  }
	  v /= -PhysCond::get_ln10kT();
	}
	cout << v;
      }
      cout << endl;
    }
  }
}

// mulsidecomp.cc ends here
