/* A regex class for picking atoms in proteins  -*- C++ -*-

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

$Id: AtomRegex.cc,v 1.5 2005/01/21 20:43:53 bashford Exp $

I implement this using C regex library provided elsewhere, preferably
the POSIX regex provided by the local system.  If that is not available,
try GNU regex, which can also implement POSIX extended syntax.  On
some systems, like SunOS 4.1, there is only Berkeley regex whose syntax
requires escaping ), (, and | with \.  That is a loser as far as I'm
concerned.  Get GNU regex and use it instead.

*/


#include "AtomRegex.h"

#include <sys/types.h>
#include <regex.h>

#ifdef REG_EXTENDED

// Assume we are POSIX equiped.  Nothing more needed

#elif defined (RE_SYNTAX_POSIX_EXTENDED)

// It seems we have old GNU regex which can do POSIX syntax.
typedef struct re_pattern_buffer regex_t;

#include <stdlib.h>

#else

TROUBLE!
  If compilation gets here there is a problem with regex.
  Maybe we dont have regex at all on this system, or maybe
  only the old Berkeley style.
  This text intended to stop compilation if reached.

#endif

// This annoying indirection is needed to allow is to decare a
// rexex_t standin in the AtomRegex.h header file without knowing how
// regex_t was really going to be defined by the system's regex.h header.

struct arx_patbuf_t {regex_t* pt;};

#include "MEAD/globals.h"

#include <sstream>
using std::istringstream;
using std::ostringstream;
using std::flush;

AtomRegex::AtomRegex(char * pattern, int resnum, string chainid)
{
  // Translate the AtomRegex into a normal Regex via strstreams.
  // DEC C++ incorrectly has istrstream::istrstream(char *)
  // rather than istrstream::istrstream(const char *) so we
  // cannot make pattern a const argument.  FIXME?
  istringstream pro_ist(pattern);
  ostringstream reg_ost;
  for (char c = pro_ist.get(); c && !pro_ist.eof(); c = pro_ist.get()) {
    if (c == '&') { // The ampersand calls for special translation
      int delta = 0;
      switch (pro_ist.peek()) {
      case '+' :
      case '-' :
	pro_ist >> delta;
	if (!(pro_ist.good() || pro_ist.eof()))
	  error ("AtomRegex constructor FAILED reading delta\n");
	break;
      }
      reg_ost << (resnum + delta);
    }
    else if (c == '@') { // special translation, substitute chainid
      reg_ost << chainid;
    }
    else  // The standard case: just echo characters
      reg_ost.put(c);
  }
  reg_ost << flush;
  string re_pat = reg_ost.str();
  // Now use regex library to compile pattern
#ifdef REG_EXTENDED
  // if REG_EXTENDED is defined, assume we are POSIX and have regcomp
  // and regexp, and regerror
  const int errmesslen = 500;
  char errmess[errmesslen];
  patbuf = new arx_patbuf_t;
  patbuf->pt = new regex_t;
  int err = regcomp(patbuf->pt, re_pat.c_str(), REG_EXTENDED);
  if (err) {
    int msize = regerror(err, patbuf->pt, errmess, errmesslen);
    cerr << "regcomp returns error condition: " << errmess << endl;
    if (msize>errmesslen) cerr << "Additional text lost." << endl;
    ::error();
  }
#elif defined (RE_SYNTAX_POSIX_EXTENDED)
  // This seems to be a GNU regex that can do posix extended syntax but with
  // gnu-style function calls;
  obscure_syntax = RE_SYNTAX_POSIX_EXTENDED;
  patbuf = new arx_patbuf_t;
  patbuf->pt = new regex_t;
  patbuf->pt->translate = 0; // No character translations;
  patbuf->pt->fastmap = 0;   // Don't make a fastmap;
  patbuf->pt->buffer = 0;    // Let regex to its own allocation;
  patbuf->pt->allocated = 0;
  char *errmess =re_compile_pattern(re_pat.c_str(), re_pat.size(), patbuf->pt);
  if (errmess) {
    cerr << "regcomp returns error condition: " << errmess << endl;
    ::error();
    delete errmess;
  }
#else
  TROUBLE! (see above);
#endif
}

AtomRegex::~AtomRegex()
{
  if (patbuf==0) return;
  if (patbuf->pt == 0) {
    delete patbuf;
    return;
  }
#ifdef REG_EXTENDED
  regfree(patbuf->pt);
#elif defined (RE_SYNTAX_POSIX_EXTENDED)
  free (patbuf->pt->buffer);
#else
  TROBLE! (see above);
#endif
  delete patbuf->pt;
  delete patbuf;
}

int AtomRegex::matches(const string& str)
{

#ifdef REG_EXTENDED
  const int errmesslen = 500;
  char errmess[errmesslen];
  int err = regexec(patbuf->pt, str.c_str(), 0, 0, 0);
  // return value of 0 is match, REG_NOMATCH no match, otherwise error.
  if (!((err==REG_NOMATCH) || (err == 0))) {
    int msize = regerror(err, patbuf->pt, errmess, errmesslen);
    cerr << "regexec returns error condition: " << errmess << endl;
    if (msize>errmesslen) cerr << "Additional text lost." << endl;
  }
  return err == 0;
#elif defined (RE_SYNTAX_POSIX_EXTENDED)
// Call the GNU match function;
// str must be char *, rather than const char * for re_match,
// so some awkard allocation and copying is needed.
  char *pstr = new char [str.length()+1];
  strcopy(pstr, str.length());
  pstr[str.length()] = 0;
  int err = re_match (patbuf->pt, pstr, str.length(), 0, 0);
  // return value >=0 indicates match, -1 is no match, -2 is error.
  delete [] pstr;
  if (err == -2) cerr << "WARNING re_match got internal error" << endl;
  else if (err < -2)
    cerr << "WARNING: re_match returned " << err
      << ", whatever that means" << endl;
  return err>=0;
#else
  TROUBLE! (see above);
#endif
}

// AtomRegex.cc ends here
