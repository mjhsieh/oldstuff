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

$Id: AtomRegex.h,v 1.4 2005/01/21 20:44:04 bashford Exp $
*/
#ifndef _AtomRegex_h
#define _AtomRegex_h 1

#include <string>
using std::string;

/* 
An AtomRegex works like a POSIX extended regex (like in egrep)
with the added feature that `@' is a standin for the chainid given
to the constructor, `&' is a standin for the residue number
given by the resnum argument to the constructor, and `&+n' and `&-n',
where `n' is an integer, are stand-ins for residue numbers resnum+n
and resnum-n, respectively.

My use of this is to have a convention that canditade atoms are
to be presented in the form, `atname:resname:resnum' so for example:

  char *arxstr = "^CA:.*:&-1$:@" // Match the CA of previous residue
  int resnum = 15;             // This is residue 15.
  string chainid = ""
  AtomRegex arx(arxtr, resnum, chainid);
  arx.matches("CA:ALA:14:");    // Yes, returns 1.
  arx.matches("CA:GLU:15:");   // No, returns 0.

As a more complex example, here is the AtomRegex used as the default
in multiflex for getting a residue, resnum and its flanking peptide
groups:

   AtomRegex dflt("^((.*:.*:&:@)|([CO]:.*:&-1:@)|((N|H|HN|CA):.*:&+1:@))$",
                  resnum);

There could be more const-ness of things here, but defects
in various libraries prevent it.
*/

struct arx_patbuf_t;

class AtomRegex {
public:
  AtomRegex(char * pattern, int resnum, string chainid);
  ~AtomRegex();
  int matches(const string&);
private:
  AtomRegex() {}
  AtomRegex(const AtomRegex&) {}
  AtomRegex& operator=(const AtomRegex&) {return *this;}
  arx_patbuf_t* patbuf;  // an indirect way of referring to a regex_t
};

#endif
