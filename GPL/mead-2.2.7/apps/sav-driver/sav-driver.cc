/* Put the SolvAccVol class through some hoops.
    Copyright (c) 1993--1995 by Donald Bashford and Tony You.

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

$Id: sav-driver.cc,v 2.9 2004/11/15 05:30:21 bashford Exp $
*/
#include "MEAD/AtomSet.h"
#include "MEAD/SolvAccVol.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/AccTag.h"
#include "MEAD/PhysCond.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>
using std::istringstream;
#include <fstream>
using std::ofstream;
#include <string>
#include <unistd.h>
#include <sys/times.h>



main(int argc, char* argv[])
{
  string molname = "";
  string lattice_filename = "";
  // Default lattice parameters:
  int grid_dim=41;
  float spacing=1.0;
  Coord grid_center_in_space;
  CubeLatSpec cls;
  int center_specified=0;
  enum FileType {ascii, binary};
  FileType sav_filetype = binary;
  enum LattOutStyle {normal, don_back, tony_back};
  LattOutStyle latt_out_style = normal;
  enum TagMethod {calc_cuberep, tag_pts, call_accessible};
  TagMethod tagging_method = calc_cuberep;
  // Parse the command line
  for (int iarg=1; iarg<argc; ++iarg) {
    if (iarg+1==argc) {  // This is the last one, must be molname
      molname = argv[iarg];
    }
    else if ((string) argv[iarg] == (string) "-dimension") {
      if (iarg+1>=argc)
	error ("main: ERROR, not enough args after -dimension");
      istringstream ist(argv[++iarg]);
      ist >> grid_dim;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -dimension option");
    }
    else if ((string) argv[iarg] == (string) "-spacing") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -spacing");
      istringstream ist(argv[++iarg]);
      ist >> spacing;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -spacing option");
    }
    else if ((string) argv[iarg] == (string) "-center") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -center");
      istringstream ist(argv[++iarg]);
      ist >> grid_center_in_space;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -center option");
      center_specified=1;
    }
    else if ((string) argv[iarg] == (string) "-solrad") {
      if (iarg+1>=argc) error ("main: ERROR, not enough args after -solrad");
      istringstream ist(argv[++iarg]);
      float f;
      ist >> f;
      if (!(ist.good() || ist.eof()))
	error ("main: FAILED trying to read value after -solrad option");
      PhysCond::set_solrad(f);
    }
    else if ((string) argv[iarg] == (string) "-OutputLattice") {
      if (iarg+1>=argc)
	error ("main: ERROR, not enough args after -OutputLattice");
      istringstream ist(argv[++iarg]);
      ist >> lattice_filename;
      if (!(ist.good() || ist.eof()))
	error("main: FAILED trying to read value after -OutputLattice option");
    }
    else if ((string) argv[iarg] == (string) "-ascii") {
      sav_filetype = ascii;
    }
    else if ((string) argv[iarg] == (string) "-don_back") {
      latt_out_style = don_back;
    }
    else if ((string) argv[iarg] == (string) "-tony_back") {
      latt_out_style = tony_back;
    }
    else if ((string) argv[iarg] == (string) "-tag_pts") {
      tagging_method = tag_pts;
    }
    else if ((string) argv[iarg] == (string) "-call_accessible") {
      tagging_method = call_accessible;
    }
    // Process -blab flags
    else if ((string) argv[iarg] == (string) "-blab1") {
      blab1pt = &cout;
    }
    else if ((string) argv[iarg] == (string) "-blab2") {
      blab1pt = &cout;
      blab2pt = &cout;
    }
    else if ((string) argv[iarg] == (string) "-blab3") {
      blab1pt = &cout;
      blab2pt = &cout;
      blab3pt = &cout;
    }
    else {
      cerr << "main: ERROR command option, " << argv[iarg]
	<< ", not recognized" << endl;
      error();
    }
  }

  if (molname == (string) "")
    ::error("main: ERROR, molname not given on command line");
  // Done with comand line processing
  // Print summary of input info from the command line:
  cout << "Starting " << argv[0] << " for molecule named " << molname << "\n";



  AtomSet ats(molname);
  ats.read();
  SolvAccVol  dats(ats);
  struct tms before, after;
  long ticks = sysconf(_SC_CLK_TCK);
  if( dats.read_top_in_binary() == 0 ) {
    blab1 << "Unable to read in SolAccVol from binary file" << endl;
    blab1 << "Try to read the information from ascii file, " << endl;
    if( dats.read_top_in_ascii() == 0 ) {
      blab1 << "Unable to read in SolAccVol information"
	<< " from ascii file either.\n"
	  << "Calling SolvAccVol::anal_calc to calculate it from atom data."
	    << endl;
      times(&before);
      dats.anal_calc();
      times(&after);
      blab1 << "Returned to main from SolvAccVol::anal_calc" << endl;
      cout << "User CPU time for SolvAccVol::anal_calc : "
	<< ((float) (after.tms_utime-before.tms_utime))/ticks
	  << " seconds." << endl;
      cout << "System CPU time for SolvAccVol::anal_calc : "
	<< ((float) (after.tms_stime-before.tms_stime))/ticks
	    << " seconds." << endl;

      switch(sav_filetype) {
      case binary:
	blab1 << "Writing SolvAccVol information to binary file" << endl;
	dats.write_top_in_binary();
	break;
      case ascii:
	blab1 << "Writing SolvAccVol information to ascii file" << endl;
        dats.write_top_in_ascii();
	break;
      }
    }
  }


// Calculate a dielectric grid
  if (center_specified) {
    cls = CubeLatSpec(grid_dim, spacing, grid_center_in_space);
    cls.resolve(ats.geom_cent(), Coord()); // Probably a no-op
  }
  else {
    cls = CubeLatSpec(grid_dim, spacing, ON_GEOM_CENT);
    cls.resolve(ats.geom_cent(), Coord());
    grid_center_in_space = cls.get_center();
  }
  cout << "Accessibility calculation for lattice: " << cls << endl;
  int nsq = grid_dim*grid_dim;
  int ncube = nsq*grid_dim;
  AccTag* acc_array = new AccTag[ncube];

  if (tagging_method==calc_cuberep) {
    blab1 << "Main calling SolvAccVol::calc_cuberep" << endl;
    times(&before);
    dats.calc_cuberep(cls, acc_array);
    times(&after);
  }
  else {
    Coord *pt = new Coord[ncube];
    float grlen = (float) (grid_dim - 1);
    float halfgrlen = grlen/2;
    Coord grid_center_in_grid (halfgrlen, halfgrlen, halfgrlen);
    Coord coord;
    for (int i = 0 ; i<grid_dim; ++i) {
      coord.x = (((float) i) - grid_center_in_grid.x)*spacing
	+ grid_center_in_space.x;
      int iterm = i*nsq;
      for (int j = 0; j<grid_dim; ++j) {
	coord.y = (((float) j) - grid_center_in_grid.y)*spacing
	  + grid_center_in_space.y;
	int ijterms = iterm + j*grid_dim;
	for (int k = 0; k<grid_dim; ++k) {
	  coord.z = (((float) k) - grid_center_in_grid.z)*spacing
	    + grid_center_in_space.z;
	  int h = ijterms + k;
	  pt[h] = coord;
	}
      }
    }
    if (tagging_method==tag_pts) {
      blab1 << "Main calling SolvAccVol::tag_points" << endl;
      times(&before);
      dats.tag_points(ncube, pt, acc_array);
      times(&after);
      blab1 << "Returned to main from SolvAccVol::tag_points:" << endl;
    }
    else if (tagging_method==call_accessible) {
      blab1 << "Main doing calls to SolvAccVol::accessible" << endl;
      times(&before);
      for (int h=0; h<ncube; ++h) {
	if (dats.accessible(pt[h]))
	  acc_array[h] = exterior;
	else
	  acc_array[h] = interior;
	blab3 << "h = " << h << endl;
      }
      times(&after);
      blab1 << "Main finished calls to SolvAccVol::accessible" << endl;
    }
  }


  cout << "User CPU time for marking the points : "
    << ((float) (after.tms_utime-before.tms_utime))/ticks
      << " seconds." << endl;
  cout << "System CPU time for marking the points : "
    << ((float) (after.tms_stime-before.tms_stime))/ticks
      << " seconds." << endl;

  int num_in=0;
  for (int h=0; h<ncube; ++h) {
    if (acc_array[h]==interior) ++num_in;
  }
  cout << "Number of inaccessible points = " << num_in << endl;
// Write grid out to a text file.

  if (lattice_filename != (string) "") {
    blab1 << "Writing out accessibility lattice to file, "
      << lattice_filename << endl;
    ofstream dcrfl(lattice_filename.c_str());
    float epsin = 1.0;
    float epsext = 80.0;
    int i=0;
    switch (latt_out_style) {
    case don_back:
      for (i = 0 ; i<grid_dim; ++i) {
	int iterm = i*nsq;
	for (int j = 0; j<grid_dim; ++j) {
	  int ijterms = iterm + j*grid_dim;
	  for (int k = 0; k<grid_dim; ++k) {
	    int h = ijterms + k;
	    if (acc_array[h] == interior)
	      dcrfl << i << " " << j << " " << k << " " << epsin << "\n";
	    else
	      dcrfl << i << " " << j << " " << k << " " << epsext << "\n";
	  }
	}
      }
      break;
    case tony_back:
      for (i = 0 ; i<grid_dim; ++i) {
	int iterm = i*nsq;
	for (int j = 0; j<grid_dim; ++j) {
	  int ijterms = iterm + j*grid_dim;
	  for (int k = 0; k<grid_dim; ++k) {
	    int h = ijterms + k;
	    if (acc_array[h] == interior)
	      dcrfl << i+1 << " " << j+1 << " " << k+1 << " " << epsin << "\n";
	  }
	}
      }
      break;
    case normal:
      {
      float grlen = (float) (grid_dim - 1);
      float halfgrlen = grlen/2;
      Coord grid_center_in_grid (halfgrlen, halfgrlen, halfgrlen);
      Coord coord;
      for (i = 0 ; i<grid_dim; ++i) {
	coord.x = (((float) i) - grid_center_in_grid.x)*spacing
	  + grid_center_in_space.x;
	int iterm = i*nsq;
	for (int j = 0; j<grid_dim; ++j) {
	  coord.y = (((float) j) - grid_center_in_grid.y)*spacing
	    + grid_center_in_space.y;
	  int ijterms = iterm + j*grid_dim;
	  for (int k = 0; k<grid_dim; ++k) {
	    coord.z = (((float) k) - grid_center_in_grid.z)*spacing
	      + grid_center_in_space.z;
	    int h = ijterms + k;
	    if (acc_array[h] == interior)
	      dcrfl << coord << endl;
	  }
	}
      }
      }
    }
  }
  delete [] acc_array;
  blab1 << "Done." << endl;
}

// sav-driver.cc ends here
