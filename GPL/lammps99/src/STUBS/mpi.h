/*
c LAMMPS 99 - Molecular Dynamics Simulator
c Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/lammps.html
c Steve Plimpton, sjplimp@sandia.gov
c
c Copyright (2003) Sandia Corporation.  Under the terms of Contract
c DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
c certain rights in this software.  This software is distributed under 
c the GNU General Public License.
c
c See the README file in the top-level LAMMPS directory.
*/

/* Dummy defs for MPI C stubs */

#define MPI_Comm int
#define MPI_Request int
#define MPI_Status int
#define MPI_Datatype int
#define MPI_Op int

#define MPI_INT 1
#define MPI_FLOAT 2
#define MPI_DOUBLE 3
#define MPI_DOUBLE_PRECISION 3
#define MPI_BYTE 4

#define MPI_SUM 1
#define MPI_MAX 2
#define MPI_MIN 3
