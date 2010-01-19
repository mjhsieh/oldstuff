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

c Dummy parameters for MPI F77 stubs

      parameter (mpi_comm_world = 0)

c recv message status

      parameter (mpi_status_size = 3)

      parameter (mpi_source = 1)
      parameter (mpi_tag = 2)
      parameter (mpi_count = 3)

c recv flags

      parameter (mpi_any_source = -1)
      parameter (mpi_any_tag = -1)

c data types and sizes

      parameter (mpi_integer = 1)
      parameter (mpi_real = 2)
      parameter (mpi_double_precision = 3)
      parameter (mpi_logical = 4)
      parameter (mpi_character = 5)

c allreduce operations

      parameter (mpi_sum = 1)
      parameter (mpi_max = 2)
      parameter (mpi_min = 3)

c timer

      double precision mpi_wtime
