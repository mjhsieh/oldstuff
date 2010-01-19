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

/* Stub versions of MPI C routines (single processor) - most do nothing */

#include "mpi.h"

/* function prototypes */

void mpi_copy_integer(int *, int *, int);
void mpi_copy_float(float *, float *, int);
void mpi_copy_double(double *, double *, int);
void mpi_copy_byte(char *, char *, int);

/* MPI Functions */

void MPI_Comm_rank(MPI_Comm comm, int *me)
{
  *me = 0;
}

void MPI_Comm_size(MPI_Comm comm, int *nprocs)
{
  *nprocs = 1;
}

void MPI_Send(void *buf, int count, MPI_Datatype datatype,
	      int dest, int tag, MPI_Comm comm)
{
  printf("MPI Stub WARNING: Should not send message to self\n");
}

void MPI_Recv(void *buf, int count, MPI_Datatype datatype,
	      int source, int tag, MPI_Comm comm, MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not recv message from self\n");
}

void MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
	       int source, int tag, MPI_Comm comm, MPI_Request *request)
{
  printf("MPI Stub WARNING: Should not recv message from self\n");
}

void MPI_Wait(MPI_Request *request, MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not wait on message from self\n");
}

void MPI_Waitany(int count, MPI_Request *request, int *index, 
		 MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not wait on message from self\n");
}

void MPI_Comm_dup(MPI_Comm comm, MPI_Comm *comm_out) { }

void MPI_Comm_free(MPI_Comm *comm) { }

/* copy values from data1 to data2 */

void MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
		   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  if (datatype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,count);
  else if (datatype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf,count);
  else if (datatype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf,count);
  else if (datatype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf,count);
}

/* copy values from data1 to data2 */

void MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		   void *recvbuf, int recvcount, MPI_Datatype recvtype,
		   MPI_Comm comm)
{
  if (sendtype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf,sendcount);
}

/* copy values from data1 to data2 */

void MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		    void *recvbuf, int *recvcounts, int *displs,
		    MPI_Datatype recvtype, MPI_Comm comm)
{
  if (sendtype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf,sendcount);
}

/* copy values from data1 to data2 */

void MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts,
			MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  if (datatype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,*recvcounts);
  else if (datatype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf,*recvcounts);
  else if (datatype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf,*recvcounts);
  else if (datatype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf,*recvcounts);
}


/*
-------------------
Added routines for data copying
-------------------
*/

void mpi_copy_integer(int *data1, int *data2, int n)
{
  int i;
  for (i = 0; i < n; i++) data2[i] = data1[i];
}

void mpi_copy_float(float *data1, float *data2, int n)
{
  int i;
  for (i = 0; i < n; i++) data2[i] = data1[i];
}

void mpi_copy_double(double *data1, double *data2, int n)
{
  int i;
  for (i = 0; i < n; i++) data2[i] = data1[i];
}

void mpi_copy_byte(char *data1, char *data2, int n)
{
  int i;
  for (i = 0; i < n; i++) data2[i] = data1[i];
}
