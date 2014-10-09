#ifndef __MPI_INDEXED_SUM_FAST__H
#define __MPI_INDEXED_SUM_FAST__H

#include <mpi.h>

void dIAdd_MPI(void *in, void* inout,
	int* len, MPI_Datatype *datatype);
void zIAdd_MPI(void *in, void* inout,
	int* len, MPI_Datatype *datatype);

void dINrm2_MPI(void *in, void* inout,
	int* len, MPI_Datatype *datatype);

int dIMPICreate(int fold, MPI_Datatype* newtype);
int zIMPICreate(int fold, MPI_Datatype* newtype);
int dIMPIScaleCreate(int fold, MPI_Datatype* newtype);

#endif
