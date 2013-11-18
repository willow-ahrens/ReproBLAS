#ifndef __MPI_SINDEXED__H
#define __MPI_SINDEXED__H

#include <mpi.h>

void sIAdd_MPI(
	void *invec, void* inoutvec,
	int* len, MPI_Datatype *datatype
);
void cIAdd_MPI(
	void *invec, void* inoutvec,
	int* len, MPI_Datatype *datatype
);
void sINrm2_MPI(
	void *invec, void* inoutvec,
	int* len, MPI_Datatype *datatype
);
int sIMPICreate      (int fold, MPI_Datatype* newtype);
int sIMPIScaleCreate (int fold, MPI_Datatype* newtype);

#endif
