/*
 *  Created   13/10/25   H.D. Nguyen
 */

#ifndef _REPRODUCIBLE_MPI_FLOAT__H_
#define _REPRODUCIBLE_MPI_FLOAT__H_

#include <mpi.h>
#include "idxd.h"

extern MPI_Datatype MPI_IDOUBLE;
extern MPI_Datatype MPI_IFLOAT;
extern MPI_Datatype MPI_ICOMPLEX;
extern MPI_Datatype MPI_IDOUBLE_COMPLEX;
extern MPI_Datatype MPI_IDOUBLE_SCALE;
extern MPI_Datatype MPI_IFLOAT_SCALE;

extern MPI_Op MPI_RSUM;
extern MPI_Op MPI_RNRM2;

extern int RMPI_Init();
extern void IAdd_MPI(void *in, void* inout,
	int* len, MPI_Datatype *datatype);
extern void INrm2_MPI(void *in, void* inout,
	int* len, MPI_Datatype *datatype);

#endif
