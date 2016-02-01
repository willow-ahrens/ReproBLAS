#ifndef IDXDMPI_H_
#define IDXDMPI_H_

#include <mpi.h>
#include "idxd.h"

MPI_Op idxdMPI_DIDIADD(const int fold);
MPI_Op idxdMPI_ZIZIADD(const int fold);
MPI_Op idxdMPI_SISIADD(const int fold);
MPI_Op idxdMPI_CICIADD(const int fold);

MPI_Op idxdMPI_DIDIADDSQ(const int fold);
MPI_Op idxdMPI_SISIADDSQ(const int fold);

MPI_Datatype idxdMPI_DOUBLE_INDEXED(const int fold);
MPI_Datatype idxdMPI_DOUBLE_COMPLEX_INDEXED(const int fold);
MPI_Datatype idxdMPI_FLOAT_INDEXED(const int fold);
MPI_Datatype idxdMPI_FLOAT_COMPLEX_INDEXED(const int fold);

MPI_Datatype idxdMPI_DOUBLE_INDEXED_SCALED(const int fold);
MPI_Datatype idxdMPI_FLOAT_INDEXED_SCALED(const int fold);

#endif
