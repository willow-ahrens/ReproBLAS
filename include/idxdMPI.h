/**
 * @file  idxdMPI.h
 * @brief idxdMPI.h defines MPI wrapper functions for indexed types and the necessary functions to perform reproducible reductions.
 *
 * This header is modeled after cblas.h, and as such functions are prefixed with character sets describing the data types they operate upon. For example, the function @c dfoo would perform the function @c foo on @c double possibly returning a @c double.
 *
 * If two character sets are prefixed, the first set of characters describes the output and the second the input type. For example, the function @c dzbar would perform the function @c bar on @c double @c complex and return a @c double.
 *
 * Such character sets are listed as follows:
 * - d - double (@c double)
 * - z - complex double (@c *void)
 * - s - float (@c float)
 * - c - complex float (@c *void)
 * - di - indexed double (#double_indexed)
 * - zi - indexed complex double (#double_complex_indexed)
 * - si - indexed float (#float_indexed)
 * - ci - indexed complex float (#float_complex_indexed)
 *
 * The goal of using indexed types is to obtain either more accurate or reproducible summation of floating point numbers. In reproducible summation, floating point numbers are split into several slices along predefined boundaries in the exponent range. The space between two boundaries is called a bin. Indexed types are composed of several accumulators, each accumulating the slices in a particular bin. The accumulators correspond to the largest consecutive nonzero bins seen so far.
 *
 * The parameter @c fold describes how many accumulators are used in the indexed types supplied to a subroutine (an indexed type with @c k accumulators  is @c k-fold). The default value for this parameter can be set in config.h. If you are unsure of what value to use for @c fold, we recommend 3. Note that the @c fold of indexed types must be the same for all indexed types that interact with each other. Operations on more than one indexed type assume all indexed types being operated upon have the same @c fold. Note that the @c fold of an indexed type may not be changed once the type has been allocated. A common use case would be to set the value of @c fold as a global macro in your code and supply it to all indexed functions that you use.
 *
 */
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
