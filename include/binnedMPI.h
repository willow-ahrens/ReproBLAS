/**
 * @file  binnedMPI.h
 * @brief binnedMPI.h defines MPI wrapper functions for binned types and the necessary functions to perform reproducible reductions.
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
 * - db - binned double (#double_binned)
 * - zb - binned complex double (#double_complex_binned)
 * - sb - binned float (#float_binned)
 * - cb - binned complex float (#float_complex_binned)
 *
 * The goal of using binned types is to obtain either more accurate or reproducible summation of floating point numbers. In reproducible summation, floating point numbers are split into several slices along predefined boundaries in the exponent range. The space between two boundaries is called a bin. Indexed types are composed of several accumulators, each accumulating the slices in a particular bin. The accumulators correspond to the largest consecutive nonzero bins seen so far.
 *
 * The parameter @c fold describes how many accumulators are used in the binned types supplied to a subroutine (an binned type with @c k accumulators  is @c k-fold). The default value for this parameter can be set in config.h. If you are unsure of what value to use for @c fold, we recommend 3. Note that the @c fold of binned types must be the same for all binned types that interact with each other. Operations on more than one binned type assume all binned types being operated upon have the same @c fold. Note that the @c fold of an binned type may not be changed once the type has been allocated. A common use case would be to set the value of @c fold as a global macro in your code and supply it to all binned functions that you use.
 *
 */
#ifndef BINNEDMPI_H_
#define BINNEDMPI_H_

#include <mpi.h>
#include "binned.h"

MPI_Op binnedMPI_DBDBADD(const int fold);
MPI_Op binnedMPI_ZBZBADD(const int fold);
MPI_Op binnedMPI_SBSBADD(const int fold);
MPI_Op binnedMPI_CBCBADD(const int fold);

MPI_Op binnedMPI_DBDBADDSQ(const int fold);
MPI_Op binnedMPI_SBSBADDSQ(const int fold);

MPI_Datatype binnedMPI_DOUBLE_BINNED(const int fold);
MPI_Datatype binnedMPI_DOUBLE_COMPLEX_BINNED(const int fold);
MPI_Datatype binnedMPI_FLOAT_BINNED(const int fold);
MPI_Datatype binnedMPI_FLOAT_COMPLEX_BINNED(const int fold);

MPI_Datatype binnedMPI_DOUBLE_BINNED_SCALED(const int fold);
MPI_Datatype binnedMPI_FLOAT_BINNED_SCALED(const int fold);

#endif
