#include <mpi.h>
#include <stdio.h>

#include <idxd.h>
#include <idxdMPI.h>

#include "../../config.h"

static MPI_Datatype types[idxd_SIMAXFOLD + 1];
static int types_initialized[idxd_SIMAXFOLD + 1]; //initializes to 0

/**
 * @brief  Get an MPI_DATATYPE representing scaled indexed single precision.
 *
 * Creates (if it has not already been created) and returns a datatype handle
 * for an MPI datatype that represents a scaled indexed single precision type.
 *
 * This method may call @c MPI_Type_contiguous() and @c MPI_Type_commit().
 * If there is an error, this method will call @c MPI_Abort().
 *
 * @param fold the fold of the indexed types
 *
 * @author Peter Ahrens
 * @date   18 Jun 2016
 */
MPI_Datatype idxdMPI_FLOAT_INDEXED_SCALED(const int fold){
  int rc;
  if(!types_initialized[fold]){
    rc = MPI_Type_contiguous(idxd_sinum(fold) + 1, MPI_FLOAT, types + fold);
    if(rc != MPI_SUCCESS){
      if (rc == MPI_ERR_TYPE) {
        fprintf(stderr, "[%s.%d] ReproBLAS error: MPI_Type_contiguous error: MPI_ERR_TYPE\n", __FILE__, __LINE__);
      } else if (rc == MPI_ERR_COUNT) {
        fprintf(stderr, "[%s.%d] ReproBLAS error: MPI_Type_contiguous error: MPI_ERR_COUNT\n", __FILE__, __LINE__);
      } else if (rc == MPI_ERR_INTERN) {
        fprintf(stderr, "[%s.%d] ReproBLAS error: MPI_Type_contiguous error: MPI_ERR_INTERN\n", __FILE__, __LINE__);
      } else {
        fprintf(stderr, "[%s.%d] ReproBLAS error: MPI_Type_contiguous error: %d\n", __FILE__, __LINE__, rc);
      }
      MPI_Abort(MPI_COMM_WORLD, rc);
      return 0;
    }
    rc = MPI_Type_commit(types + fold);
    if(rc != MPI_SUCCESS){
      if (rc == MPI_ERR_TYPE) {
        fprintf(stderr, "[%s.%d] ReproBLAS error: MPI_Type_commit error: MPI_ERR_TYPE\n", __FILE__, __LINE__);
      } else {
        fprintf(stderr, "[%s.%d] ReproBLAS error: MPI_Type_commit error: %d\n", __FILE__, __LINE__, rc);
      }
      MPI_Abort(MPI_COMM_WORLD, rc);
      return 0;
    }
    types_initialized[fold] = 1;
  }
  return types[fold];
}
