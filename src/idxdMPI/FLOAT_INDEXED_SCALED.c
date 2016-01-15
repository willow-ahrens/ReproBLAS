#include <mpi.h>
#include <stdio.h>

#include <idxd.h>
#include <idxdMPI.h>

#include "../../config.h"

static MPI_Datatype types[idxd_SIMAXFOLD + 1];
static int types_initialized[idxd_SIMAXFOLD + 1]; //initializes to 0

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
