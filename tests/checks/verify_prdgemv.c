#include <indexedBLAS.h>
#include <MPI_reproBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include "../common/test_opt.h"

#include "../common/test_matvec_fill_header.h"

static opt_option incY   = {._int.header.type       = opt_int,
                            ._int.header.short_name = 'y',
                            ._int.header.long_name  = "incY",
                            ._int.header.help       = "Y vector increment",
                            ._int.required          = 0,
                            ._int.min               = 1,
                            ._int.max               = INT_MAX,
                            ._int.value             = 1};

static opt_option FillY  = {._named.header.type       = opt_named,
                            ._named.header.short_name = 'j',
                            ._named.header.long_name  = "FillY",
                            ._named.header.help       = "Y fill type",
                            ._named.required          = 0,
                            ._named.n_names           = (int)util_vec_fill_n_names,
                            ._named.names             = (char**)util_vec_fill_names,
                            ._named.descs             = (char**)util_vec_fill_descs,
                            ._named.value             = 0};

static opt_option ScaleY = {._double.header.type       = opt_double,
                            ._double.header.short_name = 'v',
                            ._double.header.long_name  = "ScaleY",
                            ._double.header.help       = "Y scale",
                            ._double.required          = 0,
                            ._double.min               = 0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1.0};

static opt_option CondY  = {._double.header.type       = opt_double,
                            ._double.header.short_name = 'e',
                            ._double.header.long_name  = "CondY",
                            ._double.header.help       = "Y condition number",
                            ._double.required          = 0,
                            ._double.min               = 1.0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1e3};

static opt_option alpha  = {._double.header.type       = opt_double,
                            ._double.header.short_name = 'l',
                            ._double.header.long_name  = "alpha",
                            ._double.header.help       = "alpha",
                            ._double.required          = 0,
                            ._double.min               = 1.0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1e3};

static opt_option beta   = {._double.header.type       = opt_double,
                            ._double.header.short_name = 'm',
                            ._double.header.long_name  = "beta",
                            ._double.header.help       = "beta",
                            ._double.required          = 0,
                            ._double.min               = 1.0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1e3};

void wrap_prdgemv(int rank, int nprocs, const char Order,
                  const char TransA, const int M, const int N,
                  const double *A, const double alpha, const int lda,
                  const double *X, const int incX,
                  const double beta, const double *Y, const int incY){
  rblas_order_t o;
  rblas_transpose_t t;
  switch(Order){
    case 'R':
      o = rblas_Row_Major;
      break;
    default:
      o = rblas_Col_Major;
      break;
  }
  switch(TransA){
    case 'N':
      t = rblas_No_Trans;
      break;
    default:
      t = rblas_Trans;
      break;
  }
  double *Abuf;
  double *myA = (double*)malloc(M*N/nprocs*sizeof(double));
  double *myX = (double*)malloc(N/nprocs*sizeof(double));
  if(rank == 0){
    Abuf = (double*)malloc(N*M*sizeof(double));
    for(int p = 0; p < nprocs; p++){
      for(int i = 0; i < M; i++){
        for(int j = 0; j < N/nprocs; j++){
          Abuf[p * M * N/nprocs + i * N/nprocs + j] = A[i * lda + (p * N/nprocs + j)];
        }
      }
    }
  }
  MPI_Scatter(Abuf, M * N/nprocs, MPI_DOUBLE, myA, M * N/nprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(X, N/nprocs, MPI_DOUBLE, myX, M * N/nprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  prdgemv(rank, nprocs, o, t, M, N, A, lda, X, incX, Y, incY);
}

int verify_dgemv_reproducibility(char Order, char TransA, int M, int N, int NX, int NY, double alpha, double *A, int lda, double* X, int incX, double beta, double *Y, Idouble *YI, int incY, double *ref, Idouble *Iref, int max_num_blocks) {
  (void)NX;
  // GENERATE DATA
  int i;
  int num_blocks = 1;

  double *res;
  int rank;
  int nprocs;

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(rank == 0){
    //compute with unpermuted data
    res = malloc(NY * incY * sizeof(double));
    memcpy(res, Y, NY * incY * sizeof(double));
  }
  wrap_prdgemv(rank, nprocs, Order, TransA, M, N, A, alpha, lda, X, incX, beta, res, incY);
  if(rank == 0){
    for(i = 0; i < NY; i++){
      if(res[i * incY] != ref[i * incY]){
        printf("dgemv(A, X, Y) = %g != %g\n", res[i * incY], ref[i * incY]);
        return 1;
      }
    }
  }
  return 0;
}

int matvec_fill_show_help(void){
  opt_show_option(incY);
  opt_show_option(FillY);
  opt_show_option(ScaleY);
  opt_show_option(CondY);
  opt_show_option(alpha);
  opt_show_option(beta);
  return 0;
}

const char* matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &incY);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &ScaleY);
  opt_eval_option(argc, argv, &CondY);
  opt_eval_option(argc, argv, &alpha);
  opt_eval_option(argc, argv, &beta);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify dgemv reproducibility");
  return name_buffer;
}

int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, int FillA, double ScaleA, double CondA, int lda, int FillX, double ScaleX, double CondX, int incX){
  int rc = 0;
  int max_num_blocks = 1024;

  opt_eval_option(argc, argv, &incY);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &ScaleY);
  opt_eval_option(argc, argv, &CondY);
  opt_eval_option(argc, argv, &alpha);
  opt_eval_option(argc, argv, &beta);

  util_random_seed();
  int NX;
  int NY;
  char NTransA;
  switch(TransA){
    case 'n':
    case 'N':
      NX = N;
      NY = M;
      NTransA = 'T';
    break;
    default:
      NX = M;
      NY = N;
      NTransA = 'N';
    break;
  }
  int nprocs;
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double *A;
  double *X;
  double *Y;
  double *ref;
  Idouble *Iref;
  Idouble *YI = NULL;

  if(rank == 0){
    A  = util_dmat_alloc(Order, M, N, lda);
    X  = util_dvec_alloc(NX, incX);
    Y  = util_dvec_alloc(NY, incY._int.value);
    util_dmat_fill(Order, 'n', M, N, A, lda, FillA, ScaleA, CondA);
    util_dvec_fill(NX, X, incX, FillX, ScaleX, CondX);
    util_dvec_fill(NY, Y, incY._int.value, FillY._named.value, ScaleY._double.value, CondY._double.value);
    ref  = (double*)malloc(NY * incY._int.value * sizeof(double));
    Iref = NULL;
    //compute with unpermuted data
    memcpy(ref, Y, NY * incY._int.value * sizeof(double));
    memcpy(Iref, YI, NY * incY._int.value * dISize(DEFAULT_FOLD));
  }

  int *P;

  wrap_prdgemv(rank, nprocs, Order, TransA, M, N, A, alpha._double.value, lda, X, incX, beta._double.value, ref, incY._int.value);

  if(rank == 0){
    P = util_identity_permutation(NX);
    util_dvec_reverse(NX, X, incX, P, 1);
    util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);
  }

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  if(rank == 0){
    P = util_identity_permutation(NX);
    util_dvec_sort(NX, X, incX, P, 1, util_Increasing);
    util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);
  }

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  if(rank == 0){
    P = util_identity_permutation(NX);
    util_dvec_sort(NX, X, incX, P, 1, util_Decreasing);
    util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);
  }

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  if(rank == 0){
    P = util_identity_permutation(NX);
    util_dvec_sort(NX, X, incX, P, 1, util_Increasing_Magnitude);
    util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);
  }

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  if(rank == 0){
    P = util_identity_permutation(NX);
    util_dvec_sort(NX, X, incX, P, 1, util_Decreasing_Magnitude);
    util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);
  }

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  if(rank == 0){
    P = util_identity_permutation(NX);
    util_dvec_shuffle(NX, X, incX, P, 1);
    util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);
  }

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  if(rank == 0){
    P = util_identity_permutation(NX);
    util_dvec_shuffle(NX, X, incX, P, 1);
    util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);
  }

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  if(rank == 0){
    P = util_identity_permutation(NX);
    util_dvec_shuffle(NX, X, incX, P, 1);
    util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);
  }

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  if(rank == 0){
    P = util_identity_permutation(NX);
    util_dvec_shuffle(NX, X, incX, P, 1);
    util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);
  }

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  free(A);
  free(X);
  free(Y);
  free(ref);
  free(Iref);

  return rc;
}
