#include <idxdBLAS.h>
#include <MPI_reproBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "bench_matvec_fill_header.h"

static opt_option incY;
static opt_option FillY;
static opt_option RealScaleY;
static opt_option ImagScaleY;
static opt_option alpha;
static opt_option beta;

void bench_prdgemv_options_initialize(void){
  incY._int.header.type       = opt_int;
  incY._int.header.short_name = 'y';
  incY._int.header.long_name  = "incY";
  incY._int.header.help       = "Y vector increment";
  incY._int.required          = 0;
  incY._int.min               = 1;
  incY._int.max               = INT_MAX;
  incY._int.value             = 1;

  FillY._named.header.type       = opt_named;
  FillY._named.header.short_name = 'j';
  FillY._named.header.long_name  = "FillY";
  FillY._named.header.help       = "Y fill type";
  FillY._named.required          = 0;
  FillY._named.n_names           = (int)util_vec_fill_n_names;
  FillY._named.names             = (char**)util_vec_fill_names;
  FillY._named.descs             = (char**)util_vec_fill_descs;
  FillY._named.value             = 0;

  RealScaleY._double.header.type       = opt_double;
  RealScaleY._double.header.short_name = 'v';
  RealScaleY._double.header.long_name  = "RealScaleY";
  RealScaleY._double.header.help       = "Y scale";
  RealScaleY._double.required          = 0;
  RealScaleY._double.min               = 0;
  RealScaleY._double.max               = DBL_MAX;
  RealScaleY._double.value             = 1.0;

  ImagScaleY._double.header.type       = opt_double;
  ImagScaleY._double.header.short_name = 'e';
  ImagScaleY._double.header.long_name  = "ImagScaleY";
  ImagScaleY._double.header.help       = "Y condition number";
  ImagScaleY._double.required          = 0;
  ImagScaleY._double.min               = 1.0;
  ImagScaleY._double.max               = DBL_MAX;
  ImagScaleY._double.value             = 1e3;

  alpha._double.header.type       = opt_double;
  alpha._double.header.short_name = 'l';
  alpha._double.header.long_name  = "alpha";
  alpha._double.header.help       = "alpha";
  alpha._double.required          = 0;
  alpha._double.min               = 1.0;
  alpha._double.max               = DBL_MAX;
  alpha._double.value             = 1e3;

  beta._double.header.type       = opt_double;
  beta._double.header.short_name = 'm';
  beta._double.header.long_name  = "beta";
  beta._double.header.help       = "beta";
  beta._double.required          = 0;
  beta._double.min               = 1.0;
  beta._double.max               = DBL_MAX;
  beta._double.value             = 1e3;
}

int bench_matvec_fill_show_help(void){
  bench_prdgemv_options_initialize();

  opt_show_option(incY);
  opt_show_option(FillY);
  opt_show_option(RealScaleY);
  opt_show_option(ImagScaleY);
  opt_show_option(alpha);
  opt_show_option(beta);
  return 0;
}

const char* bench_matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  bench_prdgemv_options_initialize();

  opt_eval_option(argc, argv, &incY);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &RealScaleY);
  opt_eval_option(argc, argv, &ImagScaleY);
  opt_eval_option(argc, argv, &alpha);
  opt_eval_option(argc, argv, &beta);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify dgemv reproducibility");
  return name_buffer;
}
extern void blacs_get_(int *, int* ,int*);
extern void blacs_gridinit_(int*, char*, int*, int*);
extern void blacs_gridinfo_(int *, int*, int*, int*, int*);

extern void pdgemv_(char*, int*, int*, double*, double*, int*, int*, int*, double*, int*, int*, int*, int*, double*, double*, int*, int*, int*, int*);

int bench_matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, int trials){
  int rc = 0;

  bench_prdgemv_options_initialize();

  opt_eval_option(argc, argv, &incY);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &RealScaleY);
  opt_eval_option(argc, argv, &ImagScaleY);
  opt_eval_option(argc, argv, &alpha);
  opt_eval_option(argc, argv, &beta);

  util_random_seed();
  int NX;
  int NY;
  switch(TransA){
    case 'n':
    case 'N':
      NX = N;
      NY = M;
    break;
    default:
      NX = M;
      NY = N;
    break;
  }
  int nprocs;
  int rank;
  int icontext;
  int nprow;
  int npcol;
  int myprow;
  int mypcol;
  int i0;
  int i1;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  i0 = 0;
  i1 = 0;
  blacs_get_(&i0, &i1, &icontext);
  nprow = 1;
  npcol = nprocs;
  blacs_gridinit_(&icontext, "Row", &nprow, &npcol);
  blacs_gridinfo_(&icontext, &nprow, &npcol, &myprow, &mypcol);

  double *A;
  double *X;
  double *Y;
  double *ref;

  if(rank == 0){
    A  = util_dmat_alloc(Order, M, N, lda);
    X  = util_dvec_alloc(NX, incX);
    Y  = util_dvec_alloc(NY, incY._int.value);
    util_dmat_fill(Order, 'n', M, N, A, lda, FillA, RealScaleA, ImagScaleA);
    util_dvec_fill(NX, X, incX, FillX, RealScaleX, ImagScaleX);
    util_dvec_fill(NY, Y, incY._int.value, FillY._named.value, RealScaleY._double.value, ImagScaleY._double.value);
    ref  = (double*)malloc(NY * incY._int.value * sizeof(double));
    //compute with unpermuted data
    memcpy(ref, Y, NY * incY._int.value * sizeof(double));
  }else{
    A = NULL;
    X = NULL;
    Y = NULL;
    ref = NULL;
  }

  bench_prdgemv_options_initialize();

  // GENERATE DATA
  int i;

  double *res;

  double *Abuf;
  double *myA = (double*)malloc(M*(N/nprocs)*sizeof(double));
  double *myX = (double*)malloc((N/nprocs)*sizeof(double));
  int p;
  int j;
  if(rank == 0){
    Abuf = (double*)malloc(N*M*sizeof(double));
    for(p = 0; p < nprocs; p++){
      for(i = 0; i < M; i++){
        for(j = 0; j < N/nprocs; j++){
          Abuf[p * M * N/nprocs + i * N/nprocs + j] = A[i * lda + (p * N/nprocs + j)];
        }
      }
    }
  }else{
    Abuf = 0;
  }
  MPI_Scatter(Abuf, M * N/nprocs, MPI_DOUBLE, myA, M * N/nprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(X, N/nprocs, MPI_DOUBLE, myX, N/nprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int descA[9];
  descA[0] = 1;
  descA[1] = icontext;
  descA[2] = N;
  descA[3] = M;
  descA[4] = N/nprocs;
  descA[5] = M;
  descA[6] = 0;
  descA[7] = 0;
  descA[8] = N;

  int descX[9];
  descX[0] = 1;
  descX[1] = icontext;
  descX[2] = N;
  descX[3] = 1;
  descX[4] = N/nprocs;
  descX[5] = 1;
  descX[6] = 0;
  descX[7] = 0;
  descX[8] = N;

  int descY[9];
  descY[0] = 1;
  descY[1] = icontext;
  descY[2] = M;
  descY[3] = 1;
  descY[4] = M;
  descY[5] = 1;
  descY[6] = 0;
  descY[7] = 0;
  descY[8] = M;
  
  char trans = 'T';
  int nn = N/nprocs;
  i0 = 0;
  i1 = 1;

  for(i = 0; i < trials; i++){
    if(rank == 0){
      //compute with unpermuted data
      res = malloc(NY * incY._int.value * sizeof(double));
      memcpy(res, Y, NY * incY._int.value * sizeof(double));
    }else{
      res = NULL;
    }
    time_tic();
    
    pdgemv_(&trans, &nn, &M, &alpha._double.value, myA, &i1, &i1, descA, myX, &i1, &i1, descX, &i1, &beta._double.value, Y, &i1, &i1, descY, &i1);
    time_toc();
  }

  if(rank == 0){
    //TODO make generic fold testing
    int fold = 3;
    metric_load_double("time", time_read());
    metric_load_long_long("trials", (long long)trials);
    metric_load_long_long("input", (long long)N * M + N + M);
    metric_load_long_long("output", (long long)NY);
    metric_load_long_long("d_mul", (long long)N * M);
    metric_load_long_long("d_add", (long long)(3 * fold - 2) * N * M);
    metric_load_long_long("d_orb", (long long)fold * N * M);
    metric_dump();
    free(Abuf);
  }
  free(myA);
  free(myX);

  if(rank == 0){
    free(A);
    free(X);
    free(Y);
    free(ref);
  }

  MPI_Finalize();

  return rc;
}
