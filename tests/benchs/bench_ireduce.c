#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <idxd.h>
#include <MPI_idxd.h>
#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"
#include "../common/test_BLAS.h"

#include "bench_vecvec_fill_header.h"
    
int bench_vecvec_fill_show_help(void){
  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [reduce]");
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials){
  (void)argc;
  (void)argv;
  (void)FillY;
  (void)RealScaleY;
  (void)ImagScaleY;
  (void)incY;
  int rc = 0;
  int i, j;

  util_random_seed();

  int nprocs;
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  RMPI_Init();

  double *X = util_dvec_alloc(N, incX);
  double *Y = util_dvec_alloc(N, incX);
  double *IY = util_dvec_alloc(N * dinum(DEFAULT_FOLD), incX);
  double *IX = util_dvec_alloc(N * dinum(DEFAULT_FOLD), incX);

  //fill X
  util_dvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);

  time_tic();
  for(i = 0; i < trials; i++){
    for (j = 0; j < N; j++){
      didconv(DEFAULT_FOLD, X[j], IX + j * dinum(DEFAULT_FOLD));
    }
    MPI_Reduce(IX, IY, N, MPI_IDOUBLE, MPI_RSUM, 0, MPI_COMM_WORLD);
    if(rank == 0){
      for(j = 0; j < N; j++){
        Y[j] = ddiconv(DEFAULT_FOLD, IY + j * dinum(DEFAULT_FOLD));
      }
    }
  }
  time_toc();

  if(rank == 0){
    metric_load_double("time", time_read());
    metric_load_long_long("trials", (long long)trials);
    metric_load_long_long("input", (long long)1 * N);
    metric_load_long_long("output", (long long)1);
    metric_load_long_long("d_add", (long long)N);
    metric_load_long_long("d_orb", (long long)N);
    metric_dump();
  }

  MPI_Finalize();

  free(X);
  free(Y);
  return rc;
}
