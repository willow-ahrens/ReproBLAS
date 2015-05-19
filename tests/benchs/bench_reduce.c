#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
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

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY, int trials){
  (void)argc;
  (void)argv;
  (void)FillY;
  (void)ScaleY;
  (void)CondY;
  (void)incY;
  int rc = 0;
  int i;

  util_random_seed();

  int nprocs;
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double *X = util_dvec_alloc(N, incX);
  double *Y = util_dvec_alloc(N, incX);

  //fill X
  util_dvec_fill(N, X, incX, FillX, ScaleX, CondX);

  time_tic();
  for(i = 0; i < trials; i++){
    MPI_Reduce(X, Y, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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
