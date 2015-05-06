#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [zdotu]");
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY, int trials){
  (void)argc;
  (void)argv;
  int rc = 0;
  int i;
  double complex res = 0.0;

  util_random_seed();

  double complex *X = util_zvec_alloc(N, incX);
  double complex *Y = util_zvec_alloc(N, incY);

  //fill X and Y
  util_zvec_fill(N, X, incX, FillX, ScaleX, CondX);
  util_zvec_fill(N, Y, incY, FillY, ScaleY, CondY);

  time_tic();
  for(i = 0; i < trials; i++){
    CALL_ZDOTU(res, N, X, incX, Y, incY);
  }
  time_toc();

  metric_load_double("time", time_read());
  metric_load_double("res_real", creal(res));
  metric_load_double("res_imag", cimag(res));
  metric_load_long_long("trials", (long long)trials);
  metric_load_long_long("input", (long long)2 * N);
  metric_load_long_long("output", (long long)1);
  metric_load_long_long("d_fma", (long long)4 * N);
  metric_dump();

  free(X);
  free(Y);
  return rc;
}
