#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <indexedBLAS.h>
#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "bench_vecvec_fill_header.h"

int bench_vecvec_fill_show_help(void){
  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [camaxm]");
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY, int trials){
  int rc = 0;
  int i;

  float complex res = 0.0;

  util_random_seed();

  float complex *X = util_cvec_alloc(N, incX);
  float complex *Y = util_cvec_alloc(N, incY);

  //fill X and Y
  util_cvec_fill(N, X, incX, FillX, ScaleX, CondX);
  util_cvec_fill(N, Y, incY, FillY, ScaleY, CondY);

  time_tic();
  for(i = 0; i < trials; i++){
    camaxm_sub(N, X, incX, Y, incY, &res);
  }
  time_toc();

  metric_load_double("time", time_read());
  metric_load_float("res_real", crealf(res));
  metric_load_float("res_imag", cimagf(res));
  metric_load_long_long("trials", (long long)trials);
  metric_load_long_long("input", (long long)2 * N);
  metric_load_long_long("output", (long long)1);
  metric_load_long_long("s_mul", (long long)4 * N);
  metric_load_long_long("s_cmp", (long long)4 * N);
  metric_load_long_long("s_orb", (long long)4 * N);
  metric_dump();

  free(X);
  free(Y);
  return rc;
}
