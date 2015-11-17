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
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [sdot]");
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials){
  (void)argc;
  (void)argv;
  int rc = 0;
  int i;
  float res = 0.0;

  util_random_seed();

  float *X = util_svec_alloc(N, incX);
  float *Y = util_svec_alloc(N, incY);

  //fill X and Y
  util_svec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
  util_svec_fill(N, Y, incY, FillY, RealScaleY, ImagScaleY);

  time_tic();
  for(i = 0; i < trials; i++){
    CALL_SDOT(res, N, X, incX, Y, incY);
  }
  time_toc();

  double dN = (double)N;
  metric_load_double("time", time_read());
  metric_load_float("res", res);
  metric_load_double("trials", (double)trials);
  metric_load_double("input", 2.0 * dN);
  metric_load_double("output", 1.0);
  metric_load_double("normalizer", dN);
  metric_load_double("s_fma", dN);
  metric_dump();

  free(X);
  free(Y);
  return rc;
}
