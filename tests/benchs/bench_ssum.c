#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "bench_vecvec_fill_header.h"

int bench_vecvec_fill_show_help(void){
  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [ssum]");
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
  float res = 0.0;

  util_random_seed();

  float *X = util_svec_alloc(N, incX);

  //fill X
  util_svec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);

  time_tic();
  for(i = 0; i < trials; i++){
    res = 0;
    for(j = 0; j < N; j++){
      res += X[j * incX];
    }
  }
  time_toc();

  double dN = (double)N;
  metric_load_double("time", time_read());
  metric_load_float("res", res);
  metric_load_double("trials", (double)trials);
  metric_load_double("input", dN);
  metric_load_double("output", 1.0);
  metric_load_double("normalizer", dN);
  metric_load_double("s_add", dN);
  metric_dump();

  free(X);
  return rc;
}
