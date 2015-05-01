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
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [zamaxm]");
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY, int trials){
  int rc = 0;
  double complex res;

  util_random_seed();

  double complex *X = util_zvec_alloc(N, incX);
  double complex *Y = util_zvec_alloc(N, incY);

  //fill X and Y
  util_zvec_fill(N, X, incX, FillX, ScaleX, CondX);
  util_zvec_fill(N, Y, incY, FillY, ScaleY, CondY);

  time_tic();
  for(int i = 0; i < trials; i++){
    res = zamaxm(N, X, incX, Y, incY);
  }
  time_toc();

  metric_load_double("time", time_read());
  metric_load_int("trials", trials);
  metric_load_int("input", 2 * N);
  metric_load_int("output", 1);
  metric_load_int("d_mul", 4 * N);
  metric_load_int("d_cmp", 4 * N);
  metric_load_int("d_orb", 4 * N);
  metric_dump();

  free(X);
  free(Y);
  return rc;
}
