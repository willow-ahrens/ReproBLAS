#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <reproBLAS.h>
#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_perf.h"

#include "bench_vecvec_fill_header.h"

int bench_vecvec_fill_show_help(void){
  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rdznrm2]");
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY, int trials){
  int rc = 0;
  double complex res;

  util_random_seed();

  double complex *X = util_zvec_alloc(N, incX);

  //fill X
  util_zvec_fill(N, X, incX, FillX, ScaleX, CondX);

  time_tic();
  for(int i = 0; i < trials; i++){
    res = rdznrm2(N, X, incX);
  }
  time_toc();

  //TODO make generic fold testing
  int fold = 3;
  metric_load_double("time", time_read());
  metric_load_int("trials", trials);
  metric_load_int("input", N);
  metric_load_int("output", 1);
  metric_load_int("d_mul", 4 * N);
  metric_load_int("d_add", (3 * fold - 2) * 2 * N);
  metric_load_int("d_or", fold * 2 * N);
  metric_dump();

  free(X);
  return rc;
}
