#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <reproBLAS.h>
#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_perf.h"

#include "bench_vecvec_fill_header.h"

int bench_vecvec_fill_desc(void){
  char *op_names[] = {"s_mul", "s_add", "s_orb"};
  int op_counts[] = {4, 28, 12};
  perf_output_desc(3, op_names, op_counts);
  return 0;
}

int bench_vecvec_fill_show_help(void){
  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rcdotc]");
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY, int trials){
  int rc = 0;
  float complex res;

  util_random_seed();

  float complex *X = util_cvec_alloc(N, incX);
  float complex *Y = util_cvec_alloc(N, incY);

  //fill X and Y
  util_cvec_fill(N, X, incX, FillX, ScaleX, CondX);
  util_cvec_fill(N, Y, incY, FillY, ScaleY, CondY);

  time_tic();
  for(int i = 0; i < trials; i++){
    res = rcdotc(N, X, incX, Y, incY);
  }
  time_toc();

  perf_output_perf(time_read(), N, trials);

  free(X);
  free(Y);
  return rc;
}
