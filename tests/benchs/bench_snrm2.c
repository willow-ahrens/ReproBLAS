#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_perf.h"
#include "../common/test_BLAS.h"

#include "bench_vecvec_fill_header.h"

int bench_vecvec_fill_desc(void){
  char *op_names[] = {"s_add", "s_mul", "s_cmp", "s_orb"};
  int op_counts[] = {1, 2, 1, 1};
  perf_output_desc(4, op_names, op_counts);
  return 0;
}

int bench_vecvec_fill_show_help(void){
  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [snrm2]");
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY, int trials){
  int rc = 0;
  float res;

  util_random_seed();

  float *X = util_svec_alloc(N, incX);

  //fill X
  util_svec_fill(N, X, incX, FillX, ScaleX, CondX);

  time_tic();
  for(int i = 0; i < trials; i++){
    CALL_SNRM2(res, N, X, incX);
  }
  time_toc();

  perf_output_perf(time_read(), N, trials);

  free(X);
  return rc;
}
