#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_perf.h"
#include <rblas.h>
#include <IndexedFP.h>

#include "vecvec_fill_bench_header.h"

int vecvec_fill_bench_desc(void){
  char *op_names[] = {"s_add", "s_orb"};
  int op_counts[] = {14, 8};
  perf_output_desc(2, op_names, op_counts);
  return 0;
}

int vecvec_fill_bench_show_help(void){
  return 0;
}

const char* vecvec_fill_bench_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rscasum]");
  return name_buffer;
}

int vecvec_fill_bench_test(int argc, char** argv, int N, int incx, int incy, int type, double scale, double cond, int trials){
  int rc = 0;
  float complex res;
  I_float_Complex Ires;

  util_random_seed();

  float complex *x = cvec_alloc(N, incx);

  //fill x
  cvec_fill(N, x, incx, type, scale, cond);

  time_tic();
  for(int i = 0; i < trials; i++){
    res = rscasum(N, x, incx);
  }
  time_toc();

  perf_output_perf(time_read(), N, trials);

  free(x);
  return rc;
}
