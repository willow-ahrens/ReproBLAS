#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_perf.h"
#include "../common/blas_inc.h"
#include <rblas.h>
#include <IndexedFP.h>

#include "vecvec_fill_bench_header.h"

int vecvec_fill_bench_desc(void){
  char *op_names[] = {"s_add", "s_mul"};
  int op_counts[] = {4, 4};
  perf_output_desc(2, op_names, op_counts);
  return 0;
}

int vecvec_fill_bench_show_help(void){
  return 0;
}

const char* vecvec_fill_bench_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [cdotu]");
  return name_buffer;
}

int vecvec_fill_bench_test(int argc, char** argv, int N, int incx, int incy, int type, double scale, double cond, int trials){
  int rc = 0;
  float complex res;
  I_float_Complex Ires;
  float complex *x = cvec_alloc(N, incx);
  float complex *y = cvec_alloc(N, incy);

  util_random_seed();

  //fill empty space with random data to check increment
  cvec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  cvec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill x
  cvec_fill(N, x, incx, type, scale, cond);

  //fill y with -i where necessary
  cvec_fill(N, y, incy, vec_fill_CONSTANT, -_Complex_I, 1.0);

  time_tic();
  for(int i = 0; i < trials; i++){
    CALL_CDOTU(res, N, x, incx, y, incy);
  }
  time_toc();

  perf_output_perf(time_read(), N, trials);

  free(x);
  free(y);
  return rc;
}
