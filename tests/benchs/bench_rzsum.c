#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"
#include "../common/test_time.h"
#include "../common/test_perf.h"
#include <rblas.h>
#include <IndexedFP.h>

#include "vecvec_fill_bench_header.h"

#define FLOP_PER_N 14

extern const char* vecvec_fill_bench_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rzsum]");
  return name_buffer;
}

extern int vecvec_fill_bench_test(int argc, char** argv, int N, int incx, int incy, int type, int perf_unit, int trials){
  int rc = 0;
  double complex res;
  I_double_Complex Ires;
  double complex *x = zvec_alloc(N, incx);

  vec_random_seed();

  //fill empty space with random data to check increment
  zvec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);

  //fill x
  zvec_fill(N, x, incx, type, 1.0, opt_read_float(argc, argv, "-c", 1.0));

  time_tic();
  for(int i = 0; i < trials; i++){
    res = rzsum(N, x, incx);
  }
  time_toc();

  printf("%e\n", perf_output(time_read(), N, trials, FLOP_PER_N, perf_unit, perf_prec_DOUBLE));

  free(x);
  return rc;
}
