#include <stdio.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"

#include "../common/test_vecvec_fill_header.h"
#include "unit_wrapper.h"

#define NAME_SIZE 100

const char* vecvec_fill_bench_name(int argc, char** argv);
int vecvec_fill_bench_test(int argc, char** argv, int N, int incx, int incy, int type, int unit, int trials);

const char* vecvec_fill_name(int argc, char** argv){
  static char namebuf[NAME_SIZE];
  int unit = opt_read_int(argc, argv, "-u", UNIT_HERTZ);
  int trials = opt_read_int(argc, argv, "-a", 1000);
  snprintf(namebuf, NAME_SIZE * sizeof(char), "%s (%s) (%d trials)", vecvec_fill_bench_name(argc, argv), unit_name(unit), trials);
  return namebuf;
}

int vecvec_fill_test(int argc, char** argv, int N, int incx, int incy, int type){
  int unit = opt_read_int(argc, argv, "-u", UNIT_HERTZ);
  int trials = opt_read_int(argc, argv, "-a", 1000);
  int rc = vecvec_fill_bench_test(argc, argv, N, incx, incy, type, unit, trials);
  return rc;
}
