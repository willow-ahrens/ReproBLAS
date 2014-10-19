#include <stdio.h>
#include "test_opt.h"

#include "test_header.h"

#define MAX_NAME 100

const char *vecvec_name(int argc, char** argv);
int vecvec_test(int argc, char** argv, int N, int incx, int incy);

const char* name(int argc, char** argv){
  static char namebuf[MAX_NAME];
  int N = opt_read_int(argc, argv, "-N", 256);
  int incx = opt_read_int(argc, argv, "-incx", 1);
  int incy = opt_read_int(argc, argv, "-incy", 1);
  snprintf(namebuf, MAX_NAME, "%s N=%d ix=%d iy=%d", vecvec_name(argc, argv), N, incx, incy);
  return namebuf;
}

int test(int argc, char** argv){
  int N = opt_read_int(argc, argv, "-N", 256);
  int incx = opt_read_int(argc, argv, "-incx", 1);
  int incy = opt_read_int(argc, argv, "-incy", 1);
  int rc = vecvec_test(argc, argv, N, incx, incy);
  return rc;
}
