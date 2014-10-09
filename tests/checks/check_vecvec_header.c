#include <stdio.h>
#include "../common/test_opt.h"
#define MAX_NAME 100

extern const char *vecvec_name(int argc, char** argv);
extern int vecvec_check(int argc, char** argv, int N, int incx, int incy);

static char namebuf[MAX_NAME];

const char* name(int argc, char** argv){
  int N = opt_read_int(argc, argv, "-N", 256);
  int incx = opt_read_int(argc, argv, "-incx", 1);
  int incy = opt_read_int(argc, argv, "-incy", 1);
  snprintf(namebuf, MAX_NAME, "%s N=%d ix=%d iy=%d", vecvec_name(argc, argv), N, incx, incy);
  return namebuf;
}

int check(int argc, char** argv){
  int N = opt_read_int(argc, argv, "-N", 256);
  int incx = opt_read_int(argc, argv, "-incx", 1);
  int incy = opt_read_int(argc, argv, "-incy", 1);
  int rc = vecvec_check(argc, argv, N, incx, incy);
  return rc;
}
