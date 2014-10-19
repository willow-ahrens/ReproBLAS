#include <stdio.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"
#include "test_vecvec.h"

#define NAME_SIZE 100

const char* vecvec_fill_name(int argc, char** argv);
int vecvec_fill_test(int argc, char** argv, int N, int incx, int incy, int type);

const char* vecvec_name(int argc, char** argv){
  static char namebuf[NAME_SIZE];
  int type = opt_read_int(argc, argv, "-t", 0);
  snprintf(namebuf, NAME_SIZE * sizeof(char), "%s (%s)", vecvec_fill_name(argc, argv), vec_fill_name(type));
  return namebuf;
}

int vecvec_test(int argc, char** argv, int N, int incx, int incy){
  int type = opt_read_int(argc, argv, "-t", 0);
  int rc = vecvec_fill_test(argc, argv, N, incx, incy, type);
  return rc;
}
