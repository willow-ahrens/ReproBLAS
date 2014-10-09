#include <stdio.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"

#define NAME_SIZE 100

extern const char* vecvec_fill_name(int argc, char** argv);
extern int vecvec_fill_check(int argc, char** argv, int N, int incx, int incy, int type);

static char namebuf[NAME_SIZE];

const char* vecvec_name(int argc, char** argv){
  int type = opt_read_int(argc, argv, "-t", 0);
  snprintf(namebuf, NAME_SIZE * sizeof(char), "%s (%s)", vecvec_fill_name(argc, argv), vec_fill_name(type));
  return namebuf;
}

int vecvec_check(int argc, char** argv, int N, int incx, int incy){
  int type = opt_read_int(argc, argv, "-t", 0);
  int rc = vecvec_fill_check(argc, argv, N, incx, incy, type);
  return rc;
}
