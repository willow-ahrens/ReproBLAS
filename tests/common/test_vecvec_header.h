#include <stdio.h>
#include "test_opt.h"

#include "test_header.h"

const char *vecvec_name(int argc, char** argv);
int vecvec_test(int argc, char** argv, int N, int incX, int incY);

const char* name(int argc, char** argv){
  opt_option N;
  opt_option incX;
  opt_option incY;
  static char name_buffer[MAX_LINE];

  N.header.type       = opt_int;
  N.header.short_name = 'N';
  N.header.long_name  = "N_dim";
  N.header.help       = "N dimension size";
  N._int.required     = 0;
  N._int.min          = 0;
  N._int.max          = INT_MAX;
  N._int.value        = 2048;

  incX.header.type       = opt_int;
  incX.header.short_name = 'x';
  incX.header.long_name  = "incX";
  incX.header.help       = "X vector increment";
  incX._int.required     = 0;
  incX._int.min          = 1;
  incX._int.max          = INT_MAX;
  incX._int.value        = 1;

  incY.header.type       = opt_int;
  incY.header.short_name = 'y';
  incY.header.long_name  = "incY";
  incY.header.help       = "Y vector increment";
  incY._int.required     = 0;
  incY._int.min          = 1;
  incY._int.max          = INT_MAX;
  incY._int.value        = 1;

  if(help._flag.exists){
    opt_show_option(N);
    opt_show_option(incX);
    opt_show_option(incY);
    vecvec_name(argc, argv);
    return "";
  }

  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &incX);
  opt_eval_option(argc, argv, &incY);
  snprintf(name_buffer, MAX_LINE, "%s N=%d incX=%d incY=%d", vecvec_name(argc, argv), N._int.value, incX._int.value, incY._int.value);
  return name_buffer;
}

int test(int argc, char** argv){
  opt_option N;
  opt_option incX;
  opt_option incY;
  int rc;

  N.header.type       = opt_int;
  N.header.short_name = 'N';
  N.header.long_name  = "N_dim";
  N.header.help       = "N dimension size";
  N._int.required     = 0;
  N._int.min          = 0;
  N._int.max          = INT_MAX;
  N._int.value        = 2048;

  incX.header.type       = opt_int;
  incX.header.short_name = 'x';
  incX.header.long_name  = "incX";
  incX.header.help       = "X vector increment";
  incX._int.required     = 0;
  incX._int.min          = 1;
  incX._int.max          = INT_MAX;
  incX._int.value        = 1;

  incY.header.type       = opt_int;
  incY.header.short_name = 'y';
  incY.header.long_name  = "incY";
  incY.header.help       = "Y vector increment";
  incY._int.required     = 0;
  incY._int.min          = 1;
  incY._int.max          = INT_MAX;
  incY._int.value        = 1;

  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &incX);
  opt_eval_option(argc, argv, &incY);

  rc = vecvec_test(argc, argv, N._int.value, incX._int.value, incY._int.value);
  return rc;
}
