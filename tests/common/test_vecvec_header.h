#include <stdio.h>
#include "test_opt.h"
#include "test_util.h"

#include "test_header.h"

int vecvec_show_help(void);
const char *vecvec_name(int argc, char** argv);
int vecvec_test(int argc, char** argv, int N, int incX, int incY);

static opt_option N    = {._int.header.type       = opt_int,
                          ._int.header.short_name = 'N',
                          ._int.header.long_name  = "N_dim",
                          ._int.header.help       = "N dimension size",
                          ._int.required          = 0,
                          ._int.min               = 0,
                          ._int.max               = INT_MAX,
                          ._int.value             = 2048};

static opt_option incX = {._int.header.type       = opt_int,
                          ._int.header.short_name = 'x',
                          ._int.header.long_name  = "incX",
                          ._int.header.help       = "X vector increment",
                          ._int.required          = 0,
                          ._int.min               = 1,
                          ._int.max               = INT_MAX,
                          ._int.value             = 1};

static opt_option incY = {._int.header.type       = opt_int,
                          ._int.header.short_name = 'y',
                          ._int.header.long_name  = "incY",
                          ._int.header.help       = "Y vector increment",
                          ._int.required          = 0,
                          ._int.min               = 1,
                          ._int.max               = INT_MAX,
                          ._int.value             = 1};

int show_help(void){
  opt_show_option(N);
  opt_show_option(incX);
  opt_show_option(incY);
  return vecvec_show_help();
}

const char* name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &incX);
  opt_eval_option(argc, argv, &incY);
  snprintf(name_buffer, MAX_LINE, "%s N=%d incX=%d incY=%d", vecvec_name(argc, argv), N._int.value, incX._int.value, incY._int.value);
  return name_buffer;
}

int test(int argc, char** argv){
  int rc;

  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &incX);
  opt_eval_option(argc, argv, &incY);

  rc = vecvec_test(argc, argv, N._int.value, incX._int.value, incY._int.value);
  return rc;
}
