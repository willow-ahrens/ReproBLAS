#include <stdio.h>
#include "test_opt.h"
#include "test_util.h"

#include "test_header.h"

int vecvec_show_help(void);
const char *vecvec_name(int argc, char** argv);
int vecvec_test(int argc, char** argv, int N, int incX, int incY);

static opt_option N;
static opt_option incX;
static opt_option incY;

static void vecvec_options_initialize(void){
  N._int.header.type       = opt_int;
  N._int.header.short_name = 'N';
  N._int.header.long_name  = "N_dim";
  N._int.header.help       = "N dimension size";
  N._int.required          = 0;
  N._int.min               = 0;
  N._int.max               = INT_MAX;
  N._int.value             = 2048;

  incX._int.header.type       = opt_int;
  incX._int.header.short_name = '\0';
  incX._int.header.long_name  = "incX";
  incX._int.header.help       = "X vector increment";
  incX._int.required          = 0;
  incX._int.min               = 1;
  incX._int.max               = INT_MAX;
  incX._int.value             = 1;

  incY._int.header.type       = opt_int;
  incY._int.header.short_name = '\0';
  incY._int.header.long_name  = "incY";
  incY._int.header.help       = "Y vector increment";
  incY._int.required          = 0;
  incY._int.min               = 1;
  incY._int.max               = INT_MAX;
  incY._int.value             = 1;
}

int show_help(void){
  vecvec_options_initialize();

  opt_show_option(N);
  opt_show_option(incX);
  opt_show_option(incY);
  return vecvec_show_help();
}

const char* name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  vecvec_options_initialize();

  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &incX);
  opt_eval_option(argc, argv, &incY);
  snprintf(name_buffer, MAX_LINE, "%s N=%d incX=%d incY=%d", vecvec_name(argc, argv), N._int.value, incX._int.value, incY._int.value);
  return name_buffer;
}

int test(int argc, char** argv){
  int rc;

  vecvec_options_initialize();

  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &incX);
  opt_eval_option(argc, argv, &incY);

  rc = vecvec_test(argc, argv, N._int.value, incX._int.value, incY._int.value);
  return rc;
}
