#include <stdio.h>
#include "test_opt.h"

#include "test_vecvec_header.h"

int vecvec_fill_show_help(void);
const char* vecvec_fill_name(int argc, char** argv);
int vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY);

static opt_option FillX;
static opt_option ScaleX;
static opt_option CondX;
static opt_option FillY;
static opt_option ScaleY;
static opt_option CondY;

static void vecvec_fill_options_initialize(void){
  FillX._named.header.type       = opt_named;
  FillX._named.header.short_name = 'f';
  FillX._named.header.long_name  = "FillX";
  FillX._named.header.help       = "X fill type";
  FillX._named.required          = 0;
  FillX._named.n_names           = (int)util_vec_fill_n_names;
  FillX._named.names             = (char**)util_vec_fill_names;
  FillX._named.descs             = (char**)util_vec_fill_descs;
  FillX._named.value             = 0;

  ScaleX._double.header.type       = opt_double;
  ScaleX._double.header.short_name = 's';
  ScaleX._double.header.long_name  = "ScaleX";
  ScaleX._double.header.help       = "X scale";
  ScaleX._double.required          = 0;
  ScaleX._double.min               = 0;
  ScaleX._double.max               = DBL_MAX;
  ScaleX._double.value             = 1.0;

  CondX._double.header.type       = opt_double;
  CondX._double.header.short_name = 'c';
  CondX._double.header.long_name  = "CondX";
  CondX._double.header.help       = "X condition number";
  CondX._double.required          = 0;
  CondX._double.min               = 1.0;
  CondX._double.max               = DBL_MAX;
  CondX._double.value             = 1e3;

  FillY._named.header.type       = opt_named;
  FillY._named.header.short_name = 'g';
  FillY._named.header.long_name  = "FillY";
  FillY._named.header.help       = "Y fill type";
  FillY._named.required          = 0;
  FillY._named.n_names           = (int)util_vec_fill_n_names;
  FillY._named.names             = (char**)util_vec_fill_names;
  FillY._named.descs             = (char**)util_vec_fill_descs;
  FillY._named.value             = 0;

  ScaleY._double.header.type       = opt_double;
  ScaleY._double.header.short_name = 't';
  ScaleY._double.header.long_name  = "ScaleY";
  ScaleY._double.header.help       = "Y scale";
  ScaleY._double.required          = 0;
  ScaleY._double.min               = 0;
  ScaleY._double.max               = DBL_MAX;
  ScaleY._double.value             = 1.0;

  CondY._double.header.type       = opt_double;
  CondY._double.header.short_name = 'd';
  CondY._double.header.long_name  = "CondY";
  CondY._double.header.help       = "Y condition number";
  CondY._double.required          = 0;
  CondY._double.min               = 1.0;
  CondY._double.max               = DBL_MAX;
  CondY._double.value             = 1e3;
}

int vecvec_show_help(void){
  vecvec_fill_options_initialize();

  opt_show_option(FillX);
  opt_show_option(ScaleX);
  opt_show_option(CondX);
  opt_show_option(FillY);
  opt_show_option(ScaleY);
  opt_show_option(CondY);
  return vecvec_fill_show_help();
}

const char* vecvec_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  vecvec_fill_options_initialize();

  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &ScaleX);
  opt_eval_option(argc, argv, &CondX);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &ScaleY);
  opt_eval_option(argc, argv, &CondY);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s X=%s*%g Y=%s*%g", vecvec_fill_name(argc, argv), util_vec_fill_names[FillX._named.value], ScaleX._double.value, util_vec_fill_names[FillY._named.value], ScaleY._double.value);
  return name_buffer;
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY){
  int rc;

  vecvec_fill_options_initialize();

  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &ScaleX);
  opt_eval_option(argc, argv, &CondX);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &ScaleY);
  opt_eval_option(argc, argv, &CondY);
  rc = vecvec_fill_test(argc, argv, N, FillX._named.value, ScaleX._double.value, CondX._double.value, incX, FillY._named.value, ScaleY._double.value, CondY._double.value, incY);
  return rc;
}
