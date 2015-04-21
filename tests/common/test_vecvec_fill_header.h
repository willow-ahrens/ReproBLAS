#include <stdio.h>
#include "test_opt.h"

#include "test_vecvec_header.h"

int vecvec_fill_show_help(void);
const char* vecvec_fill_name(int argc, char** argv);
int vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY);

static opt_option FillX = {._named.header.type       = opt_named,
                           ._named.header.short_name = 'f',
                           ._named.header.long_name  = "FillX",
                           ._named.header.help       = "X fill type",
                           ._named.required          = 0,
                           ._named.n_names           = (int)util_vec_fill_n_names,
                           ._named.names             = (char**)util_vec_fill_names,
                           ._named.descs             = (char**)util_vec_fill_descs,
                           ._named.value             = 0};

static opt_option ScaleX = {._double.header.type       = opt_double,
                            ._double.header.short_name = 's',
                            ._double.header.long_name  = "ScaleX",
                            ._double.header.help       = "X scale",
                            ._double.required          = 0,
                            ._double.min               = 0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1.0};

static opt_option CondX  = {._double.header.type       = opt_double,
                            ._double.header.short_name = 'c',
                            ._double.header.long_name  = "CondX",
                            ._double.header.help       = "X condition number",
                            ._double.required          = 0,
                            ._double.min               = 1.0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1e3};

static opt_option FillY = {._named.header.type       = opt_named,
                           ._named.header.short_name = 'g',
                           ._named.header.long_name  = "FillY",
                           ._named.header.help       = "Y fill type",
                           ._named.required          = 0,
                           ._named.n_names           = (int)util_vec_fill_n_names,
                           ._named.names             = (char**)util_vec_fill_names,
                           ._named.descs             = (char**)util_vec_fill_descs,
                           ._named.value             = 0};

static opt_option ScaleY = {._double.header.type       = opt_double,
                            ._double.header.short_name = 't',
                            ._double.header.long_name  = "ScaleY",
                            ._double.header.help       = "Y scale",
                            ._double.required          = 0,
                            ._double.min               = 0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1.0};

static opt_option CondY  = {._double.header.type       = opt_double,
                            ._double.header.short_name = 'd',
                            ._double.header.long_name  = "CondY",
                            ._double.header.help       = "Y condition number",
                            ._double.required          = 0,
                            ._double.min               = 1.0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1e3};

int vecvec_show_help(void){
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

  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &ScaleX);
  opt_eval_option(argc, argv, &CondX);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &ScaleY);
  opt_eval_option(argc, argv, &CondY);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s X=%s Y=%s", vecvec_fill_name(argc, argv), util_vec_fill_names[FillX._named.value], util_vec_fill_names[FillY._named.value]);
  return name_buffer;
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY){
  int rc;

  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &ScaleX);
  opt_eval_option(argc, argv, &CondX);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &ScaleY);
  opt_eval_option(argc, argv, &CondY);
  rc = vecvec_fill_test(argc, argv, N, FillX._named.value, ScaleX._double.value, CondX._double.value, incX, FillY._named.value, ScaleY._double.value, CondY._double.value, incY);
  return rc;
}
