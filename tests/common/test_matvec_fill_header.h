#include <stdio.h>
#include "test_opt.h"

#include "test_matvec_header.h"

int matvec_fill_show_help(void);
const char *matvec_fill_name(int argc, char** argv);
int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, int FillA, double ScaleA, double CondA, int lda, int FillX, double ScaleX, double CondX, int incX);

static opt_option FillA         = {._named.header.type       = opt_named,
                                   ._named.header.short_name = 'f',
                                   ._named.header.long_name  = "FillA",
                                   ._named.header.help       = "A fill type",
                                   ._named.required          = 0,
                                   ._named.n_names           = (int)util_mat_fill_n_names,
                                   ._named.names             = (char**)util_mat_fill_names,
                                   ._named.descs             = (char**)util_mat_fill_descs,
                                   ._named.value             = 0};

static opt_option ScaleA        = {._double.header.type       = opt_double,
                                   ._double.header.short_name = 's',
                                   ._double.header.long_name  = "ScaleA",
                                   ._double.header.help       = "A scale",
                                   ._double.required          = 0,
                                   ._double.min               = 0,
                                   ._double.max               = DBL_MAX,
                                   ._double.value             = 1.0};

static opt_option CondA         = {._double.header.type       = opt_double,
                                   ._double.header.short_name = 'c',
                                   ._double.header.long_name  = "CondA",
                                   ._double.header.help       = "A condition number",
                                   ._double.required          = 0,
                                   ._double.min               = 1.0,
                                   ._double.max               = DBL_MAX,
                                   ._double.value             = 1e3};

static opt_option FillX         = {._named.header.type       = opt_named,
                                   ._named.header.short_name = 'g',
                                   ._named.header.long_name  = "FillX",
                                   ._named.header.help       = "X fill type",
                                   ._named.required          = 0,
                                   ._named.n_names           = (int)util_vec_fill_n_names,
                                   ._named.names             = (char**)util_vec_fill_names,
                                   ._named.descs             = (char**)util_vec_fill_descs,
                                   ._named.value             = 0};

static opt_option ScaleX        = {._double.header.type       = opt_double,
                                   ._double.header.short_name = 't',
                                   ._double.header.long_name  = "ScaleX",
                                   ._double.header.help       = "X scale",
                                   ._double.required          = 0,
                                   ._double.min               = 0,
                                   ._double.max               = DBL_MAX,
                                   ._double.value             = 1.0};

static opt_option CondX         = {._double.header.type       = opt_double,
                                   ._double.header.short_name = 'd',
                                   ._double.header.long_name  = "CondX",
                                   ._double.header.help       = "X condition number",
                                   ._double.required          = 0,
                                   ._double.min               = 1.0,
                                   ._double.max               = DBL_MAX,
                                   ._double.value             = 1e3};

int matvec_show_help(void){
  opt_show_option(FillA);
  opt_show_option(ScaleA);
  opt_show_option(CondA);
  opt_show_option(FillX);
  opt_show_option(ScaleX);
  opt_show_option(CondX);
  return matvec_fill_show_help();
}

const char* matvec_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &FillA);
  opt_eval_option(argc, argv, &ScaleA);
  opt_eval_option(argc, argv, &CondA);
  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &ScaleX);
  opt_eval_option(argc, argv, &CondX);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%s) (%s)", matvec_fill_name(argc, argv), util_mat_fill_names[FillA._named.value], util_vec_fill_names[FillX._named.value]);
  return name_buffer;
}

int matvec_test(int argc, char** argv, char Order, char TransA, int M, int N, int lda, int incX){
  int rc;

  opt_eval_option(argc, argv, &FillA);
  opt_eval_option(argc, argv, &ScaleA);
  opt_eval_option(argc, argv, &CondA);
  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &ScaleX);
  opt_eval_option(argc, argv, &CondX);
  rc = matvec_fill_test(argc, argv, Order, TransA, M, N, FillA._named.value, ScaleA._double.value, CondA._double.value, lda, FillX._named.value, ScaleX._double.value, CondX._double.value, incX);
  return rc;
}
