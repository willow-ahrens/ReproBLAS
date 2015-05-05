#include <stdio.h>
#include "test_opt.h"

#include "test_matvec_header.h"

int matvec_fill_show_help(void);
const char *matvec_fill_name(int argc, char** argv);
int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, int FillA, double ScaleA, double CondA, int lda, int FillX, double ScaleX, double CondX, int incX);

static opt_option FillA;
static opt_option ScaleA;
static opt_option CondA;
static opt_option FillX;
static opt_option ScaleX;
static opt_option CondX;

static void matvec_fill_options_initialize(void){
  FillA._named.header.type       = opt_named;
  FillA._named.header.short_name = 'f';
  FillA._named.header.long_name  = "FillA";
  FillA._named.header.help       = "A fill type";
  FillA._named.required          = 0;
  FillA._named.n_names           = (int)util_mat_fill_n_names;
  FillA._named.names             = (char**)util_mat_fill_names;
  FillA._named.descs             = (char**)util_mat_fill_descs;
  FillA._named.value             = 0;

  ScaleA._double.header.type       = opt_double;
  ScaleA._double.header.short_name = 's';
  ScaleA._double.header.long_name  = "ScaleA";
  ScaleA._double.header.help       = "A scale";
  ScaleA._double.required          = 0;
  ScaleA._double.min               = 0;
  ScaleA._double.max               = DBL_MAX;
  ScaleA._double.value             = 1.0;

  CondA._double.header.type       = opt_double;
  CondA._double.header.short_name = 'c';
  CondA._double.header.long_name  = "CondA";
  CondA._double.header.help       = "A condition number";
  CondA._double.required          = 0;
  CondA._double.min               = 1.0;
  CondA._double.max               = DBL_MAX;
  CondA._double.value             = 1e3;

  FillX._named.header.type       = opt_named;
  FillX._named.header.short_name = 'g';
  FillX._named.header.long_name  = "FillX";
  FillX._named.header.help       = "X fill type";
  FillX._named.required          = 0;
  FillX._named.n_names           = (int)util_vec_fill_n_names;
  FillX._named.names             = (char**)util_vec_fill_names;
  FillX._named.descs             = (char**)util_vec_fill_descs;
  FillX._named.value             = 0;

  ScaleX._double.header.type       = opt_double;
  ScaleX._double.header.short_name = 't';
  ScaleX._double.header.long_name  = "ScaleX";
  ScaleX._double.header.help       = "X scale";
  ScaleX._double.required          = 0;
  ScaleX._double.min               = 0;
  ScaleX._double.max               = DBL_MAX;
  ScaleX._double.value             = 1.0;

  CondX._double.header.type       = opt_double;
  CondX._double.header.short_name = 'd';
  CondX._double.header.long_name  = "CondX";
  CondX._double.header.help       = "X condition number";
  CondX._double.required          = 0;
  CondX._double.min               = 1.0;
  CondX._double.max               = DBL_MAX;
  CondX._double.value             = 1e3;
}

int matvec_show_help(void){
  matvec_fill_options_initialize();

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

  matvec_fill_options_initialize();

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

  matvec_fill_options_initialize();

  opt_eval_option(argc, argv, &FillA);
  opt_eval_option(argc, argv, &ScaleA);
  opt_eval_option(argc, argv, &CondA);
  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &ScaleX);
  opt_eval_option(argc, argv, &CondX);
  rc = matvec_fill_test(argc, argv, Order, TransA, M, N, FillA._named.value, ScaleA._double.value, CondA._double.value, lda, FillX._named.value, ScaleX._double.value, CondX._double.value, incX);
  return rc;
}
