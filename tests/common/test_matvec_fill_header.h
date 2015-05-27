#include <stdio.h>
#include "test_opt.h"

#include "test_matvec_header.h"

int matvec_fill_show_help(void);
const char *matvec_fill_name(int argc, char** argv);
int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX);

static opt_option FillA;
static opt_option RealScaleA;
static opt_option ImagScaleA;
static opt_option FillX;
static opt_option RealScaleX;
static opt_option ImagScaleX;

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

  RealScaleA._double.header.type       = opt_double;
  RealScaleA._double.header.short_name = 's';
  RealScaleA._double.header.long_name  = "RealScaleA";
  RealScaleA._double.header.help       = "A scale (real)";
  RealScaleA._double.required          = 0;
  RealScaleA._double.min               = -1 * DBL_MAX;
  RealScaleA._double.max               = DBL_MAX;
  RealScaleA._double.value             = 1.0;

  ImagScaleA._double.header.type       = opt_double;
  ImagScaleA._double.header.short_name = 'c';
  ImagScaleA._double.header.long_name  = "ImagScaleA";
  ImagScaleA._double.header.help       = "A scale (imaginary)";
  ImagScaleA._double.required          = 0;
  ImagScaleA._double.min               = -1 * DBL_MAX;
  ImagScaleA._double.max               = DBL_MAX;
  ImagScaleA._double.value             = 0.0;

  FillX._named.header.type       = opt_named;
  FillX._named.header.short_name = 'g';
  FillX._named.header.long_name  = "FillX";
  FillX._named.header.help       = "X fill type";
  FillX._named.required          = 0;
  FillX._named.n_names           = (int)util_vec_fill_n_names;
  FillX._named.names             = (char**)util_vec_fill_names;
  FillX._named.descs             = (char**)util_vec_fill_descs;
  FillX._named.value             = 0;

  RealScaleX._double.header.type       = opt_double;
  RealScaleX._double.header.short_name = 't';
  RealScaleX._double.header.long_name  = "RealScaleX";
  RealScaleX._double.header.help       = "X scale (real)";
  RealScaleX._double.required          = 0;
  RealScaleX._double.min               = -1 * DBL_MAX;
  RealScaleX._double.max               = DBL_MAX;
  RealScaleX._double.value             = 1.0;

  ImagScaleX._double.header.type       = opt_double;
  ImagScaleX._double.header.short_name = 'd';
  ImagScaleX._double.header.long_name  = "ImagScaleX";
  ImagScaleX._double.header.help       = "X scale (imag)";
  ImagScaleX._double.required          = 0;
  ImagScaleX._double.min               = -1 * DBL_MAX;
  ImagScaleX._double.max               = DBL_MAX;
  ImagScaleX._double.value             = 0.0;
}

int matvec_show_help(void){
  matvec_fill_options_initialize();

  opt_show_option(FillA);
  opt_show_option(RealScaleA);
  opt_show_option(ImagScaleA);
  opt_show_option(FillX);
  opt_show_option(RealScaleX);
  opt_show_option(ImagScaleX);
  return matvec_fill_show_help();
}

const char* matvec_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  matvec_fill_options_initialize();

  opt_eval_option(argc, argv, &FillA);
  opt_eval_option(argc, argv, &RealScaleA);
  opt_eval_option(argc, argv, &ImagScaleA);
  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &RealScaleX);
  opt_eval_option(argc, argv, &ImagScaleX);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%s) (%s)", matvec_fill_name(argc, argv), util_mat_fill_names[FillA._named.value], util_vec_fill_names[FillX._named.value]);
  return name_buffer;
}

int matvec_test(int argc, char** argv, char Order, char TransA, int M, int N, int lda, int incX){
  int rc;

  matvec_fill_options_initialize();

  opt_eval_option(argc, argv, &FillA);
  opt_eval_option(argc, argv, &RealScaleA);
  opt_eval_option(argc, argv, &ImagScaleA);
  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &RealScaleX);
  opt_eval_option(argc, argv, &ImagScaleX);
  rc = matvec_fill_test(argc, argv, Order, TransA, M, N, FillA._named.value, RealScaleA._double.value, ImagScaleA._double.value, lda, FillX._named.value, RealScaleX._double.value, ImagScaleX._double.value, incX);
  return rc;
}
