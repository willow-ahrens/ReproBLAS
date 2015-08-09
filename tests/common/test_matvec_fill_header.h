#include <stdio.h>
#include "test_opt.h"

#include "test_matvec_header.h"

int matvec_fill_show_help(void);
const char *matvec_fill_name(int argc, char** argv);
int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY);

static opt_option FillA;
static opt_option RealScaleA;
static opt_option ImagScaleA;
static opt_option FillX;
static opt_option RealScaleX;
static opt_option ImagScaleX;
static opt_option FillY;
static opt_option RealScaleY;
static opt_option ImagScaleY;
static opt_option RealAlpha;
static opt_option ImagAlpha;
static opt_option RealBeta;
static opt_option ImagBeta;

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
  ImagScaleX._double.header.help       = "X scale (imaginary)";
  ImagScaleX._double.required          = 0;
  ImagScaleX._double.min               = -1 * DBL_MAX;
  ImagScaleX._double.max               = DBL_MAX;
  ImagScaleX._double.value             = 0.0;

  FillY._named.header.type       = opt_named;
  FillY._named.header.short_name = 'j';
  FillY._named.header.long_name  = "FillY";
  FillY._named.header.help       = "Y fill type";
  FillY._named.required          = 0;
  FillY._named.n_names           = (int)util_vec_fill_n_names;
  FillY._named.names             = (char**)util_vec_fill_names;
  FillY._named.descs             = (char**)util_vec_fill_descs;
  FillY._named.value             = 0;

  RealScaleY._double.header.type       = opt_double;
  RealScaleY._double.header.short_name = 'v';
  RealScaleY._double.header.long_name  = "RealScaleY";
  RealScaleY._double.header.help       = "Y scale (real)";
  RealScaleY._double.required          = 0;
  RealScaleY._double.min               = -1*DBL_MAX;
  RealScaleY._double.max               = DBL_MAX;
  RealScaleY._double.value             = 1.0;

  ImagScaleY._double.header.type       = opt_double;
  ImagScaleY._double.header.short_name = 'e';
  ImagScaleY._double.header.long_name  = "ImagScaleY";
  ImagScaleY._double.header.help       = "Y scale (imaginary)";
  ImagScaleY._double.required          = 0;
  ImagScaleY._double.min               = -1*DBL_MAX;
  ImagScaleY._double.max               = DBL_MAX;
  ImagScaleY._double.value             = 0;

  RealAlpha._double.header.type       = opt_double;
  RealAlpha._double.header.short_name = 'l';
  RealAlpha._double.header.long_name  = "RealAlpha";
  RealAlpha._double.header.help       = "alpha (real)";
  RealAlpha._double.required          = 0;
  RealAlpha._double.min               = -DBL_MAX;
  RealAlpha._double.max               = DBL_MAX;
  RealAlpha._double.value             = 1.0;

  ImagAlpha._double.header.type       = opt_double;
  ImagAlpha._double.header.short_name = 'l';
  ImagAlpha._double.header.long_name  = "ImagAlpha";
  ImagAlpha._double.header.help       = "alpha (imaginary)";
  ImagAlpha._double.required          = 0;
  ImagAlpha._double.min               = -DBL_MAX;
  ImagAlpha._double.max               = DBL_MAX;
  ImagAlpha._double.value             = 0.0;

  RealBeta._double.header.type       = opt_double;
  RealBeta._double.header.short_name = 'l';
  RealBeta._double.header.long_name  = "RealBeta";
  RealBeta._double.header.help       = "beta (real)";
  RealBeta._double.required          = 0;
  RealBeta._double.min               = -DBL_MAX;
  RealBeta._double.max               = DBL_MAX;
  RealBeta._double.value             = 0.0;

  ImagBeta._double.header.type       = opt_double;
  ImagBeta._double.header.short_name = 'l';
  ImagBeta._double.header.long_name  = "ImagBeta";
  ImagBeta._double.header.help       = "beta (imaginary)";
  ImagBeta._double.required          = 0;
  ImagBeta._double.min               = -DBL_MAX;
  ImagBeta._double.max               = DBL_MAX;
  ImagBeta._double.value             = 0.0;
}

int matvec_show_help(void){
  matvec_fill_options_initialize();

  opt_show_option(FillA);
  opt_show_option(RealScaleA);
  opt_show_option(ImagScaleA);
  opt_show_option(FillX);
  opt_show_option(RealScaleX);
  opt_show_option(ImagScaleX);
  opt_show_option(FillY);
  opt_show_option(RealScaleY);
  opt_show_option(ImagScaleY);
  opt_show_option(RealAlpha);
  opt_show_option(ImagAlpha);
  opt_show_option(RealBeta);
  opt_show_option(ImagBeta);
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
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &RealScaleY);
  opt_eval_option(argc, argv, &ImagScaleY);
  opt_eval_option(argc, argv, &RealAlpha);
  opt_eval_option(argc, argv, &ImagAlpha);
  opt_eval_option(argc, argv, &RealBeta);
  opt_eval_option(argc, argv, &ImagBeta);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%s) (%s) (%s)", matvec_fill_name(argc, argv), util_mat_fill_names[FillA._named.value], util_vec_fill_names[FillX._named.value], util_vec_fill_names[FillY._named.value]);
  return name_buffer;
}

int matvec_test(int argc, char** argv, char Order, char TransA, int M, int N, int lda, int incX, int incY){
  int rc;

  matvec_fill_options_initialize();

  opt_eval_option(argc, argv, &FillA);
  opt_eval_option(argc, argv, &RealScaleA);
  opt_eval_option(argc, argv, &ImagScaleA);
  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &RealScaleX);
  opt_eval_option(argc, argv, &ImagScaleX);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &RealScaleY);
  opt_eval_option(argc, argv, &ImagScaleY);
  opt_eval_option(argc, argv, &RealAlpha);
  opt_eval_option(argc, argv, &ImagAlpha);
  opt_eval_option(argc, argv, &RealBeta);
  opt_eval_option(argc, argv, &ImagBeta);
  rc = matvec_fill_test(argc, argv, Order, TransA, M, N, RealAlpha._double.value, ImagAlpha._double.value, FillA._named.value, RealScaleA._double.value, ImagScaleA._double.value, lda, FillX._named.value, RealScaleX._double.value, ImagScaleX._double.value, incX, RealBeta._double.value, ImagBeta._double.value, FillY._named.value, RealScaleY._double.value, ImagScaleY._double.value, incY);
  return rc;
}
