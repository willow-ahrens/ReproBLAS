#include <float.h>
#include <stdio.h>

#include "test_opt.h"

#include "test_matmat_header.h"

int matmat_fill_show_help(void);
const char *matmat_fill_name(int argc, char** argv);
int matmat_fill_test(int argc, char** argv, char Order, char TransA, char TransB, int M, int N, int K, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillB, double RealScaleB, double ImagScaleB, int ldb, double RealBeta, double ImagBeta, int FillC, double RealScaleC, double ImagScaleC, int ldc);

static opt_option FillA;
static opt_option RealScaleA;
static opt_option ImagScaleA;
static opt_option FillB;
static opt_option RealScaleB;
static opt_option ImagScaleB;
static opt_option FillC;
static opt_option RealScaleC;
static opt_option ImagScaleC;
static opt_option RealAlpha;
static opt_option ImagAlpha;
static opt_option RealBeta;
static opt_option ImagBeta;

static void matmat_fill_options_initialize(void){
  FillA._named.header.type       = opt_named;
  FillA._named.header.short_name = '\0';
  FillA._named.header.long_name  = "FillA";
  FillA._named.header.help       = "A fill type";
  FillA._named.required          = 0;
  FillA._named.n_names           = (int)util_mat_fill_n_names;
  FillA._named.names             = (char**)util_mat_fill_names;
  FillA._named.descs             = (char**)util_mat_fill_descs;
  FillA._named.value             = 0;

  RealScaleA._double.header.type       = opt_double;
  RealScaleA._double.header.short_name = '\0';
  RealScaleA._double.header.long_name  = "RealScaleA";
  RealScaleA._double.header.help       = "A scale (real)";
  RealScaleA._double.required          = 0;
  RealScaleA._double.min               = -1 * DBL_MAX;
  RealScaleA._double.max               = DBL_MAX;
  RealScaleA._double.value             = 1.0;

  ImagScaleA._double.header.type       = opt_double;
  ImagScaleA._double.header.short_name = '\0';
  ImagScaleA._double.header.long_name  = "ImagScaleA";
  ImagScaleA._double.header.help       = "A scale (imaginary)";
  ImagScaleA._double.required          = 0;
  ImagScaleA._double.min               = -1 * DBL_MAX;
  ImagScaleA._double.max               = DBL_MAX;
  ImagScaleA._double.value             = 0.0;

  FillB._named.header.type       = opt_named;
  FillB._named.header.short_name = '\0';
  FillB._named.header.long_name  = "FillB";
  FillB._named.header.help       = "B fill type";
  FillB._named.required          = 0;
  FillB._named.n_names           = (int)util_mat_fill_n_names;
  FillB._named.names             = (char**)util_mat_fill_names;
  FillB._named.descs             = (char**)util_mat_fill_descs;
  FillB._named.value             = 0;

  RealScaleB._double.header.type       = opt_double;
  RealScaleB._double.header.short_name = '\0';
  RealScaleB._double.header.long_name  = "RealScaleB";
  RealScaleB._double.header.help       = "B scale (real)";
  RealScaleB._double.required          = 0;
  RealScaleB._double.min               = -1 * DBL_MAX;
  RealScaleB._double.max               = DBL_MAX;
  RealScaleB._double.value             = 1.0;

  ImagScaleB._double.header.type       = opt_double;
  ImagScaleB._double.header.short_name = '\0';
  ImagScaleB._double.header.long_name  = "ImagScaleB";
  ImagScaleB._double.header.help       = "B scale (imaginary)";
  ImagScaleB._double.required          = 0;
  ImagScaleB._double.min               = -1 * DBL_MAX;
  ImagScaleB._double.max               = DBL_MAX;
  ImagScaleB._double.value             = 0.0;

  FillC._named.header.type       = opt_named;
  FillC._named.header.short_name = '\0';
  FillC._named.header.long_name  = "FillC";
  FillC._named.header.help       = "C fill type";
  FillC._named.required          = 0;
  FillC._named.n_names           = (int)util_mat_fill_n_names;
  FillC._named.names             = (char**)util_mat_fill_names;
  FillC._named.descs             = (char**)util_mat_fill_descs;
  FillC._named.value             = 0;

  RealScaleC._double.header.type       = opt_double;
  RealScaleC._double.header.short_name = '\0';
  RealScaleC._double.header.long_name  = "RealScaleC";
  RealScaleC._double.header.help       = "C scale (real)";
  RealScaleC._double.required          = 0;
  RealScaleC._double.min               = -1 * DBL_MAX;
  RealScaleC._double.max               = DBL_MAX;
  RealScaleC._double.value             = 1.0;

  ImagScaleC._double.header.type       = opt_double;
  ImagScaleC._double.header.short_name = '\0';
  ImagScaleC._double.header.long_name  = "ImagScaleC";
  ImagScaleC._double.header.help       = "C scale (imaginary)";
  ImagScaleC._double.required          = 0;
  ImagScaleC._double.min               = -1 * DBL_MAX;
  ImagScaleC._double.max               = DBL_MAX;
  ImagScaleC._double.value             = 0.0;

  RealAlpha._double.header.type       = opt_double;
  RealAlpha._double.header.short_name = '\0';
  RealAlpha._double.header.long_name  = "RealAlpha";
  RealAlpha._double.header.help       = "alpha (real)";
  RealAlpha._double.required          = 0;
  RealAlpha._double.min               = -DBL_MAX;
  RealAlpha._double.max               = DBL_MAX;
  RealAlpha._double.value             = 1.0;

  ImagAlpha._double.header.type       = opt_double;
  ImagAlpha._double.header.short_name = '\0';
  ImagAlpha._double.header.long_name  = "ImagAlpha";
  ImagAlpha._double.header.help       = "alpha (imaginary)";
  ImagAlpha._double.required          = 0;
  ImagAlpha._double.min               = -DBL_MAX;
  ImagAlpha._double.max               = DBL_MAX;
  ImagAlpha._double.value             = 0.0;

  RealBeta._double.header.type       = opt_double;
  RealBeta._double.header.short_name = '\0';
  RealBeta._double.header.long_name  = "RealBeta";
  RealBeta._double.header.help       = "beta (real)";
  RealBeta._double.required          = 0;
  RealBeta._double.min               = -DBL_MAX;
  RealBeta._double.max               = DBL_MAX;
  RealBeta._double.value             = 0.0;

  ImagBeta._double.header.type       = opt_double;
  ImagBeta._double.header.short_name = '\0';
  ImagBeta._double.header.long_name  = "ImagBeta";
  ImagBeta._double.header.help       = "beta (imaginary)";
  ImagBeta._double.required          = 0;
  ImagBeta._double.min               = -DBL_MAX;
  ImagBeta._double.max               = DBL_MAX;
  ImagBeta._double.value             = 0.0;
}

int matmat_show_help(void){
  matmat_fill_options_initialize();

  opt_show_option(FillA);
  opt_show_option(RealScaleA);
  opt_show_option(ImagScaleA);
  opt_show_option(FillB);
  opt_show_option(RealScaleB);
  opt_show_option(ImagScaleB);
  opt_show_option(FillC);
  opt_show_option(RealScaleC);
  opt_show_option(ImagScaleC);
  opt_show_option(RealAlpha);
  opt_show_option(ImagAlpha);
  opt_show_option(RealBeta);
  opt_show_option(ImagBeta);
  return matmat_fill_show_help();
}

const char* matmat_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  matmat_fill_options_initialize();

  opt_eval_option(argc, argv, &FillA);
  opt_eval_option(argc, argv, &RealScaleA);
  opt_eval_option(argc, argv, &ImagScaleA);
  opt_eval_option(argc, argv, &FillB);
  opt_eval_option(argc, argv, &RealScaleB);
  opt_eval_option(argc, argv, &ImagScaleB);
  opt_eval_option(argc, argv, &FillC);
  opt_eval_option(argc, argv, &RealScaleC);
  opt_eval_option(argc, argv, &ImagScaleC);
  opt_eval_option(argc, argv, &RealAlpha);
  opt_eval_option(argc, argv, &ImagAlpha);
  opt_eval_option(argc, argv, &RealBeta);
  opt_eval_option(argc, argv, &ImagBeta);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s A(%s) B(%s) C(%s)", matmat_fill_name(argc, argv), util_mat_fill_names[FillA._named.value], util_mat_fill_names[FillB._named.value], util_mat_fill_names[FillC._named.value]);
  return name_buffer;
}

int matmat_test(int argc, char** argv, char Order, char TransA, char TransB, int M, int N, int K, int lda, int ldb, int ldc){
  int rc;

  matmat_fill_options_initialize();

  opt_eval_option(argc, argv, &FillA);
  opt_eval_option(argc, argv, &RealScaleA);
  opt_eval_option(argc, argv, &ImagScaleA);
  opt_eval_option(argc, argv, &FillB);
  opt_eval_option(argc, argv, &RealScaleB);
  opt_eval_option(argc, argv, &ImagScaleB);
  opt_eval_option(argc, argv, &FillC);
  opt_eval_option(argc, argv, &RealScaleC);
  opt_eval_option(argc, argv, &ImagScaleC);
  opt_eval_option(argc, argv, &RealAlpha);
  opt_eval_option(argc, argv, &ImagAlpha);
  opt_eval_option(argc, argv, &RealBeta);
  opt_eval_option(argc, argv, &ImagBeta);

  rc = matmat_fill_test(argc, argv, Order, TransA, TransB, M, N, K, RealAlpha._double.value, ImagAlpha._double.value, FillA._named.value, RealScaleA._double.value, ImagScaleA._double.value, lda, FillB._named.value, RealScaleB._double.value, ImagScaleB._double.value, ldb, RealBeta._double.value, ImagBeta._double.value, FillC._named.value, RealScaleC._double.value, ImagScaleC._double.value, ldc);
  return rc;
}
